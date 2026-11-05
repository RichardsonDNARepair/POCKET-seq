"""Enumerate exact-8 or one-mismatch-8 plus exact-p9–p12 sites on both strands."""
import csv
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from itertools import chain

import numpy as np
import pysam

from common import DATA, open_text, table_writer
from energy import HeightModel
import weissman

RC = str.maketrans("ACGT", "TGCA")
KEYS = ["target_chrom", "target_start0", "target_end0"]
BLOCK = 1_000_000
METRICS = ["sites", "scored_sites", "mean_eta", "sum_eta", "max_eta"]
SUMMARY = ["seed_8nt", "extension_p9_p12", "target_status", "excluded_intended_sites",
           "genome_perfect_matches", *["genome_" + k for k in METRICS],
           *["promoter_" + k for k in METRICS]]
SITE_FIELDS = ["guide", "chrom", "start0", "end0", "strand", "sequence_5to3", "genomic_seed_8nt",
               "seed_mid0", "seed_mismatches", "branch", "full_spacer_match", "E12", "GC250",
               "callable_flank_bases", "post_PAM_base", "post1_G", "eta", "target_status",
               "promoter_gene_ids", "promoter_gene_names"]


@dataclass
class Magnitudes:
    sites: int = 0
    scored: int = 0
    total: float = 0.0
    maximum: float = float("-inf")

    def add(self, value):
        self.sites += 1
        if value is not None:
            self.scored += 1
            self.total += value
            self.maximum = max(self.maximum, value)

    def values(self, prefix):
        numbers = (self.sites, self.scored, self.total/self.scored if self.scored else "NA",
                   self.total if self.scored or not self.sites else "NA",
                   self.maximum if self.scored else "NA")
        return dict(zip((prefix + "_" + k for k in METRICS), numbers))


def spacer(row):
    return row["guide_sequence"].strip().upper().replace("U", "T")


def read_guides(path):
    with open_text(path) as stream:
        lines = (line for line in stream if line.strip() and not line.startswith("#"))
        first = next(lines, "").rstrip("\r\n")
        if not first:
            raise ValueError(f"{path}: empty guide list")
        headings = next(csv.reader([first], delimiter="\t"))
        top5 = weissman.looks_like_top5(headings)
        if top5:
            records = csv.reader(chain([first], lines), delimiter="\t")
            rows = [weissman.read_row(fields, i) for i, fields in enumerate(records, 1)]
            headings = weissman.COLUMNS + weissman.TARGET_FIELDS
        elif "guide_sequence" in headings:
            if len(set(headings)) != len(headings) or set(headings) & set(SUMMARY):
                raise ValueError(f"{path}: duplicated or reserved column name")
            rows = list(csv.DictReader(lines, headings, delimiter="\t"))
        elif len(headings) == 1:
            rows = [{"guide_sequence": value.strip()} for value in chain([first], lines)]
            headings = ["guide_sequence"]
        else:
            raise ValueError(f"{path}: unknown table format; supply a guide_sequence header or original Weissman v2 Top5 rows")
    if "guide" not in headings:
        headings = ["guide", *headings]
        for i, row in enumerate(rows, 1):
            row["guide"] = f"guide_{i}"
    seen = set()
    for row in rows:
        if None in row or None in row.values():
            raise ValueError(f"{path}: inconsistent number of columns")
        name = row["guide"]
        if not name or name in seen:
            raise ValueError(f"Empty or duplicated guide ID: {name!r}")
        seen.add(name)
        if not re.fullmatch("[ACGT]{19,20}", spacer(row)):
            raise ValueError(f"{name}: supply a 19- or 20-nt spacer without PAM, in guide orientation")
        supplied = [row.get(k, "") not in {"", ".", "NA"} for k in KEYS]
        if any(supplied) and not all(supplied):
            raise ValueError(f"{name}: supply all three intended target fields or none")
        if all(supplied):
            start, end = int(row[KEYS[1]]), int(row[KEYS[2]])
            if start < 0 or end-start != 23:
                raise ValueError(f"{name}: target must be a 0-based half-open, 23-bp protospacer+PAM interval")
    if not rows:
        raise ValueError(f"{path}: no guide rows")
    return headings, rows, top5


def code(seq):
    key = 0
    for base in seq:
        key = 4*key + "ACGT".index(base)
    return key


class SeedIndex:
    def __init__(self, sequences):
        self.exact = defaultdict(list)
        self.extended = defaultdict(list)
        for seq in sorted(sequences):
            self.exact[code(seq[-8:])].append(seq)
            for i in range(8):
                for base in "ACGT":
                    if base != seq[-8:][i]:
                        variant = seq[-8:][:i] + base + seq[-8:][i+1:]
                        self.extended[code(seq[-12:-8] + variant)].append(seq)
        self.mask8 = np.zeros(4**8, dtype=bool)
        self.mask12 = np.zeros(4**12, dtype=bool)
        self.mask8[list(self.exact)] = True
        self.mask12[list(self.extended)] = True


BASE_CODE = np.full(256, 4, dtype=np.uint8)
BASE_CODE[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4)


def matches(sequence, first, stop, index):
    """Yield oriented 23-mer offsets whose seed passes either branch."""
    dna = BASE_CODE[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
    starts = np.flatnonzero((dna[first+21:stop+21] == 2) & (dna[first+22:stop+22] == 2)) + first
    if starts.size == 0:
        return
    key = np.zeros(starts.size, dtype=np.int32)
    for offset in range(8, 20):
        key = (key << 2) | (dna[starts+offset] & 3)
    key8 = key & 65535
    chosen = np.flatnonzero(index.mask8[key8] | index.mask12[key])
    for j in chosen:
        p = int(starts[j])
        site = sequence[p:p+23]
        if re.fullmatch("[ACGT]{23}", site):
            for seq in index.exact.get(key8[j], ()):
                yield p, site, seq, 0
            for seq in index.extended.get(key[j], ()):
                yield p, site, seq, 1


def scan_chrom(fa, chrom, index):
    length = fa.get_reference_length(chrom)
    for a in range(0, length-22, BLOCK):
        b = min(a+BLOCK, length-22)
        lo, hi = max(0, a-250), min(length, b+22+250)
        forward = fa.fetch(chrom, lo, hi).upper()
        for strand in ("+", "-"):
            if strand == "+":
                dna = forward
                first, stop = a-lo, b-lo
            else:
                dna = forward.translate(RC)[::-1]
                first, stop = hi-b-22, hi-a-22
            for pos, site, seq, mismatches in matches(dna, first, stop, index):
                start = lo+pos if strand == "+" else hi-pos-23
                flanks = dna[max(0, pos-250):pos] + dna[pos+23:pos+273]
                n = sum(flanks.count(base) for base in "ACGT")
                gc = (flanks.count("C")+flanks.count("G"))/n if n else None
                post = dna[pos+23:pos+24] or "N"
                yield start, strand, site, seq, mismatches, gc, n, post


def target_map(rows, fa, chosen, top5=False):
    targets = {}
    for r in rows:
        if r.get("target_chrom", "") in {"", ".", "NA"}:
            continue
        chrom, start, end = r["target_chrom"], int(r["target_start0"]), int(r["target_end0"])
        if chrom not in chosen or end > fa.get_reference_length(chrom):
            raise ValueError(f"{r['guide']}: intended interval outside requested reference contigs")
        s = fa.fetch(chrom, start, end).upper()
        seq = spacer(r)[1:] if top5 else spacer(r)
        strands = [st for st, z in [("+", s), ("-", s.translate(RC)[::-1])]
                   if z.endswith("GG") and z[:20].endswith(seq) and re.fullmatch("[ACGT]{23}", z)]
        if len(strands) != 1:
            raise ValueError(f"{r['guide']}: intended 23-mer does not uniquely match this spacer + NGG on the supplied FASTA")
        targets[r["guide"]] = (chrom, start, strands[0])
    return targets


def run(args, out, promoters, contigs):
    headings, guides, top5 = read_guides(args.guides)
    if args.emit_sites is None:
        args.emit_sites = "none" if top5 else "promoter"
    by_seq = defaultdict(list)
    for r in guides:
        by_seq[spacer(r)].append(r["guide"])
    model = HeightModel()
    genome, promoter = defaultdict(Magnitudes), defaultdict(Magnitudes)
    perfect, excluded = defaultdict(int), defaultdict(int)
    if not args.fasta.is_file() or not args.fasta.with_name(args.fasta.name+".fai").is_file():
        raise ValueError("Provide a FASTA with its .fai index (create with pysam.faidx or samtools faidx)")
    with pysam.FastaFile(str(args.fasta)) as fa, table_writer(out/"sites.tsv.gz", SITE_FIELDS) as writer:
        chosen = tuple(fa.references) if contigs is None else contigs
        if set(chosen)-set(fa.references):
            raise ValueError("Requested contigs absent from FASTA: "+", ".join(sorted(set(chosen)-set(fa.references))[:8]))
        if args.assembly == "hg19":
            sizes = json.loads((DATA/"hg19_primary_lengths.json").read_text())
            if any(c in sizes and fa.get_reference_length(c) != sizes[c] for c in chosen):
                raise ValueError("FASTA chromosome lengths disagree with hg19")
        missing = set(promoters.positions)-set(fa.references)
        if missing and contigs is None:
            raise ValueError("TSS contigs absent from FASTA: "+", ".join(sorted(missing)[:8]))
        common = set(chosen)&promoters.positions.keys()
        if not common or any(promoters.positions[c][-1] >= fa.get_reference_length(c) for c in common):
            raise ValueError("TSS coordinates or contig names disagree with the selected genome")
        if top5:
            weissman.resolve_targets(guides, fa, set(chosen), promoters, args.assembly)
        targets = target_map(guides, fa, set(chosen), top5)
        bulk = args.emit_sites == "none" and len(by_seq) >= 1000
        if bulk:
            from library_summary import summarize
            genome, promoter, perfect, excluded = summarize(fa, chosen, guides, targets, promoters, args.window, model)
        else:
            idx = SeedIndex(by_seq)
        for chrom in (() if bulk else chosen):
            for start, st, site, seq, mm, gc, n, post in scan_chrom(fa, chrom, idx):
                center = start+(15.5 if st == "+" else 6.5)
                promoters_here = promoters.near(chrom, center, args.window)
                e = model.e12(seq, site)
                eta = model.eta(e, gc, post) if gc is not None else None
                full_match = site[:20].endswith(seq)
                for name in by_seq[seq]:
                    intended = targets.get(name) == (chrom, start, st)
                    perfect[name] += int(full_match)
                    if intended:
                        excluded[name] += 1
                    else:
                        genome[name].add(eta)
                        if promoters_here:
                            promoter[name].add(eta)
                    if args.emit_sites == "none" or (args.emit_sites == "promoter" and not promoters_here):
                        continue
                    genes = sorted({(t.gene_id, t.gene_name) for t in promoters_here})
                    writer.writerow(dict(guide=name, chrom=chrom, start0=start, end0=start+23, strand=st,
                        sequence_5to3=site, genomic_seed_8nt=site[12:20], seed_mid0=center, seed_mismatches=mm,
                        branch="exact8" if mm == 0 else "one_mm_plus_exact_p9_p12", full_spacer_match=int(full_match),
                        E12=e, GC250=gc if n else "NA", callable_flank_bases=n, post_PAM_base=post,
                        post1_G=int(post == "G"), eta=eta if n else "NA",
                        target_status="intended" if intended else ("off_target" if name in targets else "unknown"),
                        promoter_gene_ids=";".join(gid for gid, _ in genes),
                        promoter_gene_names=";".join(symbol for _, symbol in genes)))
    with table_writer(out/"guide_summary.tsv", headings+SUMMARY) as writer:
        for r in guides:
            name = r["guide"]
            seq = spacer(r)
            if name in targets and excluded[name] != 1:
                raise ValueError(f"{name}: expected exactly one verified intended site")
            writer.writerow(dict(r, seed_8nt=seq[-8:], extension_p9_p12=seq[-12:-8],
                target_status="specified" if name in targets else "unknown", excluded_intended_sites=excluded[name],
                genome_perfect_matches=perfect[name], **genome[name].values("genome"), **promoter[name].values("promoter")))
    return dict(model=model.model_id, contigs_scanned=chosen, guide_rows=len(guides), unique_spacers=len(by_seq),
                original_weissman_top5=top5, exact_word_summary=bulk, emitted_site_scope=args.emit_sites)
