"""Rank promoter-proximal loci supported by distinct MACS3 replicate calls."""
import csv
import math
import heapq
from itertools import groupby
from bisect import bisect_left
from operator import attrgetter
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from common import open_text, table_writer


@dataclass(frozen=True)
class Call:
    chrom: str
    start: int
    end: int
    summit: int
    enrichment: float
    name: str
    replicate: int


def read_calls(path, replicate, qmax, contigs):
    calls = []
    seen_chroms = set()
    header = None
    with open_text(path) as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if fields[0] == "chr" and "abs_summit" in fields:
                header = fields
                continue
            try:
                if header is None:
                    if len(fields) != 10:
                        raise ValueError("expected a 10-column narrowPeak or MACS3 text peaks.xls with abs_summit")
                    chrom, a, b, name = fields[:4]
                    start, end = int(a), int(b)
                    summit = start + int(fields[9])
                    fe, logq = float(fields[6]), float(fields[8])
                else:
                    r = dict(zip(header, fields))
                    chrom, name = r["chr"], r["name"]
                    start, end, summit = int(r["start"])-1, int(r["end"]), int(r["abs_summit"])-1
                    fe, logq = float(r["fold_enrichment"]), float(r["-log10(qvalue)"])
                if not (0 <= start <= summit < end and fe > 0 and math.isfinite(fe) and math.isfinite(logq)):
                    raise ValueError("invalid interval, summit, enrichment or q-value")
            except (KeyError, ValueError) as exc:
                raise ValueError(f"{path}:{line_no}: {exc}") from exc
            seen_chroms.add(chrom)
            if (contigs is None or chrom in contigs) and logq+0.00001 >= -math.log10(qmax):
                calls.append(Call(chrom, start, end, summit, fe, name, replicate))
    return calls, seen_chroms


def components(calls):
    group = []
    right = -1
    for p in sorted(calls, key=lambda v:(v.chrom, v.start, v.end, v.name)):
        if group and (p.chrom != group[0].chrom or p.start >= right):
            yield group
            group = []
        right = max(right, p.end) if group else p.end
        group.append(p)
    if group:
        yield group


def shared_signals(group, n):
    """Strongest combination of replicate calls covering at least one common base."""
    calls = sorted((p for p in group if p.replicate < n), key=attrgetter("start"))
    active = [[] for _ in range(n)]
    best = None
    best_log = -math.inf
    serial = 0
    for start, entering in groupby(calls, key=attrgetter("start")):
        for call in entering:
            heap = active[call.replicate]
            heapq.heappush(heap, (-call.enrichment, call.name, call.end, call.summit, serial, call))
            serial += 1
        for heap in active:
            while heap and heap[0][-1].end <= start:
                heapq.heappop(heap)
        if all(active):
            strongest = [heap[0][-1] for heap in active]
            log_signal = sum(math.log(call.enrichment) for call in strongest)
            if log_signal > best_log:
                best = strongest
                best_log = log_signal
    return best


def fragment_counts(path, loci):
    index = defaultdict(list)
    for i, r in enumerate(loci):
        index[r["chrom"]].append((r["summit0"], i))
    ordered = {}
    for c, values in index.items():
        values.sort()
        ordered[c] = tuple(zip(*values))
    counts = [0]*len(loci)
    seen_chroms = set()
    with open_text(path) as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            f = line.split()
            if len(f) != 3:
                raise ValueError(f"{path}:{line_no}: expected a three-column fragment BED (chrom, start, end)")
            c, a, b = f[0], int(f[1]), int(f[2])
            if not 0 <= a < b:
                raise ValueError(f"{path}:{line_no}: invalid fragment")
            seen_chroms.add(c)
            if c in ordered:
                pos, idx = ordered[c]
                for j in range(bisect_left(pos, a-250), bisect_left(pos, b+250)):
                    counts[idx[j]] += 1
    if loci and not (seen_chroms & ordered.keys()):
        raise ValueError(f"{path}: no fragment contigs match the replicated loci")
    return counts


def load_terms(path):
    by_id, by_name = defaultdict(set), defaultdict(set)
    with open_text(path) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = set(reader.fieldnames or ())
        if "term_id" not in fields or not {"gene_id", "gene_name"} & fields:
            raise ValueError(f"{path}: annotations require term_id and gene_id or gene_name; term_name is optional")
        for i, r in enumerate(reader, 2):
            if None in r or None in r.values() or not r["term_id"]:
                raise ValueError(f"{path}:{i}: invalid annotation row")
            key = r.get("gene_id") or r.get("gene_name")
            if not key:
                raise ValueError(f"{path}:{i}: missing gene key")
            dest = by_id if r.get("gene_id") else by_name
            dest[key].add((r["term_id"], r.get("term_name", "") or r["term_id"]))
    if not by_id and not by_name:
        raise ValueError(f"{path}: empty pathway annotation")
    return by_id, by_name


def distinct_files(paths):
    identities = [(Path(p).stat().st_dev, Path(p).stat().st_ino) for p in paths]
    if len(set(identities)) != len(paths):
        raise ValueError("Replicate files must be distinct; the same file cannot supply two libraries")


def run(args, out, promoters, contigs):
    n = len(args.replicates)
    if n < 2 or not (0 < args.q <= 1):
        raise ValueError("Supply ≥2 independent replicate peak files and 0 < q ≤ 1")
    paths = [*args.replicates, *([args.pooled] if args.pooled else [])]
    distinct_files(paths)
    if args.fragments is not None or args.effective_sizes is not None:
        if len(args.fragments or []) != n or len(args.effective_sizes or []) != n:
            raise ValueError("Supply one fragment BED and one effective library size per replicate, in the same order")
        if not all(math.isfinite(v) and v > 0 for v in args.effective_sizes):
            raise ValueError("Effective library sizes must be positive finite values")
        distinct_files(args.fragments)
    terms_id, terms_name = load_terms(args.pathways)
    target = None
    if args.target:
        chrom, a, b = args.target
        target = (chrom, int(a), int(b))
        if not 0 <= target[1] < target[2] or (contigs is not None and chrom not in contigs):
            raise ValueError("Invalid intended target interval or excluded chromosome")
    calls = []
    per_file = []
    observable = set()
    for rep, path in enumerate(paths):
        part, seen = read_calls(path, rep, args.q, contigs)
        calls.extend(part)
        per_file.append(len(part))
        observable.update(seen)
        if seen and not seen.intersection(contigs if contigs is not None else promoters.genes):
            raise ValueError(f"{path}: no contig names overlap the selected reference")
    if observable and not observable.intersection(promoters.genes):
        raise ValueError("No shared contigs between peak files and TSS annotation")
    loci = []
    for group in components(calls):
        support = shared_signals(group, n)
        if support is None:
            continue
        core_start = max(p.start for p in support)
        core_end = min(p.end for p in support)
        core_calls = (p for p in group if p.start < core_end and p.end > core_start)
        best = min(core_calls, key=lambda p: (-p.enrichment, p.name, p.summit))
        a, b = group[0].start, max(p.end for p in group)
        intended = target and target[0] == best.chrom and target[1] < b and target[2] > a
        row = dict(guide=args.guide, locus_id=f"{args.guide}|{best.chrom}:{a}-{b}", chrom=best.chrom, start0=a, end0=b,
                 summit0=best.summit, summit_peak=best.name,
                 core_start0=core_start, core_end0=core_end,
                 target_status="unknown" if target is None else ("intended" if intended else "off_target"))
        for i, p in enumerate(support, 1):
            row[f"rep{i}_peak"] = p.name
            row[f"rep{i}_FE"] = p.enrichment
        loci.append(row)
    for rep in range(1, n+1):
        if args.fragments:
            counts = fragment_counts(args.fragments[rep-1], loci)
            for r, v in zip(loci, counts):
                r[f"rep{rep}_fragments"] = v
                r[f"rep{rep}_signal"] = v*1e6/args.effective_sizes[rep-1]
        else:
            for r in loci:
                r[f"rep{rep}_signal"] = r[f"rep{rep}_FE"]
    unit = "nCPM" if args.fragments else "fold_enrichment"
    for r in loci:
        vals = [r[f"rep{i}_signal"] for i in range(1, n+1)]
        r["signal_geomean"] = math.exp(sum(map(math.log, vals))/n) if all(vals) else 0.0
        r["signal_unit"] = unit
    loci.sort(key=lambda r:(-r["signal_geomean"], r["chrom"], r["start0"]))
    fields = ["guide", "locus_id", "chrom", "start0", "end0", "summit0", "summit_peak", "core_start0", "core_end0",
            "target_status", "promoter_rank", "signal_geomean", "signal_unit"]
    for i in range(1, n+1):
        fields.extend([f"rep{i}_peak", f"rep{i}_FE", f"rep{i}_signal"])
        if args.fragments:
            fields.append(f"rep{i}_fragments")
    fields += ["gene_id", "gene_name", "transcript", "TSS0", "gene_strand", "signed_distance_bp", "term_ids", "term_names"]
    rank = links = annotated = on_links = 0
    with table_writer(out / "promoter_peaks.tsv", fields) as writer, table_writer(out / "intended_promoters.tsv", fields) as controls:
        for locus in loci:
            genes = promoters.near(locus["chrom"], locus["summit0"], args.window)
            if not genes:
                continue
            intended = locus["target_status"] == "intended"
            if not intended:
                rank += 1
            for gene in genes:
                terms = sorted(terms_id.get(gene.gene_id, set()) | terms_name.get(gene.gene_name, set()))
                row = dict(locus, promoter_rank="NA" if intended else rank,
                           gene_id=gene.gene_id, gene_name=gene.gene_name, transcript=gene.transcript,
                           TSS0=gene.position, gene_strand=gene.strand,
                           signed_distance_bp=gene.distance(locus["summit0"]),
                           term_ids=";".join(term for term, _ in terms),
                           term_names=";".join(label for _, label in terms))
                if intended:
                    controls.writerow(row)
                    on_links += 1
                else:
                    writer.writerow(row)
                    links += 1
                    annotated += bool(terms)
    return dict(q_passing_calls_per_input=per_file, replicated_loci=len(loci), ranked_promoter_loci=rank,
                peak_gene_TSS_links=links, links_with_pathway_terms=annotated,
                intended_promoter_links=on_links, signal_unit=unit)
