"""Reader and hg19 target verification for the original CRISPRi v2 Top5 library."""
import re
from collections import defaultdict

COLUMNS = ["guide", "sublibrary", "subset", "original_column4", "target_gene",
           "original_column6", "oligo", "original_activity_score", "guide_sequence",
           "promoter_group"]
TARGET_FIELDS = ["target_chrom", "target_start0", "target_end0", "intended_mapping"]
LEFT = "CCACCTTGTTG"
RIGHT = "GTTTAAGAGCTAAGCTG"
RC = str.maketrans("ACGT", "TGCA")
NAME = re.compile(r"(.+)_([+-])_(\d+)\.23-(.+)")


def looks_like_top5(fields):
    return len(fields) == len(COLUMNS) and fields[2] == "Top5"


def read_row(fields, line):
    if not looks_like_top5(fields):
        raise ValueError(f"Top5 row {line}: expected ten original columns with Top5 in column 3")
    row = dict(zip(COLUMNS, fields))
    sequence = row["guide_sequence"].upper()
    if not re.fullmatch("G[ACGT]{19}", sequence):
        raise ValueError(f"Top5 row {line}: column 9 must be the 5′ G followed by 19 targeting bases")
    if row["oligo"].upper() != LEFT + sequence + RIGHT:
        raise ValueError(f"Top5 row {line}: column 9 disagrees with the spacer between the column-7 cloning arms")
    row.update(dict.fromkeys(TARGET_FIELDS, ""))
    return row


def resolve_targets(rows, fasta, contigs, promoters, assembly):
    by_gene = defaultdict(set)
    for chrom, entries in promoters.genes.items():
        if chrom in contigs:
            for gene in entries:
                by_gene[gene.gene_name].add(chrom)
    for row in rows:
        if row["target_gene"] == "negative_control":
            row["intended_mapping"] = "non_targeting_control"
            continue
        row["intended_mapping"] = "unresolved"
        name = NAME.fullmatch(row["guide"])
        if assembly != "hg19" or not name or name[1] != row["target_gene"]:
            continue
        plus = name[2] == "-"
        start = int(name[3]) - (22 if plus else 0)
        if start < 0:
            continue
        hint = by_gene[row["target_gene"]]
        suffix = row["guide_sequence"].upper()[1:]
        for chromosomes in (hint, contigs - hint):
            hits = []
            for chrom in chromosomes:
                if start + 23 > fasta.get_reference_length(chrom):
                    continue
                site = fasta.fetch(chrom, start, start + 23).upper()
                if not plus:
                    site = site.translate(RC)[::-1]
                if site[1:20] == suffix and site.endswith("GG") and re.fullmatch("[ACGT]{23}", site):
                    hits.append(chrom)
            if hits:
                if len(hits) == 1:
                    row.update(target_chrom=hits[0], target_start0=start, target_end0=start + 23,
                               intended_mapping="verified_hg19_19nt_NGG")
                break
