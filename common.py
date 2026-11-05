"""Shared coordinate and table handling."""
import csv
import gzip
from bisect import bisect_left, bisect_right
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"


def open_text(path, mode="rt"):
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    return opener(path, mode, encoding="utf-8-sig" if "r" in mode else "utf-8", newline="")


@contextmanager
def table_writer(path, fields):
    with open_text(path, "wt") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        yield writer


def table_rows(path, required):
    with open_text(path) as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        if len(set(fields)) != len(fields) or not set(required).issubset(fields):
            raise ValueError(f"{path}: expected unique headers including {', '.join(required)}")
        for line, row in enumerate(reader, 2):
            if None in row or None in row.values():
                raise ValueError(f"{path}:{line}: inconsistent number of columns")
            yield row


@dataclass(frozen=True, order=True)
class TSS:
    position: int
    gene_id: str
    gene_name: str
    strand: str
    transcript: str

    def distance(self, site):
        return (site - self.position) * (1 if self.strand == "+" else -1)


class Promoters:
    def __init__(self, path):
        rows = defaultdict(set)
        for row in table_rows(path, ["chrom", "TSS0", "strand", "gene_id", "gene_name"]):
            pos = int(row["TSS0"])
            if pos < 0 or row["strand"] not in {"+", "-"} or not all(row[k] for k in ("chrom", "gene_id", "gene_name")):
                raise ValueError(f"{path}: invalid TSS row {row}")
            rows[row["chrom"]].add(TSS(pos, row["gene_id"], row["gene_name"], row["strand"], row.get("refseq_transcript", "")))
        if not rows:
            raise ValueError(f"{path}: empty TSS annotation")
        self.genes = {chrom: sorted(genes) for chrom, genes in rows.items()}
        self.positions = {chrom: [t.position for t in genes] for chrom, genes in self.genes.items()}

    def near(self, chrom, site, window):
        positions = self.positions.get(chrom, ())
        start = bisect_left(positions, site - window)
        end = bisect_right(positions, site + window)
        return self.genes.get(chrom, ())[start:end]


def read_contigs(path, assembly):
    if path is None:
        return tuple(f"chr{n}" for n in [*range(1, 23), "X", "Y"]) if assembly == "hg19" else None
    if path == "all":
        return None
    with open_text(path) as stream:
        names = [line.split()[0] for line in stream if line.strip() and not line.startswith("#")]
    if not names or len(set(names)) != len(names):
        raise ValueError(f"{path}: expected unique contig names")
    return tuple(names)
