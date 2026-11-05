#!/usr/bin/env python3
"""Convert matching GO GAF and OBO releases into a biological-process lookup."""
import argparse
from pathlib import Path

from common import open_text, table_writer


def term_names(path):
    names = {}
    current = {}
    with open_text(path) as stream:
        for line in stream:
            line = line.strip()
            if line.startswith("[") or not line:
                if current.get("namespace") == "biological_process" and current.get("is_obsolete") != "true":
                    names[current["id"]] = current["name"]
                current = {}
            elif ": " in line:
                field, value = line.split(": ", 1)
                if field in {"id", "name", "namespace", "is_obsolete"}:
                    current[field] = value
    if current.get("namespace") == "biological_process" and current.get("is_obsolete") != "true":
        names[current["id"]] = current["name"]
    return names


def make_mapping(gaf, obo, out, taxon):
    names = term_names(obo)
    seen = set()
    with open_text(gaf) as stream:
        for i, line in enumerate(stream, 1):
            if line.startswith("!") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 15:
                raise ValueError(f"{gaf}:{i}: incomplete GAF record")
            symbol, relation, term = fields[2:5]
            organism = fields[12].split("|")[0]
            if (symbol.strip() and organism == f"taxon:{taxon}" and fields[8] == "P" and "NOT" not in relation.split("|")
                    and term in names and term != "GO:0008150"):
                seen.add((symbol, term))
    if not seen:
        raise ValueError("No biological-process associations in the supplied files for this taxon")
    if out.exists():
        raise ValueError(f"Output exists: {out}")
    with table_writer(out, ["gene_name", "term_id", "term_name"]) as writer:
        for gene, term in sorted(seen):
            writer.writerow(dict(gene_name=gene, term_id=term, term_name=names[term]))
    return len(seen)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gaf", type=Path, required=True)
    parser.add_argument("--obo", type=Path, required=True)
    parser.add_argument("--taxon", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    print(make_mapping(args.gaf, args.obo, args.out, args.taxon), "associations")
