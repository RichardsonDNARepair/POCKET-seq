"""Frozen sequence energy and height model used in the manuscript."""
import csv
import json

from common import DATA


class HeightModel:
    def __init__(self):
        with (DATA / "algorithm_parameters.tsv").open() as stream:
            steps = {
                row["term"]: round(1_000_000 * float(row["value"]))
                for row in csv.DictReader(stream, delimiter="\t")
                if row["model"] == "NN37_table"
            }
        pairs = [a + b for a in "ACGT" for b in "ACGT"]
        self.dna = {pair: steps["DNA_DNA_" + pair] for pair in pairs}
        self.hybrid = {pair: steps["RNA_DNA_" + pair] for pair in pairs}
        model = json.loads((DATA / "frozen_model.json").read_text())
        self.weights = model["coefficients"]
        self.intercept = model["intercept"]
        self.model_id = model["model_id"]

    def e12(self, guide, site):
        genomic, spacer = site[8:20], guide[-12:]
        total = -2_000_000 * (genomic[-1] != spacer[-1])
        for i in range(11):
            pair = genomic[i:i+2]
            total += self.dna[pair] + 500_000
            if pair == spacer[i:i+2]:
                total -= self.hybrid[pair]
        cents = (abs(total) + 5_000) // 10_000
        return (-cents if total < 0 else cents) / 100

    def eta(self, energy, gc, post):
        w = self.weights
        return (self.intercept + w["E12"] * energy + w["GC250"] * gc
                + w["E12_hinge6"] * max(0, energy - 6) + w["post1_G"] * (post == "G"))
