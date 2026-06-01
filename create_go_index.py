import argparse
import sys

import pandas as pd

parser = argparse.ArgumentParser(description="process user in")
parser.add_argument(
    "-a", "--annotation", help="downloaded from geneontology.org in .gaf format"
)
parser.add_argument(
    "-o", "--ontology", help="downloaded from geneontology.org in .obo format"
)

args = parser.parse_args()
annotation_fp = args.annotation
ontology_fp = args.ontology


def load_annotation(annotation_fp):
    data = open(r"" + annotation_fp, "r")
    annotation_raw = data.readlines()
    annotation_columns = [
        "Database",
        "Object ID",
        "Object Symbol",
        "Qualifier",
        "GO ID",
        "Reference",
        "Evidence Code",
        "With/From",
        "Aspect",
    ]
    t = []
    data = []
    c = 0
    for i in range(len(annotation_raw)):
        z = []
        if "!" not in annotation_raw[i]:
            t.append(annotation_raw[i].split("\t"))
            for j in range(9):
                z.append(t[c][j])
            data.append(z)
            c += 1
    annotation = pd.DataFrame(data=data, columns=annotation_columns)
    annotation = annotation[["Object Symbol", "GO ID"]]
    return annotation


def load_ontology(ontology_fp):
    data = open(r"" + ontology_fp, "r")
    ontology_raw = data.readlines()
    ontology_columns = ["GO ID", "name", "namespace"]
    t = []
    for i in range(len(ontology_raw)):
        if "[Term]" in ontology_raw[i]:
            temp_id = ontology_raw[i + 1].split("\n")[0].split("id: ")[1]
            temp_name = ontology_raw[i + 2].split("\n")[0].split("name: ")[1]
            temp_namespace = ontology_raw[i + 3].split("\n")[0].split("namespace: ")[1]
            t.append([temp_id, temp_name, temp_namespace])
    ontology = pd.DataFrame(data=t, columns=ontology_columns)
    return ontology


annotation = load_annotation(annotation_fp)
ontology = load_ontology(ontology_fp)
index = annotation.merge(ontology, how="left", on="GO ID")
index.to_csv("go_index.txt", sep="\t", mode="w", header=True, index=False)
