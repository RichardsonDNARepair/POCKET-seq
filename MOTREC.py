"""
This script takes MACS3 callpeak output and calculates the closest feature from the input GFF3 file on the
corresponding chromosome. This is done by taking the midpoint of the MACS3 peak and matching it to the start
of the nearest feature in the GFF3 file. The value used to calculate the closest feature corresponds to the
strand the feature is located on.
"""

# imports
import argparse
import sys

import numpy as np
import pandas as pd


# read gff3 as df
def process_gff3_to_df(filepath):
    file = open(r"" + filepath, "r")
    txt = file.readlines()
    columns = [
        "chr",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]
    t = []
    for i in range(len(txt)):
        if "#" not in txt[i]:
            t.append(txt[i].split("\t"))

    file.close()
    df = pd.DataFrame(data=t, columns=columns)
    return df


# read macs3 output into pandas df
def process_macs3_to_df(fp):
    file = open(r"" + fp, "r")
    txt = file.readlines()
    columns = [
        "chr",
        "start",
        "end",
        "length",
        "abs_summit",
        "pileup",
        "-log10(pvalue)",
        "fold_enrichment",
        "-log10(qvalue)",
        "name",
    ]
    t = []
    for i in range(len(txt)):
        if "#" not in txt[i] and txt[i] != "\n":
            t.append(txt[i].split("\n")[0].split("\t"))
    file.close()
    df = pd.DataFrame(data=t, columns=columns)
    df = df.drop([0], axis="index")
    df[["start", "end", "fold_enrichment"]] = df[
        ["start", "end", "fold_enrichment"]
    ].apply(pd.to_numeric)
    df = df.reset_index(drop=True)
    return df


# get gene name from info and keep specified feature class only
def prepare_gff3(gff3, to_subset, ftc):
    if to_subset != "all":
        gff3 = gff3.loc[gff3["type"] == to_subset]
        gff3 = gff3.reset_index(drop=True)
    gene_name = []
    for i in range(len(gff3)):
        name = gff3.at[i, "attributes"].split(ftc)
        name = name[1].split(";")
        gene_name.append(name[0])
    gff3["gene"] = gene_name
    return gff3


# grab all the unique chromosomes in df
def create_chrom_list(df):
    e = pd.unique(df["chr"]).tolist()
    return e


# subset the data
def subset(df, chr):
    # subset data
    subset = df.query("chr == @chr")  # subset the data based on chromosome
    subset = subset.reset_index()  # reset old indices
    # get feature start
    subset = find_start_index(subset)
    return subset


# checks if feature is on top or bottom and gives an index to map to
def find_start_index(df):
    e = []
    for i in range(len(df)):
        if df.at[i, "strand"] == "+":
            e.append(float(df.at[i, "start"]))
        else:
            e.append(float(df.at[i, "end"]))
    df["start_point"] = e
    return df


# calc midpoint
def midpoint(df):
    e = []
    for i in range(len(df)):
        e.append((df.at[i, "end"] + df.at[i, "start"]) / 2)
    df["midpoint"] = e
    return df


# find nearest function
def find_nearest(arr, value):
    arr = np.asarray(arr)
    d = np.abs(arr - value)
    idx = (d).argmin()
    return idx


# heavy lifter
def map(to_map, subsets, bflag):
    # empty arrays
    closest_feature = []
    distance = []
    b = []

    # iterate through df
    for i in range(len(to_map)):
        # if chromosome in genome:
        if to_map.at[i, "chr"] in subsets[0]:
            # get index of correct chromosome in df of chromosome subsets
            chrom_index = subsets[0].index(to_map.at[i, "chr"])
            # get index of closest feature in chromosome df
            closest_feature_index = find_nearest(
                subsets[1][chrom_index][["start_point"]].to_numpy(),
                to_map.at[i, "midpoint"],
            )
            # take index and return gene symbol
            closest_feature.append(
                subsets[1][chrom_index].at[closest_feature_index, "gene"]
            )
            # calculate distance between feature start in chromosome and peak midpoint
            d = (
                np.abs(
                    subsets[1][chrom_index].at[closest_feature_index, "start_point"]
                    - to_map.at[i, "midpoint"]
                )
                / 1000
            )
            distance.append(d)

            # add a binary flag if distance is less than desired amount
            if d <= bflag:
                b.append(1)
            else:
                b.append(np.nan)

        # if chromosome not in genome: skip peak
        else:
            closest_feature.append(np.nan)
            distance.append(np.nan)
            b.append(np.nan)

    # add arrays to dataframe
    to_map["closest feature"] = closest_feature
    to_map["distance"] = distance
    to_map["within bflag"] = b

    # return dataframe
    return to_map


def add_phenotype(df, p_score, f_class):
    if ".csv" in p_score:
        p_df = pd.read_csv(p_score, delimiter=",")
    else:
        p_df = pd.read_csv(p_score, delimiter="\t")
    e = []
    for i in range(len(df)):
        gene = df.at[i, "closest feature"]
        p = p_df[p_df[f_class] == gene]
        if len(p) != 0:
            p = p.reset_index(drop=True)
            score = p.at[0, "epsilon"]
            e.append(score)
        else:
            e.append(np.nan)
    df["epsilon"] = e

    return df


# gene ontology mapping function
def go_mechanism(df, term, index_fp):
    # load data
    go_dataframe = pd.read_csv(index_fp, delimiter="\t")
    go_subset = go_dataframe.query("name == @term")  # subset GO based on desired term
    ctr = 0
    e = []
    for i in range(len(df)):  # iterate through df
        gene = df.at[i, "closest feature"]
        found = (
            go_subset["Object Symbol"].eq(gene).any()
        )  # find if feature present in subset
        # mark down results
        if found == True:
            e.append(1)
            ctr += 1
        else:
            e.append(np.nan)
    # save and return results
    df["go=" + term] = e
    return df, ctr


# argparse functionality
parser = argparse.ArgumentParser(description="process user in")
parser.add_argument(
    "--no_macs3",
    help="for files that are not MACS3 output. File should be a .xlsx with at minimum columns labeled 'chr', 'start', and 'end'",
    action="store_true",
)
parser.add_argument(
    "-g", "--gff3", help="input gff3 file (list of features by location)"
)
parser.add_argument("-p", "--peaks", help="input MACS3 peaks output.")
parser.add_argument(
    "-t",
    "--type",
    help="type of feature to subset in gff3 file. Default = transcript",
    default="transcript",
)
parser.add_argument(
    "-o",
    "--output",
    default="none",
    help="output file name. Not recommended if inputting more than one peaks file",
)
parser.add_argument(
    "-f",
    "--bflag",
    default=1,
    help="add a binary flag to all mapped peaks within a certain kb distance. Default = 1kb",
)
parser.add_argument(
    "-i",
    "--index",
    default="none",
    help="This represents the go index, a merge of the gene ontology annotation and ontology from geneontology.org. This can be created with the create_go_index.py script",
)
parser.add_argument(
    "--go",
    default="none",
    help="Gene Ontology term to search for. Will flag all matching results. Case sensitive. Default = no go search",
)
parser.add_argument(
    "--feature_class",
    default="symbol",
    help="type of classification for closest feature. Default = symbol. Options: 'symbol', 'ensg'",
)
parser.add_argument(
    "--phenotype",
    default="none",
    help="add filepath to phenotype data and annotate peaks. Tab separated text or comma separated are both acceptable. file should contain label named consistently with --feature class tag (symbol or ensg) and phenotype values column named epsilon",
)

args = parser.parse_args()
gff3_fp = args.gff3
peaks_fp = args.peaks
no_macs3 = args.no_macs3
bflag = float(args.bflag)
go = args.go
go_index = args.index
out = args.output
phenotype = args.phenotype
classification = args.feature_class
if classification == "symbol":
    ftc = "gene_name="  # feature class
elif classification == "ensg":
    ftc = "gene_id="
    go = "none"
    sys.stdout.write("Warning: ensg and go functionality non-compatible; go disabled")
    sys.stdout.write("\n")
else:
    sys.stdout.write("Warning: Invalid feature class; symbol used")
    sys.stdout.write("\n")
    ftc = "gene_name="

# read in all features in the human genome
sys.stdout.write("Loading in gff3...")
sys.stdout.write("\n")
features = process_gff3_to_df(gff3_fp)

sys.stdout.write("Preparing gff3...")
sys.stdout.write("\n")
features = prepare_gff3(features, args.type, ftc)


# subset features based on chromosome
sys.stdout.write("Subsetting Data...")
sys.stdout.write("\n")
chr_list = create_chrom_list(features)
df_list = []
for i in chr_list:
    df_list.append(subset(features, i))
subsets = [chr_list, df_list]


# grab file to map peaks back to
sys.stdout.write("Loading in peaks...")
sys.stdout.write("\n")

if no_macs3 == True:
    to_map = pd.read_excel(peaks_fp)
else:
    to_map = process_macs3_to_df(peaks_fp)

# prepare output name
if out != "none":
    o = out
else:
    o = peaks_fp + "_annotated"

# create the annotation
sys.stdout.write("Mapping...")
sys.stdout.write("\n")
to_map = midpoint(to_map)
to_map = map(to_map, subsets, bflag)
to_map = to_map.drop(columns=["midpoint"])

# annotate with phenotype score
if phenotype != "none":
    sys.stdout.write("Annotating with phenotype score...")
    sys.stdout.write("\n")
    to_map = add_phenotype(to_map, phenotype, classification)

# gene ontology
if (go != "none") and (go_index != "none"):
    sys.stdout.write("Conducting Gene Ontology Search...")
    sys.stdout.write("\n")
    to_map, counter = go_mechanism(to_map, go, go_index)
elif (go != "none") or (go_index != "none"):
    sys.stdout.write("Warning: either go term or go index unspecified")
    sys.stdout.write("\n")

to_map.to_csv(o + ".txt", sep="\t", mode="w", header=True, index=False)

# output file
sys.stdout.write("Finished...")
sys.stdout.write("\n")
