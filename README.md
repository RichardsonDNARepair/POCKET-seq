# POCKET-seq analysis pipeline

MOTREC (Mapping Off-Targets in Transcription Regulation Experiments with Callpeak Analysis) is a Python command-line tool for annotating ChIP-seq peak calls from POCKET-seq experiments with their nearest genomic features, phenotype scores, and Gene Ontology terms. It was developed to identify on- and off-target binding events in CRISPR transcriptional regulation experiments (CRISPRi).

## How It Works

MOTREC takes MACS3 peak calls and a genome annotation (GFF3) as input. For each peak, it:

1. Calculates the peak midpoint
2. Identifies the nearest annotated genomic feature on the correct strand
3. Reports the closest feature name and distance in kb
4. Optionally flags peaks within a user-defined distance threshold
5. Optionally annotates peaks with phenotype scores from a prior screen (e.g. CRISPRi epsilon scores)
6. Optionally flags peaks whose nearest feature belongs to a specified Gene Ontology term

This allows researchers to cross-reference ChIP-seq binding locations with functional genomic 
data to separate biologically meaningful off-target events from background noise.

## Scripts

| Script | Description |
|---|---|
| `MOTREC.py` | Main pipeline — peak annotation, phenotype integration, GO flagging |
| `create_go_index.py` | Preprocesses GO annotation + ontology files into a merged index for use with `--go` |

## Usage

### Basic — annotate peaks with nearest genomic feature
```bash
python MOTREC.py -g genome.gff3 -p peaks.xls
```

### Full — with GO annotation, phenotype scores, and custom distance flag
```bash
python MOTREC.py \
  -g gencode.hg38.gff3 \
  -p MACS3_peaks.xls \
  -o MOTREC_output \
  -i go_index.txt \
  --go "DNA repair" \
  --phenotype screen_epsilon_scores.csv \
  -f 5 \
  --type transcript \
  --feature_class symbol
```

### Creating the GO index (required for `--go`)
```bash
python create_go_index.py \
  -a goa_human.gaf \
  -o go-basic.obo
```
GO annotation (`.gaf`) and ontology (`.obo`) files can be downloaded from [geneontology.org](https://geneontology.org/docs/download-go-annotations/).

## Arguments

### MOTREC.py

| Flag | Description | Default |
|---|---|---|
| `-g, --gff3` | Genome annotation in GFF3 format (e.g. from [GENCODE](https://www.gencodegenes.org/)) | required |
| `-p, --peaks` | MACS3 callpeak output (`.xls`) | required |
| `-o, --output` | Output filename (without extension) | `<peaks_input>_annotated` |
| `-t, --type` | Feature type to use from GFF3 (`transcript`, `gene`, `CDS`, etc.) | `transcript` |
| `-f, --bflag` | Distance threshold in kb for binary proximity flag | `1` |
| `--feature_class` | Feature identifier type: `symbol` or `ensg` | `symbol` |
| `--phenotype` | Path to phenotype score file (tab- or comma-separated, must contain `epsilon` column) | none |
| `-i, --index` | Path to GO index created by `create_go_index.py` | none |
| `--go` | GO term to flag (case-sensitive, requires `--index`) | none |
| `--no_macs3` | Use for non-MACS3 peak files (expects `.xlsx` with `chr`, `start`, `end` columns) | false |

### create_go_index.py

| Flag | Description |
|---|---|
| `-a, --annotation` | GO annotation file in `.gaf` format |
| `-o, --ontology` | GO ontology file in `.obo` format |


## Output

A tab-separated `.txt` file containing all original peak columns plus:

| Column | Description |
|---|---|
| `closest feature` | Gene symbol or Ensembl ID of nearest annotated feature |
| `distance` | Distance in kb from peak midpoint to feature start (strand-aware) |
| `within bflag` | `1` if distance ≤ bflag threshold, else `NaN` |
| `epsilon` | Phenotype score from screen data (if `--phenotype` provided) |
| `go=<term>` | `1` if nearest feature is annotated with specified GO term, else `NaN` |

## Notes

- Strand-aware: for `+` strand features, distance is calculated to the feature start; for `−` strand, to the feature end
- GO search requires `--feature_class symbol` (incompatible with `--ensg`)
- For large genomes, runtime scales with the number of peaks × chromosomes; human genome runs in a few minutes on a standard laptop

## Authors

Richardson Lab at the University of California Santa Barbara
https://richardsonlab.mcdb.ucsb.edu/
