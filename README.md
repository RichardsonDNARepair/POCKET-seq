# POCKET-seq

POCKET-seq scores candidate dCas9-KRAB binding sites from guide sequences and ranks observed promoter binding from replicate MACS3 peak calls. This release uses the manuscript’s 18,789 hg19 RefSeq Select TSSs and frozen four-component η score.

## Getting started

The examples below assume a terminal opened in this directory. The peak example runs with Python 3 alone. To scan guide sequences or run the tests, install the tested dependencies in a Python 3.12 environment.

```sh
python3 -m pip install -r requirements-tested.txt
python3 -m unittest discover -s tests -v
```

**Score guide sequences.** Supply a plain text file with one 19- or 20-nt guide per line (5′ to 3′, without PAM), or a TSV with a `guide_sequence` column. Provide the hg19 FASTA and its adjacent `.fai` index. This example uses the five in-house guides and writes all genomic candidates with their scores and promoter annotations.

```sh
python3 pocketseq.py guides --assembly hg19 --fasta /path/to/hg19.fa \
  --guides examples/guides.tsv --emit-sites all --out guide_results
```

The scanner requires NGG on either strand and either an exact match at p1–p8, or one mismatch at p1–p8 with exact matches at p9–p12. Position p1 is nearest the PAM. The guide score uses sequence alone. Promoter annotation follows scoring. The manuscript’s Figure 3G applies a separate top-5% selection within each guide’s eligible branched candidates.

**Rank replicated ChIP peaks.** This example ranks promoter peaks by replicate fold enrichment and annotates genes with GO Biological Process terms. To reproduce normalized height in nCPM, also supply filtered fragment BEDs and effective library sizes in replicate order.

```sh
python3 pocketseq.py peaks --assembly hg19 --guide SOCS3g1 \
  --replicates examples/SOCS3g1_rep1.narrowPeak examples/SOCS3g1_rep2.narrowPeak \
  --pooled examples/SOCS3g1_pooled.narrowPeak --target chr17 76356103 76356126 \
  --out peak_results
```

RNA evidence was joined separately for the manuscript. The full usage guide covers alternative assemblies, custom TSS and pathway tables, intended-site handling, optional Weissman Top5 input, and all output fields. [Usage guide](USAGE.txt)

## Annotation sources

The bundled TSS table uses the representative RefSeq transcript start for genes that mapped unambiguously to the GENCODE v19 RNA annotation on hg19 primary nuclear chromosomes. [NCBI annotation release 105.20220307](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Homo_sapiens/105.20220307/) · [GENCODE v19](https://www.gencodegenes.org/human/release_19.html)

The bundled GO lookup is derived from the Gene Ontology Consortium’s human annotations (2023-04-02) and GO Basic ontology (2023-04-01). [Gene Ontology Consortium](https://geneontology.org/)
