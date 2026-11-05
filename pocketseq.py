#!/usr/bin/env python3
"""POCKET-seq: rank observed promoter binding or score guide-sequence candidates."""
import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory

from common import DATA, Promoters, read_contigs


def parser():
    base = argparse.ArgumentParser(description=__doc__)
    modes = base.add_subparsers(dest="command", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--assembly", required=True, help="Reference build label; hg19 enables bundled human defaults")
    common.add_argument("--tss", type=Path, help="TSS TSV[.gz]; required except for hg19")
    common.add_argument("--contigs", help="Contig list file, or 'all'; hg19 defaults to chr1–22,X,Y")
    common.add_argument("--window", type=int, default=1000, help="Inclusive distance to TSS (default 1000 bp)")
    common.add_argument("--out", type=Path, required=True, help="New result directory")
    g = modes.add_parser("guides", parents=[common], help="Genome-wide branched seed search and four-factor scoring")
    g.add_argument("--guides", type=Path, required=True, help="Spacer list, guide_sequence TSV, or original Weissman v2 Top5 file")
    g.add_argument("--fasta", type=Path, required=True, help="Indexed genome FASTA")
    g.add_argument("--emit-sites", choices=["all", "promoter", "none"], help="Detail rows to write (default: none for original Top5, promoter otherwise); summaries cover both scopes")
    p = modes.add_parser("peaks", parents=[common], help="Rank replicate-supported MACS3 promoter peaks for one guide")
    p.add_argument("--guide", required=True, help="Guide/sample label")
    p.add_argument("--replicates", type=Path, nargs="+", required=True, help="Separate narrowPeak or MACS3 text peaks.xls files from independent replicates")
    p.add_argument("--pooled", type=Path, help="Optional pooled calls used for boundaries and representative summit")
    p.add_argument("--q", type=float, default=.05, help="Maximum replicate-call q-value (default 0.05)")
    p.add_argument("--fragments", nargs="+", type=Path, help="One QC-filtered fragment BED per replicate for nCPM")
    p.add_argument("--effective-sizes", nargs="+", type=float, help="Effective retained-fragment library sizes in replicate order")
    p.add_argument("--target", nargs=3, metavar=("CHROM", "START0", "END0"), help="Intended interval in BED coordinates")
    p.add_argument("--pathways", type=Path, help="Gene-to-term TSV[.gz]; required except for hg19")
    return base


def main(argv=None):
    cli = parser()
    args = cli.parse_args(argv)
    try:
        if args.window < 0:
            raise ValueError("window must be nonnegative")
        if args.tss is None:
            if args.assembly != "hg19":
                raise ValueError("Supply --tss for this assembly")
            args.tss = DATA/"RefSeq_Select_hg19.tsv.gz"
        if args.command == "peaks" and args.pathways is None:
            if args.assembly != "hg19":
                raise ValueError("Supply --pathways for this organism")
            args.pathways = DATA/"human_GO_BP.tsv.gz"
        out = args.out
        if out.exists():
            raise ValueError(f"Output already exists: {out}")
        promoters = Promoters(args.tss)
        contigs = read_contigs(args.contigs, args.assembly)
        if contigs is not None and not (set(contigs)&promoters.genes.keys()):
            raise ValueError("Selected contigs and TSS annotation have no names in common")
        if args.command == "guides":
            from guide_scan import run
        else:
            from peak_rank import run
        out.parent.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=out.parent, prefix=".pocketseq-") as tmp:
            stage = Path(tmp)/"result"
            stage.mkdir()
            summary = run(args, stage, promoters, contigs)
            info = dict(completed_utc=datetime.now(timezone.utc).isoformat(), working_directory=str(Path.cwd()),
                        arguments=vars(args), **summary)
            (stage/"run.json").write_text(json.dumps(info, default=str, indent=2)+"\n", encoding="utf-8")
            if out.exists():
                raise ValueError(f"Output already exists: {out}")
            os.rename(stage, out)
        print(out.resolve())
    except (OSError, ValueError, ImportError) as exc:
        cli.exit(2, f"pocketseq: {exc}\n")


if __name__ == "__main__":
    main()
