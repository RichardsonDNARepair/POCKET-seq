import csv
import gzip
import io
import math
import random
import sys
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import guide_scan
from energy import HeightModel
from pocketseq import main

GUIDE = "GGCTGCGTGCGGGGCCGAAG"
RC = str.maketrans("ACGT", "TGCA")


def read_rows(path):
    with (gzip.open(path, "rt") if path.suffix == ".gz" else path.open()) as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


class Workflows(unittest.TestCase):
    def setUp(self):
        temp = TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.tss = self.root / "tss.tsv"
        self.tss.write_text("chrom\tTSS0\tstrand\tgene_id\tgene_name\n"
                            "scaffold\t1100\t+\tid1\tGene1\nscaffold\t1101\t-\tid2\tGene2\n")
        self.terms = self.root / "terms.tsv"
        self.terms.write_text("gene_id\tterm_id\tterm_name\nid1\tPW1\tExample pathway\n")

    def run_cmd(self, mode, out, *args):
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            main([mode, "--assembly", "custom", "--tss", str(self.tss), "--out", str(self.root/out), *map(str, args)])

    def test_both_strands_boundaries_and_counts(self):
        rng = random.Random(28)
        ref = list("".join(rng.choices("ACGT", k=2600)))
        one = list(GUIDE)
        one[-4] = "C"
        self.assertNotEqual(one[-4], GUIDE[-4])
        bad = one.copy()
        bad[-9] = "A" if GUIDE[-9] != "A" else "T"
        inserts = [(0, GUIDE+"CGG"), (85, GUIDE+"AGG"), (190, ("".join(one)+"TGG").translate(RC)[::-1]),
                   (360, "".join(bad)+"TGG"), (1100, GUIDE+"TGG"), (1500, "N"+GUIDE[1:]+"AGG"),
                   (2577, (GUIDE+"AGG").translate(RC)[::-1])]
        for pos, seq in inserts:
            ref[pos:pos+23] = seq
        ref = "".join(ref)
        fa = self.root / "ref.fa"
        fa.write_text(">scaffold\n"+ref+"\n")
        pysam.faidx(str(fa))
        inp = self.root / "guides.tsv"
        inp.write_text("guide\tguide_sequence\tbatch\ttarget_chrom\ttarget_start0\ttarget_end0\n"
                       f"g\t{GUIDE}\tA\tscaffold\t85\t108\n"
                       f"copy\t{GUIDE.lower().replace('t','u')}\tB\t\t\t\nshort\t{GUIDE[1:]}\tB\t\t\t\n"
                       "zero\tAAAAAAAAAAAAAAAAAAAA\tB\t\t\t\n")
        with patch.object(guide_scan, "BLOCK", 43):
            self.run_cmd("guides", "result", "--fasta", fa, "--guides", inp, "--emit-sites", "all")
        summary = {r["guide"]: r for r in read_rows(self.root/"result/guide_summary.tsv")}
        sites = { (int(r["start0"]),r["strand"]): r for r in read_rows(self.root/"result/sites.tsv.gz") if r["guide"]=="g" }
        gold = {}
        for start in range(len(ref)-22):
            forward = ref[start:start+23]
            for st, site in [("+",forward),("-",forward.translate(RC)[::-1])]:
                if set(site)-set("ACGT") or not site.endswith("GG"):
                    continue
                mm = sum(a != b for a,b in zip(site[12:20], GUIDE[-8:]))
                if mm==0 or (mm==1 and site[8:12] == GUIDE[-12:-8]):
                    gold[start,st] = site
        self.assertEqual(set(gold), set(sites))
        self.assertIn((190,"-"), gold)
        self.assertNotIn((360,"+"), gold)
        self.assertEqual(len(sites), len([r for r in read_rows(self.root/"result/sites.tsv.gz") if r["guide"]=="g"]))
        promoter, genome = [], []
        for (start,st), seq in gold.items():
            row = sites[start,st]
            self.assertEqual(row["sequence_5to3"], seq)
            flank = ref[max(0,start-250):start]+ref[start+23:start+273]
            ngc=(flank.count("G")+flank.count("C"))/sum(flank.count(b) for b in "ACGT")
            self.assertAlmostEqual(float(row["GC250"]),ngc)
            after = (ref[start+23:start+24] if st=="+" else ref[max(0,start-1):start].translate(RC)) or "N"
            self.assertEqual(row["post_PAM_base"],after)
            mid = start+(15.5 if st=="+" else 6.5)
            near = [f"id{i+1}" for i,p in enumerate([1100,1101]) if abs(mid-p)<=1000]
            self.assertEqual(row["promoter_gene_ids"],";".join(near))
            if (start,st)!=(85,"+"):
                genome.append(float(row["eta"]))
                if near:promoter.append(float(row["eta"]))
        for label, vals in [("genome",genome),("promoter",promoter)]:
            self.assertEqual(int(summary["g"][label+"_sites"]),len(vals))
            self.assertAlmostEqual(float(summary["g"][label+"_sum_eta"]),sum(vals))
            self.assertAlmostEqual(float(summary["g"][label+"_mean_eta"]),sum(vals)/len(vals))
        self.assertEqual(summary["g"]["batch"],"A")
        self.assertEqual(summary["copy"]["guide_sequence"], GUIDE.lower().replace("t", "u"))
        self.assertEqual(int(summary["copy"]["genome_sites"]), len(genome)+1)
        self.assertEqual(summary["copy"]["genome_sum_eta"],summary["short"]["genome_sum_eta"])
        self.assertEqual(summary["zero"]["genome_sites"],"0")
        self.assertEqual(summary["zero"]["genome_mean_eta"],"NA")
        self.assertEqual(summary["zero"]["genome_sum_eta"],"0.0")
        reference_model = HeightModel()
        self.assertEqual(reference_model.e12(GUIDE,GUIDE+"CGG"),9.51)
        self.assertAlmostEqual(reference_model.eta(9.51,.752,"C"),5.7009811201)

    def test_replication_targets_go_and_units(self):
        def narrow(a,b,name,fe,summit,q=4):
            return f"scaffold\t{a}\t{b}\t{name}\t100\t.\t{fe}\t6\t{q}\t{summit-a}\n"
        p1=self.root/"a.narrowPeak.gz"
        with gzip.open(p1,"wt") as f:
            f.write(narrow(90,120,"a",9,100)+narrow(495,505,"b",20,500)+narrow(800,810,"x",15,805)+narrow(1300,1310,"q",20,1305,.8))
        p2=self.root/"b_peaks.xls"
        p2.write_text("# MACS3\nchr\tstart\tend\tabs_summit\tfold_enrichment\t-log10(qvalue)\tname\n"
                      "scaffold\t101\t125\t111\t4\t4\tc\n"
                      "scaffold\t497\t507\t501\t3\t4\td\n"
                      "scaffold\t811\t820\t816\t15\t4\ty\n"
                      "scaffold\t1301\t1310\t1306\t20\t4\tz\n")
        pooled=self.root/"pooled.narrowPeak"
        pooled.write_text(narrow(803,815,"bridge",50,812))
        args=("--guide","g","--replicates",p1,p2,"--pooled",pooled,"--pathways",self.terms,"--target","scaffold","490","513")
        self.run_cmd("peaks","fe",*args)
        rows=read_rows(self.root/"fe/promoter_peaks.tsv")
        self.assertEqual(len(rows),1)
        self.assertEqual((rows[0]["summit0"],rows[0]["signed_distance_bp"]),("100","-1000"))
        self.assertAlmostEqual(float(rows[0]["signal_geomean"]),6)
        self.assertEqual(rows[0]["term_names"],"Example pathway")
        self.assertEqual(len(read_rows(self.root/"fe/intended_promoters.tsv")),2)
        bed1=self.root/"fr1.bed";bed2=self.root/"fr2.bed"
        bed1.write_text("scaffold\t349\t350\nscaffold\t350\t351\nscaffold\t351\t352\n")
        bed2.write_text("scaffold\t0\t1\n")
        self.run_cmd("peaks","counts",*args,"--fragments",bed1,bed2,"--effective-sizes","1000000","1000000")
        row=read_rows(self.root/"counts/promoter_peaks.tsv")[0]
        self.assertEqual((row["rep1_fragments"],row["rep2_fragments"],row["signal_unit"]),("2","1","nCPM"))
        self.assertAlmostEqual(float(row["signal_geomean"]),math.sqrt(2))
        with self.assertRaises(SystemExit):self.run_cmd("peaks","invalid","--guide","g","--replicates",p1,p1,"--pathways",self.terms)
        self.assertFalse((self.root/"invalid").exists())

    def test_three_replicates_need_one_shared_base(self):
        paths = [self.root / f"r{i}.narrowPeak" for i in range(3)]
        for i, ((a, b), path) in enumerate(zip([(100, 120), (110, 130), (120, 140)], paths)):
            path.write_text(f"scaffold\t{a}\t{b}\tchain{i}\t100\t.\t9\t6\t4\t5\n"
                            f"scaffold\t500\t520\tcore{i}\t100\t.\t8\t6\t4\t5\n")
        pooled = self.root / "pooled.narrowPeak"
        pooled.write_text("scaffold\t519\t560\tbridge\t100\t.\t2\t6\t4\t2\n"
                          "scaffold\t550\t570\tsatellite\t100\t.\t100\t6\t4\t16\n")
        self.run_cmd("peaks", "triple", "--guide", "g", "--replicates", *paths, "--pooled", pooled,
                     "--pathways", self.terms)
        rows = read_rows(self.root / "triple/promoter_peaks.tsv")
        self.assertEqual({r["summit0"] for r in rows}, {"505"})
        self.assertEqual({r["core_start0"] for r in rows}, {"500"})
        self.assertEqual({r["locus_id"] for r in rows}, {"g|scaffold:500-570"})

    def test_empty_and_invalid_input(self):
        files=[self.root/"one.narrowPeak",self.root/"two.narrowPeak"]
        for f in files:f.write_text("")
        self.run_cmd("peaks","empty","--guide","g","--replicates",*files,"--pathways",self.terms)
        self.assertEqual(read_rows(self.root/"empty/promoter_peaks.tsv"),[])
        with self.assertRaises(SystemExit):
            self.run_cmd("peaks","empty","--guide","g","--replicates",*files,"--pathways",self.terms)
        files[0].write_text("scaffold\t0\t25\tno_summit\t0\t.\t3\t7\t2\t-1\n")
        with self.assertRaises(SystemExit):
            self.run_cmd("peaks","malformed","--guide","g","--replicates",*files,"--pathways",self.terms)
        self.assertFalse((self.root/"malformed").exists())


if __name__=="__main__":
    unittest.main()
