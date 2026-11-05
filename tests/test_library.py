import csv
import random
import sys
import unittest
from collections import defaultdict
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import Promoters
from energy import HeightModel
from guide_scan import code, read_guides, spacer, target_map
from library_summary import candidate_words, summarize
import weissman


def source_row(name, gene, sequence):
    return [name, "h3", "Top5", "", gene, "0", weissman.LEFT + sequence + weissman.RIGHT,
            "", sequence, "P1"]


class Library(unittest.TestCase):
    def test_branch_words_all_positions(self):
        s = "GCGGGGCCGAAG"
        words = candidate_words(np.array([code(s)], dtype=np.int32))[0]
        expected = {code("".join(prefix) + s[-8:]) for prefix in product("ACGT", repeat=4)}
        for i in range(4, 12):
            expected.update(code(s[:i] + base + s[i+1:]) for base in "ACGT" if base != s[i])
        self.assertEqual(len(words), 280)
        self.assertEqual(set(words), expected)

    def test_top5_sequence_and_bulk_summary(self):
        with TemporaryDirectory() as work:
            root = Path(work)
            fasta = root / "test.fa"
            tss = root / "tss.tsv"
            tss.write_text("chrom\tTSS0\tstrand\tgene_id\tgene_name\nchrTest\t500\t+\tID\tGene\n")
            promoters = Promoters(tss)
            rng = random.Random(33)
            bases = list("".join(rng.choices("ACGT", k=1600)))
            seq = "GGCTGCGTGCGGGGCCGAAG"
            t = "A" + seq[1:] + "TGG"
            for i, site in [(0,t),(100,t),(500,t.translate(weissman.RC)[::-1]),(900,seq+"CGG"),
                            (1577,seq[:17]+"T"+seq[18:]+"AGG")]:
                bases[i:i+23] = site
            genome = "".join(bases)
            fasta.write_text(">chrTest\n"+genome+"\n")
            pysam.faidx(str(fasta))
            f = root / "v2.txt"
            raw = [source_row("Gene_-_122.23-P1","Gene","G"+seq[1:].lower()),
                   source_row("Gene_+_500.23-P1","Gene",seq),
                   source_row("non-targeting_00000","negative_control",seq)]
            with f.open("w") as h:
                csv.writer(h, delimiter="\t").writerows(raw)
            _, rows, top5 = read_guides(f)
            self.assertTrue(top5)
            self.assertEqual(rows[0]["guide_sequence"],raw[0][8])
            self.assertEqual([spacer(r) for r in rows], [seq]*3)
            by_strand = [("+", genome), ("-", genome.translate(weissman.RC)[::-1])]
            model = HeightModel()
            with pysam.FastaFile(str(fasta)) as fa:
                weissman.resolve_targets(rows, fa, {"chrTest"}, promoters, "hg19")
                targets = target_map(rows, fa, {"chrTest"}, top5)
                self.assertEqual([targets[r["guide"]][1:] for r in rows[:2]], [(100,"+"),(500,"-")])
                self.assertEqual(rows[2]["intended_mapping"],"non_targeting_control")
                genome_summary, promoter_summary, perfect_counts, excluded = summarize(
                    fa, ["chrTest"], rows, targets, promoters, 100, model)
            for row in rows:
                name = row["guide"]
                expected = defaultdict(list)
                perfect = 0
                for strand, ref in by_strand:
                    for i in range(len(ref)-22):
                        site = ref[i:i+23]
                        if not site.endswith("GG"):
                            continue
                        mm = sum(x != y for x,y in zip(site[12:20],seq[-8:]))
                        if not (mm == 0 or (mm == 1 and site[8:12] == seq[-12:-8])):
                            continue
                        perfect += site[:20] == seq
                        start = i if strand == "+" else len(ref)-i-23
                        if targets.get(name) == ("chrTest",start,strand):
                            continue
                        flank = ref[max(0,i-250):i] + ref[i+23:i+273]
                        gc = sum(flank.count(x) for x in "GC")/len(flank)
                        score = model.eta(model.e12(seq,site),gc,ref[i+23:i+24])
                        expected["genome"].append(score)
                        mid = start + (15.5 if strand == "+" else 6.5)
                        if abs(mid-500) <= 100:
                            expected["promoter"].append(score)
                self.assertEqual(perfect_counts[name],perfect)
                self.assertEqual(excluded[name],int(name in targets))
                for scope, values in [("genome",genome_summary),("promoter",promoter_summary)]:
                    scores = expected[scope]
                    result = values[name]
                    self.assertEqual(result.sites,len(scores))
                    self.assertEqual(result.scored,len(scores))
                    self.assertAlmostEqual(result.total,sum(scores))
                    if scores:
                        self.assertAlmostEqual(result.maximum, max(scores))
                    else:
                        self.assertEqual(result.values(scope)[scope+"_max_eta"], "NA")
            raw[0][6] = raw[0][6][1:]
            with f.open("w") as h:
                csv.writer(h, delimiter="\t").writerows(raw)
            with self.assertRaisesRegex(ValueError, "cloning arms"):
                read_guides(f)
