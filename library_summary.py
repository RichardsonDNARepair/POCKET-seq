"""Exact library-scale summary using sufficient statistics for each genomic 12-mer."""
from collections import defaultdict
from itertools import product

import numpy as np

from guide_scan import BASE_CODE, BLOCK, RC, Magnitudes, code, spacer


class WordStats:
    def __init__(self):
        size = 4**12
        self.count = np.zeros(size, dtype=np.int64)
        self.scored = np.zeros(size, dtype=np.int64)
        self.total = np.zeros(size)
        self.best = np.zeros(size)
        self.runner_up = np.zeros(size)

    def add(self, words, context, valid):
        if not words.size:
            return
        # Shift for maxima only. Zero then means no valid observation.
        strength = np.where(valid, 1.0 + context, 0.0)
        order = np.lexsort((strength, words))
        words, strength = words[order], strength[order]
        first = np.r_[0, np.flatnonzero(np.diff(words)) + 1]
        end = np.r_[first[1:], len(words)]
        keys = words[first]
        self.count[keys] += end - first
        self.scored[keys] += np.add.reduceat(valid[order].astype(np.int32), first)
        self.total[keys] += np.add.reduceat(np.where(valid, context, 0.0)[order], first)
        top = strength[end - 1]
        second = np.where(end - first >= 2, strength[np.maximum(end - 2, 0)], 0.0)
        old = self.best[keys]
        self.runner_up[keys] = np.maximum.reduce((self.runner_up[keys], second, np.minimum(old, top)))
        self.best[keys] = np.maximum(old, top)


def pam_blocks(fa, chrom, positions, window, weights):
    length = fa.get_reference_length(chrom)
    for start in range(0, length - 22, BLOCK):
        stop = min(start + BLOCK, length - 22)
        lo, hi = max(0, start - 250), min(length, stop + 272)
        ref = fa.fetch(chrom, lo, hi).upper()
        for plus in (True, False):
            oriented = ref if plus else ref.translate(RC)[::-1]
            dna = BASE_CODE[np.frombuffer(oriented.encode('ascii'), dtype=np.uint8)]
            called = np.r_[np.int32(0), np.cumsum(dna < 4, dtype=np.int32)]
            gc = np.r_[np.int32(0), np.cumsum((dna == 1) | (dna == 2), dtype=np.int32)]
            a, b = (start - lo, stop - lo) if plus else (hi - stop - 22, hi - start - 22)
            i = a + np.flatnonzero((dna[a+21:b+21] == 2) & (dna[a+22:b+22] == 2))
            i = i[called[i+23] - called[i] == 23]
            if not i.size:
                continue
            left, right = np.maximum(0, i-250), np.minimum(len(dna), i+273)
            n = called[i] - called[left] + called[right] - called[i+23]
            g = gc[i] - gc[left] + gc[right] - gc[i+23]
            post_g = (i+23 < len(dna)) & (dna[np.minimum(i+23, len(dna)-1)] == 2)
            context = weights["GC250"] * (g / np.maximum(1, n)) + weights["post1_G"] * post_g
            mid = lo + i + 15.5 if plus else hi - i - 23 + 6.5
            promoter = (np.searchsorted(positions, mid + window, side="right")
                        > np.searchsorted(positions, mid - window))
            whole = np.zeros(len(i), dtype=np.int64)
            for offset in range(20):
                whole = (whole << 2) | dna[i+offset]
            yield whole, context, n > 0, promoter


def count_perfect(query, whole, counts):
    if len(query):
        ix = np.searchsorted(query, whole)
        keep = (ix < len(query)) & (query[np.minimum(ix, len(query)-1)] == whole)
        np.add.at(counts, ix[keep], 1)


def candidate_words(guides):
    suffix = guides & 65535
    exact = (np.arange(256, dtype=np.int32) << 16) | suffix[:, None]
    changed = []
    for p, difference in product(range(8), (1, 2, 3)):
        changed.append(guides ^ (difference << (2*p)))
    return np.column_stack([exact, *changed])


def score_offsets(guides, words, model):
    pairs = [a+b for a in 'ACGT' for b in 'ACGT']
    dna = np.array([model.dna[k] for k in pairs])
    hybrid = np.array([model.hybrid[k] for k in pairs])
    total = -2_000_000 * ((guides[:, None] & 3) != (words & 3)).astype(np.int64)
    for shift in range(0, 22, 2):
        genomic = (words >> shift) & 15
        matched = genomic == ((guides[:, None] >> shift) & 15)
        total += dna[genomic] + 500_000 - matched * hybrid[genomic]
    e = np.sign(total) * ((np.abs(total) + 5_000) // 10_000) / 100
    return model.intercept + model.weights["E12"] * e + model.weights["E12_hinge6"] * np.maximum(0, e-6)


def target_values(fa, target, seq, model, promoters, window):
    chrom, start, strand = target
    lo = max(0, start - 250)
    ref = fa.fetch(chrom, lo, min(fa.get_reference_length(chrom), start + 273)).upper()
    i = start - lo
    site = ref[i:i+23]
    flanks = ref[:i] + ref[i+23:]
    if strand == '+':
        post = ref[i+23:i+24]
    else:
        post = ref[max(0,i-1):i].translate(RC)
        site = site.translate(RC)[::-1]
    n = sum(flanks.count(b) for b in "ACGT")
    context = None
    eta = None
    if n:
        gc = (flanks.count("G") + flanks.count("C")) / n
        context = model.weights["GC250"] * gc + model.weights["post1_G"] * (post == "G")
        eta = model.eta(model.e12(seq, site), gc, post)
    near = promoters.near(chrom, start + (15.5 if strand == '+' else 6.5), window)
    return code(site[8:20]), context, eta, bool(near)


def summarize(fa, contigs, rows, targets, promoters, window, model):
    genome, promoter = WordStats(), WordStats()
    queries = {n: np.array(sorted({code(spacer(r)) for r in rows if len(spacer(r)) == n}), dtype=np.int64)
               for n in (19, 20)}
    matches = {n: np.zeros(len(q), dtype=np.int64) for n, q in queries.items()}
    for chrom in contigs:
        positions = np.array(promoters.positions.get(chrom, ()), dtype=np.int64)
        for whole, context, valid, near in pam_blocks(fa, chrom, positions, window, model.weights):
            words = (whole & 0xffffff).astype(np.int32)
            genome.add(words, context, valid)
            promoter.add(words[near], context[near], valid[near])
            for n, query in queries.items():
                count_perfect(query, whole & ((1 << (2*n))-1), matches[n])
    perfect_by_length = {n: dict(zip(q.tolist(), matches[n].tolist())) for n, q in queries.items()}
    summaries = [{}, {}]
    perfect, excluded = {}, defaultdict(int)
    for a in range(0, len(rows), 256):
        batch = rows[a:a+256]
        guides = np.array([code(spacer(r)[-12:]) for r in batch], dtype=np.int32)
        words = candidate_words(guides)
        offsets = score_offsets(guides, words, model)
        removals = {}
        for j, r in enumerate(batch):
            seq, name = spacer(r), r["guide"]
            perfect[name] = perfect_by_length[len(seq)][code(seq)]
            if name in targets:
                removals[j] = target_values(fa, targets[name], seq, model, promoters, window)
                excluded[name] = 1
        for scope, stats in enumerate((genome, promoter)):
            count = stats.count[words].sum(axis=1)
            n = stats.scored[words]
            scored = n.sum(axis=1)
            total = (n * offsets + stats.total[words]).sum(axis=1)
            best = stats.best[words]
            for j, (word, context, eta, near) in removals.items():
                if scope and not near:
                    continue
                ix, = np.flatnonzero(words[j] == word)
                if not stats.count[word]:
                    raise ValueError("Verified target missing from the genomic 12-mer index")
                count[j] -= 1
                if eta is not None:
                    total[j] -= eta
                    scored[j] -= 1
                    if best[j, ix] == 1 + context:
                        best[j, ix] = stats.runner_up[word]
            maximum = np.where(best > 0, offsets + best - 1, -np.inf).max(axis=1)
            for j, r in enumerate(batch):
                summaries[scope][r["guide"]] = Magnitudes(int(count[j]), int(scored[j]), float(total[j]), float(maximum[j]))
    return summaries[0], summaries[1], perfect, excluded
