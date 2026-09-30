"""Stop-codon statistics of a genome or bin, computed over six frames.

The idea
--------
A codon that terminates translation never appears in frame inside a gene; a
codon that has been reassigned to an amino acid does.  Scribonia therefore
looks at long open runs between candidate stops and asks, for every one of
TAA / TAG / TGA, whether the run interiors are depleted of it.

Every non-empty subset H of {TAA, TAG, TGA} is a *hypothesis*.  Under H, a
run is a maximal stop-to-stop stretch in one of the six frames.  Random
sequence produces runs whose lengths are geometric with rate p_H (the
trinucleotide frequency of the codons in H), so the number of runs of at
least L codons expected by chance is n_runs * (1 - p_H)^L.  The cutoff L(H)
is raised, within [L_MIN, L_MAX], until that expectation is below 5 % of the
runs actually observed: what remains ("long runs") is then dominated by
genes.  If no such L exists, the hypothesis has no excess of long runs, which
by itself says that H treats at least one sense codon as a stop.

For a hypothesis H and a codon c not in H, the interiors of long H-runs
(non-coding flanks trimmed: 2 / p_H codons at both ends, at most a quarter
of the run) give:

  D_c(H)   in-frame occurrences of c divided by the count expected from the
           position-specific nucleotide composition of the same interiors
           (~0 for a stop, ~1 for a freely used sense codon)
  R_c(H)   in-frame occurrences divided by the mean of the two shifted
           readings of the same nucleotides (composition-free variant)
  SR_c(H)  fraction of long H-runs that stay at least L(H) long once c is
           treated as a stop too: a sense codon chops genes, a stop does not
  SRx_c(H) log SR_c(H) minus the value random sequence would give

plus, per hypothesis, E(H) = log excess of long runs over the random
expectation at L(H), which is high for every subset of the true stop set
and low for any hypothesis that treats a sense codon as a stop.

Long runs are not all genes: in AT-rich or intron-rich genomes many are
non-coding stretches that merely lack the codons of H, and their interiors
make a stop look like a sense codon.  A run "looks coding" when the
nucleotide composition of its interior depends on the codon position (G
statistic of position x nucleotide above the chi-square(6) 0.1 % point),
which genes show whatever their stop codons.  Hence

  codfrac(H)  share of long H-runs that look coding
  Dcod_c(H)   D_c(H) over the coding-looking runs only

The per-contig vote block repeats D_c under the two-codon hypotheses for
every contig of >= 5 kb that carries enough long runs, weighted by length;
disagreement between contigs is the contamination signal.

Only numpy is needed.
"""

from itertools import combinations
import math

import numpy as np

from .fasta import reverse_complement
from .tables import CODONS, format_stop_set

__all__ = [
    "HYPOTHESES", "FEATURE_NAMES", "extract_features", "feature_vector",
    "hypothesis_label", "pair_without",
]

# T=3, A=0, G=2  ->  16*T + 4*x + y
CODON_INDEX = {"TAA": 16 * 3 + 0, "TAG": 16 * 3 + 2, "TGA": 16 * 3 + 4 * 2}
NUC_OF = {"A": 0, "C": 1, "G": 2, "T": 3}

# All non-empty subsets of the three codons, singletons first.
HYPOTHESES = tuple(
    frozenset(h)
    for k in (1, 2, 3)
    for h in combinations(CODONS, k)
)

L_MIN, L_MAX = 100, 1000  # bounds of the hypothesis-specific length cutoff
NULL_SHARE = 0.05         # random long runs must be below this share
TRIM_FLANK_FACTOR = 2.0   # trim 2 / p_H codons from each end of a long run
TRIM_MAX_FRACTION = 0.25  # but never more than this fraction of the run
VOTE_MIN_CONTIG = 5_000   # contigs shorter than this do not vote; longer
                          # ones vote when they carry VOTE_MIN_RUNS long runs
                          # (real metagenomic bins: median contig 10.5 kb, only
                          # 28 % of bins have one of >= 20 kb, all >= 5 kb)
VOTE_MIN_RUNS = 3         # a contig with fewer long runs does not vote
# A long run "looks coding" when the nucleotide composition of its interior
# depends on the codon position: G statistic of position x nucleotide above
# the 0.1 % point of chi-square with 6 degrees of freedom.
CODING_G = 22.46
CHUNK_BP = 2_000_000      # chromosomes are cut into pieces of this size,
                          # contigs are scanned in batches of about this size
MIN_RUNS_FOR_EVIDENCE = 20


def hypothesis_label(h):
    return "".join(c for c in CODONS if c in h)


def pair_without(c):
    """The two-codon hypothesis that leaves c out."""
    return frozenset(set(CODONS) - {c})


# 25a + 5b + c of a trinucleotide (N = 4) -> codon index 0..63, 64 with an N
_TRI_LUT = np.full(125, 64, dtype=np.uint8)
for _a in range(4):
    for _b in range(4):
        for _c in range(4):
            _TRI_LUT[25 * _a + 5 * _b + _c] = 16 * _a + 4 * _b + _c


def _tri_index(s):
    """Codon index (0..63) of every trinucleotide; 64 if it contains an N."""
    s = np.asarray(s, dtype=np.uint8)
    return _TRI_LUT[25 * s[:-2] + 5 * s[1:-1] + s[2:]]


class _Strand:
    """One strand of a batch of units (contigs or pieces of chromosomes),
    scanned in one go instead of unit by unit.

    The units are concatenated, each followed by N padding that makes the
    next one start at a multiple of 3 (at least 3 N, so no codon spans two
    units).  Frame f of every unit is then a slice of frame f of the batch:
    unit u occupies codons base[u] .. base[u] + n_codons[f][u] of
    `codons[f]`, and the padding in between holds no stop codon.
    """

    __slots__ = ("strand", "tri", "codons", "hits", "base", "n_codons",
                 "unit_cid")

    def __init__(self, seq, lengths, unit_cid):
        padded = lengths + (-lengths) % 3 + 3
        offs = np.concatenate(([0], np.cumsum(padded)[:-1]))
        starts = np.concatenate(([0], np.cumsum(lengths)[:-1]))
        total = int(padded.sum())
        self.strand = np.full(total + 2, 4, dtype=np.uint8)
        self.strand[np.arange(len(seq)) + np.repeat(offs - starts, lengths)] = seq
        self.tri = _tri_index(self.strand)         # one per position, 64 in padding
        self.codons = self.tri.reshape(-1, 3).T.copy()     # row f: frame f
        self.hits = [{c: cf == CODON_INDEX[c] for c in CODONS} for cf in self.codons]
        self.base = offs // 3
        self.n_codons = [(lengths - f) // 3 for f in range(3)]   # units of >= 3 bp
        self.unit_cid = unit_cid


def _strands(seqs, batch_bp=CHUNK_BP):
    """Yield (_Strand, forward) for both strands of batches of units."""
    cids, units, total = [], [], 0

    def batch():
        lengths = np.array([len(u) for u in units], dtype=np.int64)
        seq = np.concatenate(units)
        cid = np.array(cids, dtype=np.int64)
        yield _Strand(seq, lengths, cid), True
        yield _Strand(reverse_complement(seq), lengths[::-1].copy(), cid[::-1].copy()), False

    for cid, unit in _iter_units(seqs, batch_bp):
        if len(unit) < 3:
            continue
        cids.append(cid)
        units.append(unit)
        total += len(unit)
        if total >= batch_bp:
            yield from batch()
            cids, units, total = [], [], 0
    if units:
        yield from batch()


def _runs(is_stop, base, n_codons):
    """Maximal stop-to-stop runs in one frame of a batch.

    Every unit contributes its own runs, as if it were scanned alone: its
    stops + 1 of them, the first starting at its first codon, the last ending
    at its last codon without a terminator.  Returns starts, ends (codon
    indices, end exclusive), the unit of every run, the stop positions and
    the run each stop terminates."""
    stops = np.flatnonzero(is_stop)
    unit_end = base + n_codons
    unit_of_stop = np.searchsorted(unit_end, stops, side="right")
    n_stop = np.bincount(unit_of_stop, minlength=len(base))
    first = np.concatenate(([0], np.cumsum(n_stop + 1)[:-1]))
    first_stop = np.concatenate(([0], np.cumsum(n_stop)[:-1]))
    term_run = first[unit_of_stop] + np.arange(len(stops)) - first_stop[unit_of_stop]
    n = len(stops) + len(base)
    starts = np.empty(n, dtype=np.int64)
    ends = np.empty(n, dtype=np.int64)
    starts[first] = base
    starts[term_run + 1] = stops + 1
    ends[term_run] = stops
    ends[first + n_stop] = unit_end
    unit = np.repeat(np.arange(len(base)), n_stop + 1)
    return starts, ends, unit, stops, term_run


class _Hist:
    """Run-length histogram of one hypothesis (pass A).  Lengths above L_MAX
    share the last bin; `lensum` keeps their true lengths so that coverage
    (bp inside runs of at least L codons) stays exact."""

    __slots__ = ("counts", "lensum", "n_runs")

    def __init__(self):
        self.counts = np.zeros(L_MAX + 1, dtype=np.int64)
        self.lensum = np.zeros(L_MAX + 1, dtype=np.int64)
        self.n_runs = 0

    def add(self, lengths):
        capped = np.minimum(lengths, L_MAX)
        self.counts += np.bincount(capped, minlength=L_MAX + 1)
        self.lensum += np.bincount(capped, weights=lengths, minlength=L_MAX + 1).astype(np.int64)
        self.n_runs += len(lengths)

    def n_at_least(self, L):
        return int(self.counts[L:].sum())

    def codons_at_least(self, L):
        return int(self.lensum[L:].sum())


def _choose_L(hist, p_h):
    """Cutoff for one hypothesis: (L, found, E).

    L is the smallest length in [L_MIN, L_MAX] at which the random
    expectation is below NULL_SHARE of the observed count (found=True).  If
    no length qualifies, L is the one with the largest log excess E among
    lengths that still leave MIN_RUNS_FOR_EVIDENCE runs (found=False), so
    that E stays a continuous measure of how far the hypothesis is from
    producing gene-like runs.
    """
    if hist.n_runs == 0:
        return L_MAX, False, 0.0
    tail = np.cumsum(hist.counts[::-1])[::-1]     # tail[L] = n_at_least(L)
    best_L, best_E = L_MAX, -1e9
    for L in range(L_MIN, L_MAX + 1):
        n_obs = tail[L]
        if n_obs < MIN_RUNS_FOR_EVIDENCE and L > L_MIN:
            break
        n_null = hist.n_runs * ((1.0 - p_h) ** L)
        E = math.log((n_obs + 1.0) / (n_null + 1.0))
        if n_null <= NULL_SHARE * n_obs and n_obs >= MIN_RUNS_FOR_EVIDENCE:
            return L, True, E
        if E > best_E:
            best_L, best_E = L, E
    return best_L, False, best_E


class _Acc:
    """Interior statistics of one hypothesis (pass B)."""

    __slots__ = ("n_long", "long_len_sum", "interior_codons", "n0", "n_shift",
                 "expected", "terminators", "gc3_num", "gc3_den",
                 "n_interior", "n_coding", "n0_cod", "exp_cod")

    def __init__(self):
        self.n_long = 0
        self.long_len_sum = 0
        self.interior_codons = 0
        self.n0 = {c: 0 for c in CODONS}
        self.n_shift = {c: 0 for c in CODONS}
        self.expected = {c: 0.0 for c in CODONS}
        self.terminators = {c: 0 for c in CODONS}
        self.gc3_num = 0
        self.gc3_den = 0
        self.n_interior = 0
        self.n_coding = 0
        self.n0_cod = {c: 0 for c in CODONS}
        self.exp_cod = {c: 0.0 for c in CODONS}


def _is_stop(hits, h):
    out = None
    for c in h:
        out = hits[c].copy() if out is None else (out | hits[c])
    return out


def _add_to(arr, cid, weights=None):
    """arr[cid] += weights (or 1), repeated contig ids summed."""
    if len(cid):
        arr += np.bincount(cid, weights=weights, minlength=len(arr)).astype(arr.dtype)


def _pass_b(st, L_of, trim_of, acc, cstats, L_common):
    """Interior statistics of one strand of a batch.  cstats: per-contig
    arrays, cstats[h] = {"n_long", "cov", "n0": {c}, "exp": {c}}."""
    # Prefix sums per residue class r of the batch coordinates: codons and
    # nucleotides at positions r, r + 3, ...  Position j of the codons of
    # frame f, at codon index k, is element k + (f + j) // 3 of residue
    # class (f + j) % 3.
    n_el = len(st.codons[0])
    cs_codon = {c: [np.concatenate(([0], np.cumsum(st.hits[r][c], dtype=np.int32)))
                    for r in range(3)]
                for c in CODONS}
    nuc = [st.strand[r::3][:n_el] for r in range(3)]
    cs_nuc = [[np.concatenate(([0], np.cumsum(nuc[r] == x, dtype=np.int32)))
               for r in range(3)]
              for x in range(4)]
    for h in HYPOTHESES:
        L = L_of[h]
        trim = trim_of[h]
        a = acc[h]
        cst = cstats[h]
        for f in range(3):
            starts, ends, unit, stops, term_run = _runs(
                _is_stop(st.hits[f], h), st.base, st.n_codons[f])
            lengths = ends - starts
            cid = st.unit_cid[unit]
            cov_mask = lengths >= L_common
            _add_to(cst["cov"], cid[cov_mask], lengths[cov_mask])
            long_mask = lengths >= L
            if not long_mask.any():
                continue
            ls, le, ll = starts[long_mask], ends[long_mask], lengths[long_mask]
            lcid = cid[long_mask]
            a.n_long += len(ll)
            a.long_len_sum += int(ll.sum())
            _add_to(cst["n_long"], lcid)

            # Terminators of long runs (the last run of a unit has none).
            term = st.codons[f][stops[long_mask[term_run]]]
            for c in h:
                a.terminators[c] += int(np.count_nonzero(term == CODON_INDEX[c]))

            # Interiors: trim the flanks.
            t = np.minimum(trim, (ll * TRIM_MAX_FRACTION).astype(np.int64))
            ia = ls + t
            ib = le - t
            il = ib - ia
            keep = il > 0
            ia, ib, il, lcid = ia[keep], ib[keep], il[keep], lcid[keep]
            if len(il) == 0:
                continue
            n_int = int(il.sum())
            a.interior_codons += n_int

            comp = np.empty((3, 4, len(il)), dtype=np.int64)
            for j in range(3):
                r, d = (f + j) % 3, (f + j) // 3
                for x in range(4):
                    cn = cs_nuc[x][r]
                    comp[j, x] = cn[ib + d] - cn[ia + d]
            il_f = il.astype(np.float64)
            a.gc3_num += int(comp[2, 1].sum() + comp[2, 2].sum())
            a.gc3_den += n_int

            # Codon-position dependence of the composition (3-periodicity):
            # genes have it whatever their stop codons, random or non-coding
            # stretches that happen to be long runs do not.
            obs = comp.astype(np.float64)
            marg = obs.sum(axis=0, keepdims=True) / 3.0
            with np.errstate(divide="ignore", invalid="ignore"):
                terms = np.where(obs > 0, obs * np.log(obs / marg), 0.0)
            G = 2.0 * terms.sum(axis=(0, 1))
            coding = G > CODING_G
            a.n_interior += len(il)
            a.n_coding += int(coding.sum())

            for c in CODONS:
                if c in h:
                    continue
                x1, x2, x3 = (NUC_OF[ch] for ch in c)
                exp_c = comp[0, x1] * comp[1, x2] * comp[2, x3] / (il_f * il_f)
                cs = cs_codon[c]
                n0, n1, n2 = (cs[(f + j) % 3][ib + (f + j) // 3]
                              - cs[(f + j) % 3][ia + (f + j) // 3] for j in range(3))
                a.n0[c] += int(n0.sum())
                a.n_shift[c] += int(n1.sum() + n2.sum())
                a.expected[c] += float(exp_c.sum())
                a.n0_cod[c] += int(n0[coding].sum())
                a.exp_cod[c] += float(exp_c[coding].sum())
                _add_to(cst["n0"][c], lcid, n0)
                _add_to(cst["exp"][c], lcid, exp_c)


def _iter_units(seqs, chunk_bp=CHUNK_BP):
    for i, s in enumerate(seqs):
        if len(s) <= chunk_bp:
            yield i, s
        else:
            for start in range(0, len(s), chunk_bp):
                yield i, s[start:start + chunk_bp]


def _n50(lengths):
    if not lengths:
        return 0
    ls = sorted(lengths, reverse=True)
    half = sum(ls) / 2.0
    acc = 0
    for l in ls:
        acc += l
        if acc >= half:
            return l
    return ls[-1]


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------

def _feature_names():
    names = []
    for h in HYPOTHESES:
        lab = hypothesis_label(h)
        names += [f"E_{lab}", f"found_{lab}", f"L_{lab}", f"lognlong_{lab}",
                  f"meanlen_{lab}", f"frac_{lab}", f"codfrac_{lab}"]
        for c in CODONS:
            if c in h:
                names.append(f"termshare_{c}_{lab}")
            else:
                names += [f"D_{c}_{lab}", f"R_{c}_{lab}",
                          f"SR_{c}_{lab}", f"SRx_{c}_{lab}", f"Dcod_{c}_{lab}"]
    names += ["gc", "gc3", "p_TAA", "p_TAG", "p_TGA",
              "log_size_bp", "log_n_contigs", "log_n50", "log_n_orf_eval"]
    names += [f"covratio_{c}" for c in CODONS]
    for c in CODONS:
        names += [f"vote_stop_{c}", f"vote_sense_{c}", f"vote_n_{c}",
                  f"vote_iqr_{c}"]
    return tuple(names)


FEATURE_NAMES = _feature_names()


def extract_features(seqs, per_contig=False):
    """Compute the feature dictionary of one bin.

    seqs: iterable of uint8 arrays (A,C,G,T = 0..3, N = 4), one per contig.
    Returns dict with FEATURE_NAMES keys plus 'summary' (the statistics the
    rule-based classifier and the output table use) and, with
    per_contig=True, 'contigs' (per-contig D values).
    """
    seqs = [np.asarray(s, dtype=np.uint8) for s in seqs]
    lengths = [len(s) for s in seqs]
    size_bp = int(sum(lengths))

    # Pass A: composition and run-length histograms.  Nothing is kept for
    # pass B: the batches are cheap to rebuild.
    tri_counts = np.zeros(65, dtype=np.int64)
    nuc_counts = np.zeros(4, dtype=np.int64)
    hists = {h: _Hist() for h in HYPOTHESES}
    for st, forward in _strands(seqs):
        if forward:
            nuc_counts += np.bincount(st.strand, minlength=5)[:4]
        tri_counts += np.bincount(st.tri, minlength=65)
        for h in HYPOTHESES:
            for f in range(3):
                starts, ends, _, _, _ = _runs(_is_stop(st.hits[f], h),
                                              st.base, st.n_codons[f])
                hists[h].add(ends - starts)

    total_tri = max(1, int(tri_counts[:64].sum()))
    p = {c: tri_counts[CODON_INDEX[c]] / total_tri for c in CODONS}
    acgt = nuc_counts.sum()
    gc = float((nuc_counts[1] + nuc_counts[2]) / acgt) if acgt else 0.0

    L_of, found_of, E_of, trim_of = {}, {}, {}, {}
    for h in HYPOTHESES:
        p_h = sum(p[c] for c in h)
        L_of[h], found_of[h], E_of[h] = _choose_L(hists[h], p_h)
        trim_of[h] = int(round(TRIM_FLANK_FACTOR / p_h)) if p_h > 0 else L_MAX
    # One common cutoff for coverage comparisons: the largest among the
    # single-codon hypotheses with an excess of long runs.
    found_singles = [frozenset({c}) for c in CODONS if found_of[frozenset({c})]]
    L_common = max(L_of[h] for h in found_singles) if found_singles else L_MAX

    # Pass B: interiors, plus per-contig statistics for the votes.
    acc = {h: _Acc() for h in HYPOTHESES}
    n_seqs = len(seqs)
    cstats = {h: {"n_long": np.zeros(n_seqs, dtype=np.int64),
                  "cov": np.zeros(n_seqs, dtype=np.int64),
                  "n0": {c: np.zeros(n_seqs, dtype=np.int64) for c in CODONS},
                  "exp": {c: np.zeros(n_seqs) for c in CODONS}}
              for h in HYPOTHESES}
    for st, _ in _strands(seqs):
        _pass_b(st, L_of, trim_of, acc, cstats, L_common)
    contig_stats = {
        cid: {h: {"n_long": int(cstats[h]["n_long"][cid]),
                  "cov": int(cstats[h]["cov"][cid]),
                  "n0": {c: int(cstats[h]["n0"][c][cid]) for c in CODONS},
                  "exp": {c: float(cstats[h]["exp"][c][cid]) for c in CODONS}}
              for h in HYPOTHESES}
        for cid in range(n_seqs) if lengths[cid] >= VOTE_MIN_CONTIG
    }

    feats = {}
    summary = {"p": dict(p), "hyp": {}, "gc": gc, "size_bp": size_bp,
               "n_contigs": len(seqs), "n50": _n50(lengths)}
    best_E, best_h = -1e9, None
    for h in HYPOTHESES:
        lab = hypothesis_label(h)
        a, hist = acc[h], hists[h]
        L = L_of[h]
        p_h = sum(p[c] for c in h)
        n_long = hist.n_at_least(L)
        E = E_of[h]
        feats[f"E_{lab}"] = E
        feats[f"found_{lab}"] = 1.0 if found_of[h] else 0.0
        feats[f"L_{lab}"] = float(L)
        feats[f"lognlong_{lab}"] = math.log1p(n_long)
        feats[f"meanlen_{lab}"] = (a.long_len_sum / a.n_long) if a.n_long else 0.0
        feats[f"frac_{lab}"] = (3.0 * a.long_len_sum / (2.0 * size_bp)) if size_bp else 0.0
        feats[f"codfrac_{lab}"] = (a.n_coding / a.n_interior) if a.n_interior else 0.0
        hs = {"L": L, "found": found_of[h], "n_runs": hist.n_runs,
              "n_long": n_long, "E": E}
        if found_of[h] and E > best_E:
            best_E, best_h = E, h
        n_term = sum(a.terminators.values())
        for c in CODONS:
            if c in h:
                feats[f"termshare_{c}_{lab}"] = (a.terminators[c] / n_term) if n_term else 0.0
                continue
            D = (a.n0[c] + 0.5) / (a.expected[c] + 0.5)
            R = (a.n0[c] + 0.5) / (0.5 * a.n_shift[c] + 0.5)
            n_long_plus = hists[frozenset(h | {c})].n_at_least(L)
            SR = (n_long_plus + 0.5) / (n_long + 0.5)
            null_SR = (1.0 - p[c]) ** L
            SRx = math.log(SR) - math.log(null_SR) if null_SR > 0 else 0.0
            feats[f"D_{c}_{lab}"] = D
            feats[f"R_{c}_{lab}"] = R
            feats[f"SR_{c}_{lab}"] = SR
            feats[f"SRx_{c}_{lab}"] = SRx
            feats[f"Dcod_{c}_{lab}"] = (a.n0_cod[c] + 0.5) / (a.exp_cod[c] + 0.5)
            hs[c] = {"D": D, "R": R, "SR": SR, "SRx": SRx,
                     "n0": a.n0[c], "expected": a.expected[c]}
        summary["hyp"][lab] = hs

    if best_h is None:
        best_h = HYPOTHESES[-1]
    ab = acc[best_h]
    gc3 = (ab.gc3_num / ab.gc3_den) if ab.gc3_den else gc
    feats["gc"] = gc
    feats["gc3"] = gc3
    for c in CODONS:
        feats[f"p_{c}"] = p[c]
    feats["log_size_bp"] = math.log1p(size_bp)
    feats["log_n_contigs"] = math.log1p(len(seqs))
    feats["log_n50"] = math.log1p(summary["n50"])

    # Evidence: for every codon, the largest number of long runs among
    # hypotheses with an excess that do not contain it.
    n_orf_eval = min(
        max([hists[h].n_at_least(L_of[h]) for h in HYPOTHESES
             if c not in h and found_of[h]] + [0])
        for c in CODONS
    )
    feats["log_n_orf_eval"] = math.log1p(n_orf_eval)
    summary["n_orf_eval"] = n_orf_eval

    # Coverage ratio: the share of the genome inside c-free runs of at
    # least L_common codons, relative to the best single codon.  Runs
    # delimited by a true stop cover every gene, whatever the codon's usage
    # as a terminator; runs delimited by a sense codon used at rate q per
    # codon cover only the exp(-q * length) genes that happen to lack it.
    # One common cutoff (the largest among the single-codon hypotheses with
    # an excess) keeps the three coverages comparable.
    cov = {c: (3.0 * hists[frozenset({c})].codons_at_least(L_common) / (2.0 * size_bp))
           if size_bp else 0.0 for c in CODONS}
    cov_max = max(cov.values())
    for c in CODONS:
        feats[f"covratio_{c}"] = (cov[c] / cov_max) if cov_max > 0 else 0.0
    summary["covratio"] = {c: feats[f"covratio_{c}"] for c in CODONS}
    summary["cov"] = cov
    summary["best_hypothesis"] = hypothesis_label(best_h)
    summary["gc3"] = gc3

    # Per-contig votes: every contig of >= VOTE_MIN_CONTIG bp gets its own
    # call with the bin-level cutoffs (codon_scores, the same rule the
    # bin-level classifier uses); the bin then knows whether its contigs
    # agree.  Length-weighted.
    per_contig_rows = []
    weight_stop = {c: 0.0 for c in CODONS}
    weight_sense = {c: 0.0 for c in CODONS}
    weight_n = {c: 0 for c in CODONS}
    set_weight = {}
    D_values = {c: [] for c in CODONS}
    for cid, cs in contig_stats.items():
        chyp = {}
        for h in HYPOTHESES:
            lab = hypothesis_label(h)
            hs = {"n_long": cs[h]["n_long"], "found": found_of[h]}
            for c in CODONS:
                if c in h:
                    continue
                hs[c] = {"D": (cs[h]["n0"][c] + 0.5) / (cs[h]["exp"][c] + 0.5)}
            chyp[lab] = hs
        ccov = {c: cs[frozenset({c})]["cov"] for c in CODONS}
        cmax = max(ccov.values())
        ccovratio = {c: (ccov[c] / cmax) if cmax > 0 else 0.0 for c in CODONS}
        scores = codon_scores(chyp, ccovratio, min_runs=VOTE_MIN_RUNS)
        w = float(lengths[cid])
        row = {"contig": cid, "length": lengths[cid]}
        called = []
        voted = False
        for c in CODONS:
            pc, n, _note = scores[c]
            lab = hypothesis_label(pair_without(c))
            row[f"D_{c}"] = chyp[lab][c]["D"] if chyp[lab]["n_long"] >= VOTE_MIN_RUNS else float("nan")
            row[f"nlong_{c}"] = chyp[lab]["n_long"]
            row[f"p_{c}"] = pc
            if n >= VOTE_MIN_RUNS:
                voted = True
                weight_n[c] += 1
                if pc >= 0.5:
                    weight_stop[c] += w
                    called.append(c)
                else:
                    weight_sense[c] += w
                if not math.isnan(row[f"D_{c}"]):
                    D_values[c].append(row[f"D_{c}"])
            elif pc >= 0.5:
                called.append(c)
        if voted:
            key = format_stop_set(called)
            set_weight[key] = set_weight.get(key, 0.0) + w
            row["stop_set"] = key
        else:
            row["stop_set"] = "NA"
        per_contig_rows.append(row)

    summary["votes"] = {}
    for c in CODONS:
        tot = weight_stop[c] + weight_sense[c]
        stop_frac = (weight_stop[c] / tot) if tot else 0.0
        sense_frac = (weight_sense[c] / tot) if tot else 0.0
        if len(D_values[c]) >= 2:
            q1, q3 = np.percentile(D_values[c], [25, 75])
            iqr = float(q3 - q1)
        else:
            iqr = 0.0
        feats[f"vote_stop_{c}"] = stop_frac
        feats[f"vote_sense_{c}"] = sense_frac
        feats[f"vote_n_{c}"] = float(weight_n[c])
        feats[f"vote_iqr_{c}"] = iqr
        summary["votes"][c] = {"n": weight_n[c], "stop_frac": stop_frac,
                               "sense_frac": sense_frac, "iqr": iqr}
    total_w = sum(set_weight.values())
    if total_w > 0:
        maj_key = max(set_weight, key=set_weight.get)
        summary["majority"] = {
            "stop_set": maj_key,
            "weight_frac": set_weight[maj_key] / total_w,
            "voting_bp_frac": total_w / size_bp if size_bp else 0.0,
            "sets": {k: v / total_w for k, v in set_weight.items()},
        }
    else:
        summary["majority"] = None

    out = {name: float(feats[name]) for name in FEATURE_NAMES}
    out["summary"] = summary
    if per_contig:
        out["contigs"] = per_contig_rows
    return out


# Rule-based scoring shared by the bin-level classifier (rules.py) and the
# per-contig votes above.  D: interior depletion (< 0.2 stop-like, > 0.5
# sense-like); covratio: coverage of long c-free runs relative to the best
# codon (~1 for a stop, exp(-q * gene length) for a sense codon used at rate
# q).  Interiors under a two-codon hypothesis are clean; under a single-codon
# hypothesis adjacent genes merge into one run and D is diluted, so the
# coverage ratio carries more weight there.
D_MID, D_SLOPE = 0.35, 14.0
COV_MID, COV_SLOPE = 0.5, 8.0
W_PAIR = (0.6, 0.4)
W_SINGLE = (0.35, 0.65)


def _sigmoid(x):
    # Stable for any x: per-contig D of a short contig can be far above 1.
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    z = math.exp(x)
    return z / (1.0 + z)


def evidence_for(hyp, c, min_runs=MIN_RUNS_FOR_EVIDENCE):
    """(hypothesis label, stats, kind) used to judge codon c; kind in
    {'pair', 'single', None}.  `hyp` is {label: {"n_long", "found", c: {"D"}}}."""
    pair = hypothesis_label(pair_without(c))
    hs = hyp[pair]
    if hs["found"] and hs["n_long"] >= min_runs:
        return pair, hs, "pair"
    singles = [(hyp[d]["n_long"], d, hyp[d]) for d in CODONS if d != c
               and hyp[d]["found"] and hyp[d]["n_long"] >= min_runs]
    if singles:
        _, lab, hs = max(singles, key=lambda t: t[0])
        return lab, hs, "single"
    return None, None, None


def codon_scores(hyp, covratio, min_runs=MIN_RUNS_FOR_EVIDENCE):
    """{codon: (p_stop, evidence_runs, note)} from hypothesis statistics."""
    out = {}
    have_cov = any(hyp[d]["found"] and hyp[d]["n_long"] >= min_runs for d in CODONS)
    for c in CODONS:
        own = hyp[c]
        cr = covratio.get(c, 0.0)
        p_cov = _sigmoid((cr - COV_MID) * COV_SLOPE) if have_cov else 0.5
        lab, hs, kind = evidence_for(hyp, c, min_runs)
        if hs is not None:
            D = hs[c]["D"]
            p_D = _sigmoid(-(D - D_MID) * D_SLOPE)
            w_d, w_c = W_PAIR if kind == "pair" else W_SINGLE
            p = w_d * p_D + w_c * p_cov
            out[c] = (p, hs["n_long"], f"H={lab}:n={hs['n_long']},D={D:.2f},cov={cr:.2f}")
        elif have_cov:
            out[c] = (p_cov, own["n_long"], f"H={c}-only:n={own['n_long']},cov={cr:.2f}")
        else:
            out[c] = (0.5, 0, "no evidence")
    return out


def feature_vector(feats):
    """Ordered numpy vector for the model."""
    return np.array([feats[n] for n in FEATURE_NAMES], dtype=np.float64)
