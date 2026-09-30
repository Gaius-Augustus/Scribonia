"""Synthetic-genome tests of the feature extractor and the rule classifier.

No real data: genomes come from scribonia.synthetic.make_genome, which writes
random genes under a chosen stop set.  The assertions pin down the stop-set
call, the ordering of the depletion statistics, runtime, and the behaviour on
tiny and contaminated bins.
"""

import time

import numpy as np
import pytest

from scribonia.fasta import encode, read_fasta_str
from scribonia.features import CODONS, FEATURE_NAMES, extract_features, feature_vector
from scribonia.rules import classify_summary
from scribonia.synthetic import fragment, make_genome

STD = frozenset(CODONS)


def _genome(stop_set, size=300_000, seed=1, **kw):
    contigs = make_genome(stop_set, size_bp=size, seed=seed, **kw)
    return [encode(seq) for _, seq in contigs]


def _call(seqs):
    feats = extract_features(seqs)
    return feats, classify_summary(feats["summary"])


@pytest.mark.parametrize("stop_set,kw", [
    (STD, dict(gc=0.45)),
    (STD, dict(gc=0.22)),                      # AT-rich: high TAA background
    (STD, dict(gc=0.64)),
    (STD, dict(gc=0.40, intron_len=80, introns_per_gene=4)),
    (frozenset({"TGA"}), dict(gc=0.30)),                       # table 6
    (frozenset({"TGA"}), dict(gc=0.45, reassigned_usage=0.2)),  # rare Gln codons
    (frozenset({"TAA", "TAG"}), dict(gc=0.45)),                # table 10
    (frozenset({"TAA", "TAG"}), dict(gc=0.45, reassigned_usage=0.2)),
    (frozenset({"TAA", "TGA"}), dict(gc=0.40)),                # table 15
    (frozenset({"TAG"}), dict(gc=0.40)),
])
def test_stop_set_is_recovered(stop_set, kw):
    feats, res = _call(_genome(stop_set, **kw))
    assert res["stop_set"] == stop_set, res["notes"]
    for c in CODONS:
        if c in stop_set:
            assert res["p_stop"][c] > 0.6, (c, res["notes"][c])
        else:
            assert res["p_stop"][c] < 0.4, (c, res["notes"][c])


def test_depletion_orders_stop_before_sense():
    """Under the two-codon hypothesis without it, a true stop is depleted
    from run interiors and a reassigned codon is not."""
    feats, _ = _call(_genome(frozenset({"TGA"}), gc=0.30))
    hyp = feats["summary"]["hyp"]
    # TAA and TAG are sense: inside long TGA-runs they occur freely.
    assert hyp["TGA"]["found"]
    assert hyp["TGA"]["TAA"]["D"] > 0.4
    assert hyp["TGA"]["TAG"]["D"] > 0.4
    feats, _ = _call(_genome(STD, gc=0.45))
    hyp = feats["summary"]["hyp"]
    for c in CODONS:
        pair = "".join(d for d in CODONS if d != c)
        assert hyp[pair]["found"]
        assert hyp[pair][c]["D"] < 0.3


def test_coverage_ratio():
    feats, _ = _call(_genome(frozenset({"TGA"}), gc=0.30))
    cr = feats["summary"]["covratio"]
    assert cr["TGA"] == pytest.approx(1.0)
    assert cr["TAA"] < 0.3 and cr["TAG"] < 0.3


def test_feature_vector_is_complete_and_finite():
    feats, _ = _call(_genome(STD))
    v = feature_vector(feats)
    assert v.shape == (len(FEATURE_NAMES),)
    assert np.all(np.isfinite(v))


def test_runtime_under_three_seconds_per_300kb():
    seqs = _genome(STD)
    t0 = time.time()
    extract_features(seqs)
    assert time.time() - t0 < 3.0


def test_tiny_bin_has_low_confidence():
    feats, res = _call(_genome(STD, size=100_000))
    assert res["confidence"] == "low"


def test_fragmented_genome_still_called():
    contigs = make_genome(frozenset({"TGA"}), size_bp=600_000, seed=3, gc=0.32)
    pieces = fragment(contigs, np.random.default_rng(0), median_len=15_000)
    seqs = [encode(seq) for _, seq in pieces]
    assert len(seqs) > 10
    feats, res = _call(seqs)
    assert res["stop_set"] == frozenset({"TGA"}), res["notes"]


def test_contamination_flag_on_mixture():
    """70 % table-6 ciliate plus 30 % standard-code contigs: majority call
    with the contamination flag raised by disagreeing contig votes."""
    cil = make_genome(frozenset({"TGA"}), size_bp=700_000, seed=5, gc=0.30)
    std = make_genome(STD, size_bp=300_000, seed=6, gc=0.45)
    rng = np.random.default_rng(1)
    pieces = fragment(cil, rng, median_len=40_000) + fragment(std, rng, median_len=40_000)
    seqs = [encode(seq) for _, seq in pieces]
    feats, res = _call(seqs)
    assert res["stop_set"] == frozenset({"TGA"}), res["notes"]
    assert res["contamination_flag"], feats["summary"]["votes"]


def test_pure_bin_has_no_contamination_flag():
    contigs = make_genome(STD, size_bp=1_000_000, seed=7, gc=0.45)
    pieces = fragment(contigs, np.random.default_rng(2), median_len=40_000)
    feats, res = _call([encode(seq) for _, seq in pieces])
    assert not res["contamination_flag"], feats["summary"]["votes"]


def test_ns_are_ignored():
    seqs = _genome(STD)
    s = seqs[0].copy()
    s[1000:1100] = 4
    feats, res = _call([s])
    assert res["stop_set"] == STD


def test_fasta_reader_roundtrip():
    text = ">a desc\nACGTNacgt\nAC\n>b\nTTTT\n"
    recs = list(read_fasta_str(text))
    assert [n for n, _ in recs] == ["a", "b"]
    assert recs[0][1].tolist() == [0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1]


def test_scores_survive_extreme_depletion_values():
    """Per-contig D of a short contig can be huge; scoring must not overflow."""
    from scribonia.features import codon_scores
    hyp = {}
    for lab in ("TAA", "TAG", "TGA", "TAATAG", "TAATGA", "TAGTGA", "TAATAGTGA"):
        hyp[lab] = {"n_long": 50, "found": True}
        for c in CODONS:
            if c not in lab:
                hyp[lab][c] = {"D": 1e4}
    scores = codon_scores(hyp, {c: 0.5 for c in CODONS})
    assert all(0.0 <= scores[c][0] <= 1.0 for c in CODONS)


def test_batched_runs_match_contig_by_contig_scan():
    """Contigs are scanned in batches (concatenated with N padding); every
    contig must still get exactly the runs and terminators of a scan of that
    contig alone, including tiny contigs and stops at the contig ends."""
    from scribonia.features import HYPOTHESES, _Strand, _runs, _tri_index
    rng = np.random.default_rng(0)
    lengths = [3, 4, 5, 6, 7, 8, 9] + list(rng.integers(3, 60, size=200))
    # A, G, T only: stops everywhere, many adjacent ones and zero-length runs.
    units = [rng.choice(np.array([0, 2, 3], dtype=np.uint8), size=n) for n in lengths]
    units[3][:] = encode("TAATAG")
    st = _Strand(np.concatenate(units), np.array(lengths, dtype=np.int64),
                 np.arange(len(units)))
    for h in HYPOTHESES:
        for f in range(3):
            is_stop = np.zeros(len(st.codons[f]), dtype=bool)
            for c in h:
                is_stop |= st.hits[f][c]
            starts, ends, unit, stops, term_run = _runs(is_stop, st.base, st.n_codons[f])
            term = np.full(len(starts), -1)
            term[term_run] = st.codons[f][stops]
            for u, seq in enumerate(units):
                tri = _tri_index(seq)[f::3]
                own = np.flatnonzero(np.isin(tri, [16 * 3 + {"TAA": 0, "TAG": 2, "TGA": 8}[c]
                                                   for c in h]))
                sel = unit == u
                assert (starts[sel] - st.base[u]).tolist() == [0] + (own + 1).tolist()
                assert (ends[sel] - st.base[u]).tolist() == own.tolist() + [len(tri)]
                assert term[sel].tolist() == tri[own].tolist() + [-1]
