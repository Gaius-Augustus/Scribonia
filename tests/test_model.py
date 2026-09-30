"""Bagged model: fit, predict, save/load, report.  Needs scikit-learn and is
skipped without it."""

import numpy as np
import pytest

pytest.importorskip("sklearn")

from scribonia.features import FEATURE_NAMES  # noqa: E402
from scribonia.model import StopSetModel, evaluation_report  # noqa: E402


def _toy(n=600, seed=0):
    """Features where D of a codon under the pair hypothesis says stop (low)
    or sense (high); 12 genera, one stop set each."""
    rng = np.random.default_rng(seed)
    names = list(FEATURE_NAMES)
    sets = [(True, True, True), (False, False, True), (True, True, False)]
    X = rng.normal(0.5, 0.1, size=(n, len(names)))
    Y = np.zeros((n, 3), dtype=bool)
    groups, meta = [], []
    for i in range(n):
        g = i % 12
        y = sets[g % 3]
        Y[i] = y
        for k, (c, pair) in enumerate((("TAA", "TAGTGA"), ("TAG", "TAATGA"), ("TGA", "TAATAG"))):
            X[i, names.index(f"D_{c}_{pair}")] = rng.normal(0.1 if y[k] else 0.9, 0.1)
        groups.append(f"G{g}")
        meta.append({"size_bp": 1e6, "contamination": 0.0, "species": f"G{g} sp",
                     "ambiguous": g == 11})
    return X, Y, np.array(groups), meta


def test_bagged_fit_predict_roundtrip(tmp_path):
    X, Y, groups, meta = _toy()
    m = StopSetModel.fit(X, Y, groups=groups, n_bags=3, min_samples_leaf=20)
    assert all(len(m.heads[c]) >= 1 for c in m.heads)
    P = m.predict_proba(X)
    assert ((P >= 0.5) == Y).mean() > 0.95
    path = tmp_path / "m.joblib"
    m.save(path)
    assert np.allclose(StopSetModel.load(path).predict_proba(X), P)


def test_calibration_option_and_report():
    X, Y, groups, meta = _toy(seed=1)
    m = StopSetModel.fit(X, Y, groups=groups, n_bags=2, min_samples_leaf=20, calibrate=True)
    assert all(m.calibrators[c] is not None for c in m.calibrators)
    rep = evaluation_report(m, X, Y, meta)
    assert rep["n_ambiguous"] == sum(g == "G11" for g in groups)
    assert "G11 sp" in rep["ambiguous_species"] and "G11 sp" not in rep["per_species"]
    assert rep["exact_match"] > 0.9
