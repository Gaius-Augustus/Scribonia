"""Trained classifier: a bag of gradient-boosting models per codon.

scikit-learn is imported lazily so that feature extraction and the rule-based
fallback work in an environment without it (a workflow host with numpy only,
the test suite).  The shipped model lives in scribonia/models/ and is loaded by
name; a missing file makes `load_model` return None and the CLI falls back to
the rules.
"""

import json
import os
import sys

import numpy as np

from .features import CODONS, FEATURE_NAMES
from .rules import confidence_of
from .tables import STOP_SET_TO_TABLE, canonical_table, format_stop_set

MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models")
DEFAULT_MODEL = "scribonia_v3.joblib"

# Features the model does not see: genome composition identifies the taxon
# (and so the training genome) rather than the code, and let the model call
# the AT-richest held-out ciliate by the AT-rich recoded genomes it resembled.
EXCLUDED_FEATURES = ("gc", "gc3", "p_TAA", "p_TAG", "p_TGA")

# Stop sets the model can call: those with an NCBI table.  The three heads
# are decoded jointly (most probable of these sets under independent heads),
# so that a call is always a usable code.
CALLABLE_SETS = tuple(sorted(STOP_SET_TO_TABLE, key=lambda s: (-len(s), sorted(s))))


def decode_joint(P):
    """(n, 3) head probabilities -> (n, 3) booleans, the most probable
    stop set among CALLABLE_SETS."""
    P = np.clip(np.asarray(P, dtype=np.float64), 1e-6, 1 - 1e-6)
    S = np.array([[c in s for c in CODONS] for s in CALLABLE_SETS], dtype=float)
    score = np.log(P) @ S.T + np.log(1 - P) @ (1 - S).T
    return S[score.argmax(axis=1)].astype(bool)


# Monotonic constraints: P(stop) must not increase with D or R and must not
# decrease with SRx or the single-codon excess.  Applied per head to the
# features that concern that head's codon.
def _monotonic_constraints(codon, names=FEATURE_NAMES):
    cst = []
    for name in names:
        if name.startswith((f"D_{codon}_", f"R_{codon}_", f"Dcod_{codon}_")):
            cst.append(-1)
        elif name.startswith(f"SRx_{codon}_") or name == f"E_{codon}":
            cst.append(1)
        else:
            cst.append(0)
    return cst


class StopSetModel:
    """Three binary heads (codon is a stop), each the mean of a bag of
    gradient-boosting models trained on bootstrap samples of genera.

    Bagging over genera, not calibration, is what keeps the calls stable:
    the training set has few genera per non-standard class, and a single
    model's decision boundary moves with whichever genomes the simulation
    happened to draw.  Isotonic calibration on leave-genus-out predictions
    (the v1 default) inherits the systematic errors of those predictions
    and shifted the 0.5 boundary on held-out taxa; it is kept as an option.
    """

    def __init__(self, heads=None, calibrators=None, version="unversioned",
                 feature_names=FEATURE_NAMES, used_features=None):
        # heads[c]: float (constant), a fitted classifier, or a list of them
        self.heads = heads or {}
        self.calibrators = calibrators or {}
        self.version = version
        self.feature_names = tuple(feature_names)
        # Columns of the full feature vector the heads were trained on.
        self.used_features = tuple(used_features or feature_names)
        self._cols = [self.feature_names.index(n) for n in self.used_features]

    # -- training ---------------------------------------------------------
    @classmethod
    def fit(cls, X, Y, groups=None, version="scribonia-dev",
            max_iter=300, max_depth=4, learning_rate=0.05,
            min_samples_leaf=20, l2_regularization=0.0, n_bags=10,
            calibrate=False, n_folds=5, random_state=0,
            exclude=EXCLUDED_FEATURES):
        """X: (n, n_features); Y: (n, 3) booleans (TAA, TAG, TGA is stop);
        groups: genus labels (bootstrap units for the bags and folds for the
        optional calibration)."""
        from sklearn.ensemble import HistGradientBoostingClassifier

        used = tuple(n for n in FEATURE_NAMES if n not in set(exclude))
        X = np.asarray(X, dtype=np.float64)[:, [FEATURE_NAMES.index(n) for n in used]]
        Y = np.asarray(Y, dtype=bool)
        groups = np.asarray(groups) if groups is not None else np.arange(len(X)) // 50
        genera = np.unique(groups)
        rows = {g: np.flatnonzero(groups == g) for g in genera}
        rng = np.random.default_rng(random_state)
        bags = []
        for _ in range(max(1, n_bags)):
            pick = rng.choice(genera, size=len(genera), replace=True) if n_bags > 1 else genera
            bags.append(np.concatenate([rows[g] for g in pick]))

        def fit_head(c, idx, y):
            members = []
            for b in bags:
                i = b[np.isin(b, idx)] if idx is not None else b
                if len(set(y[i])) < 2:
                    continue
                members.append(HistGradientBoostingClassifier(
                    max_iter=max_iter, max_depth=max_depth,
                    learning_rate=learning_rate,
                    min_samples_leaf=min_samples_leaf,
                    l2_regularization=l2_regularization,
                    monotonic_cst=_monotonic_constraints(c, used),
                    random_state=random_state).fit(X[i], y[i]))
            return members

        heads, calibrators = {}, {}
        for k, c in enumerate(CODONS):
            y = Y[:, k]
            if len(set(y)) < 2:
                heads[c] = float(y.mean())
                calibrators[c] = None
                continue
            heads[c] = fit_head(c, None, y) or float(y.mean())
            calibrators[c] = None
            if calibrate:
                calibrators[c] = cls._calibrator(X, y, groups, c, n_folds, fit_head, random_state)
        return cls(heads, calibrators, version=version, used_features=used)

    @staticmethod
    def _calibrator(X, y, groups, c, n_folds, fit_head, random_state):
        from sklearn.isotonic import IsotonicRegression
        from sklearn.model_selection import GroupKFold, KFold
        oof = np.zeros(len(y))
        if len(set(groups)) >= n_folds:
            splitter = GroupKFold(n_splits=n_folds).split(X, y, groups)
        else:
            splitter = KFold(n_splits=n_folds, shuffle=True,
                             random_state=random_state).split(X, y)
        for tr, te in splitter:
            members = fit_head(c, tr, y)
            oof[te] = (np.mean([m.predict_proba(X[te])[:, 1] for m in members], axis=0)
                       if members else float(y[tr].mean()))
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0.001, y_max=0.999)
        return iso.fit(oof, y.astype(float))

    # -- inference --------------------------------------------------------
    def predict_proba(self, X):
        """X: full feature vectors (FEATURE_NAMES order)."""
        X = np.asarray(X, dtype=np.float64)
        if X.ndim == 1:
            X = X[None, :]
        X = X[:, self._cols]
        out = np.zeros((len(X), 3))
        for k, c in enumerate(CODONS):
            h = self.heads[c]
            if isinstance(h, float):
                out[:, k] = h
                continue
            members = h if isinstance(h, list) else [h]
            p = np.mean([m.predict_proba(X)[:, 1] for m in members], axis=0)
            cal = self.calibrators.get(c)
            if cal is not None:
                p = cal.predict(p)
            out[:, k] = np.clip(p, 0.001, 0.999)
        return out

    def classify(self, feats, summary):
        """Result dict in the same shape as rules.classify_summary."""
        from .features import feature_vector
        p = self.predict_proba(feature_vector(feats))[0]
        probs = {c: float(p[k]) for k, c in enumerate(CODONS)}
        n_orf_eval = int(summary.get("n_orf_eval", 0))
        call = decode_joint(p[None, :])[0]
        stop_set = frozenset(c for c, b in zip(CODONS, call) if b)
        ambiguous = [c for c in CODONS
                     if 0.25 <= probs[c] <= 0.75 and n_orf_eval >= 200]
        conf = confidence_of({c: (probs[c], n_orf_eval, "") for c in CODONS},
                             n_orf_eval)
        if summary.get("size_bp", 0) < 500_000 or n_orf_eval < 50:
            conf = "low"
        minority, iqr_max = 0.0, 0.0
        for c in CODONS:
            v = summary["votes"][c]
            if v["n"] >= 2:
                minority = max(minority, min(v["stop_frac"], v["sense_frac"]))
                iqr_max = max(iqr_max, v["iqr"])
        return {
            "p_stop": probs,
            "evidence": {c: n_orf_eval for c in CODONS},
            "notes": {c: "" for c in CODONS},
            "stop_set": stop_set,
            "stop_set_str": format_stop_set(stop_set),
            "table": canonical_table(stop_set),
            "confidence": conf,
            "ambiguous": ambiguous,
            "contamination_flag": minority > 0.15 or iqr_max > 0.4,
            "minority_vote_frac": minority,
            "classifier_version": self.version,
        }

    # -- persistence ------------------------------------------------------
    def save(self, path):
        import joblib
        joblib.dump({"heads": self.heads, "calibrators": self.calibrators,
                     "version": self.version,
                     "feature_names": list(self.feature_names),
                     "used_features": list(self.used_features)}, path)

    @classmethod
    def load(cls, path):
        import joblib
        d = joblib.load(path)
        if tuple(d["feature_names"]) != FEATURE_NAMES:
            raise ValueError(
                f"{path} was trained on a different feature set "
                f"({len(d['feature_names'])} vs {len(FEATURE_NAMES)} features)")
        return cls(d["heads"], d["calibrators"], d["version"],
                   d["feature_names"], d.get("used_features"))


def load_model(path=None):
    """Load the named model or the shipped default; None if not present.

    The shipped default needs scikit-learn and joblib; without them the
    caller falls back to the rules with a warning.  A model named explicitly
    must load."""
    default = path is None
    path = path or os.path.join(MODEL_DIR, DEFAULT_MODEL)
    if not os.path.exists(path):
        if not default:
            raise FileNotFoundError(f"model file not found: {path}")
        return None
    try:
        return StopSetModel.load(path)
    except ImportError as e:
        if not default:
            raise
        print(f"[scribonia] WARNING: {DEFAULT_MODEL} needs scikit-learn and joblib "
              f"({e}); using the rule-based classifier", file=sys.stderr)
        return None


def evaluation_report(model, X, Y, meta, out_json=None):
    """Metrics per codon, per size/contamination stratum and per species.
    meta: list of dicts with 'size_bp', 'contamination' and, from
    `scribonia simulate`, 'species' and 'ambiguous'.  Genomes labelled
    ambiguous (stop codons also read through as sense) are left out of the
    main metrics and reported apart: a stop set does not describe them."""
    from sklearn.metrics import roc_auc_score
    P = model.predict_proba(X)
    Y = np.asarray(Y, dtype=bool)
    amb = np.array([bool(m.get("ambiguous", False)) for m in meta])
    species = np.array([m.get("species", "") for m in meta])
    pred = decode_joint(P)
    exact = (pred == Y).all(axis=1)

    def calls_of(sel):
        calls = {}
        for row in pred[sel]:
            key = format_stop_set(frozenset(c for c, b in zip(CODONS, row) if b)) or "NA"
            calls[key] = calls.get(key, 0) + 1
        return dict(sorted(calls.items(), key=lambda t: -t[1]))

    main = ~amb
    rep = {"n": int(main.sum()), "n_ambiguous": int(amb.sum()),
           "exact_match": float(exact[main].mean()) if main.any() else float("nan"),
           "per_codon": {}, "strata": {}, "per_species": {}, "ambiguous_species": {}}
    for k, c in enumerate(CODONS):
        try:
            auc = float(roc_auc_score(Y[main, k], P[main, k]))
        except ValueError:
            auc = float("nan")
        rep["per_codon"][c] = {
            "auroc": auc,
            "accuracy": float((pred[main, k] == Y[main, k]).mean()) if main.any() else float("nan"),
        }
    size = np.array([m["size_bp"] for m in meta], dtype=float)
    contam = np.array([m.get("contamination", 0.0) for m in meta], dtype=float)
    std = Y.all(axis=1)
    for lo, hi in ((0, 5e5), (5e5, 1e6), (1e6, 5e6), (5e6, 1e12)):
        for clo, chi in ((0, 0.0001), (0.0001, 0.1), (0.1, 0.3)):
            sel = main & (size >= lo) & (size < hi) & (contam >= clo) & (contam < chi)
            if sel.sum() == 0:
                continue
            key = f"size[{int(lo)},{int(hi)})_contam[{clo},{chi})"
            false_ns = float((~pred[sel & std]).any(axis=1).mean()) \
                if (sel & std).sum() else float("nan")
            rep["strata"][key] = {
                "n": int(sel.sum()),
                "exact_match": float(exact[sel].mean()),
                "false_nonstandard_on_table1": false_ns,
            }
    for s in sorted(set(species)):
        sel = species == s
        dest = rep["ambiguous_species"] if amb[sel].any() else rep["per_species"]
        dest[s] = {"n": int(sel.sum()), "exact_match": float(exact[sel].mean()),
                   "calls": calls_of(sel),
                   "mean_p_stop": {c: float(P[sel, k].mean()) for k, c in enumerate(CODONS)}}
    if out_json:
        with open(out_json, "w") as fh:
            json.dump(rep, fh, indent=2)
    return rep
