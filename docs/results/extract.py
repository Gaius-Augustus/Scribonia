"""Summarise the evaluation runs of 2026-09-26/27 into the small TSVs next to
this file (needs numpy and scikit-learn; run where the runs live, with
SCRIBONIA_TRAIN_DIR and SCRIBONIA_DEV_DIR set).  docs/figures/make_figures.py
draws the figures from the TSVs.

Inputs, under SCRIBONIA_TRAIN_DIR:
  models/scribonia_v2.joblib, models/scribonia_v3.joblib    official models
  scribonia_v3/chunks/heldout_*.npz    v3 held-out bins (real contig lengths, >= 1 Mb)
under SCRIBONIA_DEV_DIR:
  v7/heldout.npz           held-out bins, real contig lengths, 0.3-30 Mb
  v7s/heldout.npz          held-out bins, real contig lengths, 10-300 kb
  v10/realtest.npz         bins of the v3 training genomes, real contig lengths, >= 1 Mb
  v10_v2/res/cv_oof_P.npy  leave-genus-out predictions for those bins, models fit
  v10_v3/res/cv_oof_P.npy    on v2 or v3 training bins of the other genera
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from scribonia.features import FEATURE_NAMES           # noqa: E402
from scribonia.model import CODONS, StopSetModel, decode_joint  # noqa: E402

DB = os.environ["SCRIBONIA_TRAIN_DIR"]
DEV = os.environ["SCRIBONIA_DEV_DIR"]
OUT = os.path.dirname(os.path.abspath(__file__))
# "Perkinsus": the genome that the results tables list under that name is the
# yeast Candida tropicalis (manifest accession mix-up, corrected 2026-09-27).
NOT_PROTIST = ("Danio", "Drosophila", "Caenorhabditis", "Arabidopsis", "Aspergillus",
               "Saccharomyces", "Candida", "Encephalitozoon", "Perkinsus")


def load(paths):
    ds = [np.load(p, allow_pickle=True) for p in paths]
    return (np.concatenate([d["X"] for d in ds]), np.concatenate([d["Y"] for d in ds]).astype(bool),
            [m for d in ds for m in d["meta"]])


def label(row):
    return ",".join(c for c, b in zip(CODONS, row) if b)


def write(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    print(name, len(rows))


def per_species(test, meta, Y, call, extra=()):
    sp = np.array([m["species"] for m in meta])
    amb = np.array([bool(m["ambiguous"]) for m in meta])
    rows = []
    for s in sorted(set(sp)):
        q = sp == s
        calls = {}
        for r in call[q]:
            calls[label(r)] = calls.get(label(r), 0) + 1
        top = ";".join(f"{k}:{v}" for k, v in sorted(calls.items(), key=lambda t: -t[1]))
        exact = "NA" if amb[q][0] else round(float((call[q] == Y[q]).all(1).mean()), 4)
        rows.append([test, s, label(Y[q][0]), "yes" if amb[q][0] else "no", int(q.sum()), exact, top])
    return rows


def main():
    v2 = StopSetModel.load(f"{DB}/models/scribonia_v2.joblib")
    v3 = StopSetModel.load(f"{DB}/models/scribonia_v3.joblib")

    # 1. accuracy against bin size (v2 model, real contig lengths)
    X, Y, meta = load([f"{DEV}/v7s/heldout.npz", f"{DEV}/v7/heldout.npz"])
    keep = np.array([not m["ambiguous"] for m in meta])
    X, Y, meta = X[keep], Y[keep], [m for m, k in zip(meta, keep) if k]
    ok = (decode_joint(v2.predict_proba(X)) == Y).all(1)
    size = np.array([m["size_bp"] for m in meta], float)
    std = Y.all(1)
    rows = []
    for lo, hi in ((1e4, 3e4), (3e4, 1e5), (1e5, 3e5), (3e5, 1e6), (1e6, 3e6), (3e6, 1e7), (1e7, 3e7)):
        q = (size >= lo) & (size < hi)
        if q.sum():
            call_std = decode_joint(v2.predict_proba(X[q & std])).all(1) if (q & std).any() else np.array([])
            rows.append([int(lo), int(hi), int(q.sum()), round(float(ok[q].mean()), 4),
                         round(float((~call_std).mean()), 4) if call_std.size else "NA"])
    write("size_curve_v2.tsv", ["size_lo", "size_hi", "n_bins", "exact", "false_nonstandard"], rows)

    # 2. leave-genus-out, v2 vs v3 training genomes, same test bins
    X, Y, meta = load([f"{DEV}/v10/realtest.npz"])
    real = np.array(["recoded" not in m["species"] and not m["ambiguous"] for m in meta])
    rows, summ = [], []
    for tag in ("v2", "v3"):
        P = np.load(f"{DEV}/v10_{tag}/res/cv_oof_P.npy")
        call = decode_joint(P)
        rows += per_species(f"genus_cv_{tag}_data", [m for m, r in zip(meta, real) if r], Y[real], call[real])
        prot = real & np.array([not m["species"].startswith(NOT_PROTIST + ("Acetabularia",)) for m in meta])
        s = Y[prot].all(1)
        summ.append([f"genus_cv_{tag}_data", int(prot.sum()), round(float((call[prot] == Y[prot]).all(1).mean()), 4),
                     round(float((~call[prot][s].all(1)).mean()), 4), round(float(call[prot][~s].all(1).mean()), 4)])
    write("genus_cv_per_species.tsv", ["test", "species", "stop_set", "ambiguous", "n_bins", "exact", "calls"], rows)
    write("genus_cv_summary.tsv", ["test", "n_bins_protist", "exact", "false_nonstandard", "nonstandard_missed"], summ)

    # 3. held-out species, v3 model
    import glob
    X, Y, meta = load(sorted(glob.glob(f"{DB}/scribonia_v3/chunks/heldout_*.npz")))
    write("heldout_v3_per_species.tsv", ["test", "species", "stop_set", "ambiguous", "n_bins", "exact", "calls"],
          per_species("heldout_v3", meta, Y, decode_joint(v3.predict_proba(X))))

    # 4. two depletion features per bin (subsample of the leave-genus-out test bins)
    X, Y, meta = load([f"{DEV}/v10/realtest.npz"])
    ix = [FEATURE_NAMES.index(n) for n in ("D_TAA_TGA", "D_TAG_TGA", "D_TGA_TAATAG")]
    real = np.flatnonzero([("recoded" not in m["species"]) for m in meta])
    pick = np.random.default_rng(0).choice(real, size=min(3000, real.size), replace=False)
    write("depletion_bins.tsv", ["species", "stop_set", "ambiguous", "D_TAA_TGA", "D_TAG_TGA", "D_TGA_TAATAG"],
          [[meta[i]["species"], label(Y[i]), "yes" if meta[i]["ambiguous"] else "no",
            *[round(float(X[i, j]), 4) for j in ix]] for i in sorted(pick)])


if __name__ == "__main__":
    main()
