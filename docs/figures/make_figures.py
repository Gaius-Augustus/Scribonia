"""Draw the README / docs figures from docs/results/*.tsv (needs matplotlib).

    python3 docs/figures/make_figures.py
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
COLORS = {"TAA,TAG,TGA": "#4c72b0", "TGA": "#dd8452", "TAA,TAG": "#55a868",
          "TAA,TGA": "#c44e52", "TAG": "#8172b3"}
NAMES = {"TAA,TAG,TGA": "standard {TAA,TAG,TGA}", "TGA": "{TGA} (TAA/TAG sense)",
         "TAA,TAG": "{TAA,TAG} (TGA sense)", "TAA,TGA": "{TAA,TGA} (TAG sense)"}
# "Perkinsus": the genome that the results tables list under that name is the
# yeast Candida tropicalis (manifest accession mix-up, corrected 2026-09-27).
NOT_PROTIST = ("Danio", "Drosophila", "Caenorhabditis", "Arabidopsis", "Aspergillus",
               "Saccharomyces", "Candida", "Encephalitozoon", "Perkinsus")


def read(name):
    with open(os.path.join(RES, name)) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def depletion():
    rows = read("depletion_bins.tsv")
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    for key in ("TAA,TAG,TGA", "TGA", "TAA,TAG"):
        r = [x for x in rows if x["stop_set"] == key and x["ambiguous"] == "no"]
        ax.scatter([float(x["D_TAA_TGA"]) for x in r], [float(x["D_TGA_TAATAG"]) for x in r],
                   s=6, alpha=0.5, color=COLORS[key], label=f"{NAMES[key]} (n={len(r)})")
    ax.set_xlabel("TAA depletion D (runs cut at TGA only)")
    ax.set_ylabel("TGA depletion D (runs cut at TAA, TAG only)")
    ax.set_xlim(-0.05, 1.6)
    ax.set_ylim(-0.05, 1.6)
    ax.axhline(0.5, color="grey", lw=0.5, ls=":")
    ax.axvline(0.5, color="grey", lw=0.5, ls=":")
    ax.text(0.02, 0.02, "stop", color="grey", fontsize=8)
    ax.text(1.2, 0.02, "sense", color="grey", fontsize=8)
    ace = [x for x in rows if x["species"].startswith("Acetabularia")]
    if ace:
        ax.annotate("Acetabularia (labelled {TGA},\nlooks standard; not trained on)",
                    xy=(0.22, 0.97), xytext=(0.28, 1.35), fontsize=7,
                    arrowprops=dict(arrowstyle="->", color="grey", lw=0.8))
    ax.set_title("In-frame use inside long open runs, one dot per bin", fontsize=10)
    ax.legend(fontsize=7, loc="upper right", markerscale=2)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "depletion.png"), dpi=300)


def size_curve():
    rows = read("size_curve_v2.tsv")
    mid = [(float(r["size_lo"]) * float(r["size_hi"])) ** 0.5 / 1e6 for r in rows]
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    ax.plot(mid, [float(r["exact"]) for r in rows], "o-", color="#4c72b0", label="exact stop set")
    ax.plot(mid, [float(r["false_nonstandard"]) for r in rows], "s--", color="#c44e52",
            label="standard bins called non-standard")
    ax.axvline(1.0, color="grey", lw=1, ls=":")
    ax.text(1.08, 0.62, "intended use:\nbins >= 1 Mb", fontsize=7, color="grey")
    ax.set_xscale("log")
    ax.set_xlabel("bin size (Mb), real contig lengths")
    ax.set_ylabel("fraction of held-out bins")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Held-out species against bin size (scribonia_v2)", fontsize=10)
    ax.legend(fontsize=7, loc="center right", bbox_to_anchor=(1.0, 0.3))
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "size_curve.png"), dpi=300)


def genus_cv():
    rows = [r for r in read("genus_cv_per_species.tsv")
            if r["ambiguous"] == "no" and not r["species"].startswith(NOT_PROTIST)]
    by = {}
    for r in rows:
        by.setdefault(r["species"], {})[r["test"]] = (float(r["exact"]), r["stop_set"])
    both = {s: v for s, v in by.items() if len(v) == 2}
    shown = sorted((s for s, v in both.items() if min(e for e, _ in v.values()) < 0.995),
                   key=lambda s: both[s]["genus_cv_v3_data"][0])
    perfect = len(both) - len(shown)
    fig, ax = plt.subplots(figsize=(6.2, 0.32 * len(shown) + 1.3))
    for i, s in enumerate(shown):
        e2, key = both[s]["genus_cv_v2_data"]
        e3, _ = both[s]["genus_cv_v3_data"]
        ax.plot([e2, e3], [i, i], color="lightgrey", lw=2, zorder=1)
        ax.scatter([e2], [i], marker="o", facecolor="white", edgecolor=COLORS[key], zorder=2)
        ax.scatter([e3], [i], marker="o", color=COLORS[key], zorder=3)
    ax.set_yticks(range(len(shown)))
    ax.set_yticklabels([f"{s} [{both[s]['genus_cv_v2_data'][1]}]" for s in shown], fontsize=7)
    ax.set_xlabel("exact stop set with the genus left out\n(bins >= 1 Mb, real contig lengths)", fontsize=8)
    ax.set_title(f"Protists below 0.995 (open: v2 genomes, filled: v3 genomes);\n"
                 f"{perfect} further protist species >= 0.995 with both", fontsize=9)
    ax.set_xlim(-0.02, 1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "genus_cv.png"), dpi=300)


if __name__ == "__main__":
    depletion()
    size_curve()
    genus_cv()
