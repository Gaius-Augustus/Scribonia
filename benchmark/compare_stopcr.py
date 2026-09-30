"""Scribonia's per-codon calls next to stopCR's (Chen et al. 2023, doi:10.1093/molbev/msad064).

stopCR counts two things per codon: occurrences inside conserved coding
regions (codon meaning) and occurrences right after a conserved C-terminus
(termination share).  A codon with an amino-acid meaning and a termination
share of at least TERM_MIN percent is bifunctional (e.g. stichotrich TGA =
Arg with ~98 % of terminations).  Scribonia predicts which codons terminate,
so it should call stop for 'stop' and 'bifunctional' codons and sense for
'sense' ones; p_stop should rise with the termination share.

Codon classes (TERM_MIN is ours: the shares are digitised, +-2 pp):
  stop          stopCR meaning '*'
  bifunctional  amino-acid meaning, termination share >= TERM_MIN
  sense         amino-acid meaning, termination share <  TERM_MIN

Samples are joined to benchmark/chen2023_stopcr.tsv by assembly accession
(number only: GCA/GCF prefix and version ignored); stopCR rows on
transcriptomes (SRR, CRA) have no counterpart and are skipped.

Inputs as for benchmark/score.py (samples.tsv and <results>/<tool>/<sample>/).
Outputs in OUTDIR:
  stopcr_per_codon.tsv  one row per accession x kind x codon
  stopcr_summary.tsv    per kind and class: Scribonia stop-call share, median
                        p_stop and d; Spearman rho of p_stop vs termination share

Usage: python3 -m benchmark.compare_stopcr SAMPLES.tsv RESULTS_DIR OUTDIR \
           [--stopcr benchmark/chen2023_stopcr.tsv]
"""

import argparse
import csv
import os
import re
import statistics

from benchmark.score import read_codetta

CODONS = ("TAA", "TAG", "TGA")
TERM_MIN = 5.0
CLASSES = ("stop", "bifunctional", "sense")
DEFAULT_STOPCR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "chen2023_stopcr.tsv")


def acc_key(acc):
    """'GCF_000165425.1' -> '000165425'; None for non-assembly accessions."""
    m = re.match(r"GC[AF]_(\d+)", (acc or "").strip())
    return m.group(1) if m else None


def read_stopcr(path):
    """{accession key: row} for the stopCR rows with an assembly accession."""
    with open(path) as fh:
        rows = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
        return {acc_key(r["accession"]): r for r in rows if acc_key(r["accession"])}


def codon_class(meaning, term_pct):
    if meaning == "*":
        return "stop"
    return "bifunctional" if term_pct >= TERM_MIN else "sense"


def _num(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def read_scribonia_row(path):
    with open(path) as fh:
        return next(csv.DictReader(fh, delimiter="\t"))


def _median(xs):
    xs = [x for x in xs if x == x]
    return round(statistics.median(xs), 4) if xs else "NA"


def _ranks(xs):
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2
        i = j + 1
    return r


def spearman(x, y):
    pairs = [(a, b) for a, b in zip(x, y) if a == a and b == b]
    if len(pairs) < 3:
        return "NA"
    rx, ry = _ranks([p[0] for p in pairs]), _ranks([p[1] for p in pairs])
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return round(num / den, 4) if den else "NA"


def compare(samples_tsv, results, outdir, stopcr_tsv=DEFAULT_STOPCR):
    stopcr = read_stopcr(stopcr_tsv)
    with open(samples_tsv) as fh:
        samples = [s for s in csv.DictReader(fh, delimiter="\t") if acc_key(s["accession"]) in stopcr]

    # (accession, kind) -> per-sample Scribonia rows, Codetta inference of the genome run
    groups, codetta = {}, {}
    for s in samples:
        key = (s["accession"], s["kind"])
        g = groups.setdefault(key, {"species": s["species"], "label": s["stop_set"], "rows": []})
        try:
            g["rows"].append(read_scribonia_row(os.path.join(results, "scribonia", s["sample"], "out.tsv")))
        except (OSError, StopIteration):
            pass
        if s["kind"] == "genome":
            try:
                inf = read_codetta(os.path.join(results, "codetta", s["sample"], "out.genetic_code.out"))
                codetta[key] = {c: inf[c][0] for c in CODONS}
            except (OSError, ValueError):
                pass

    per = []
    for (acc, kind), g in sorted(groups.items()):
        if not g["rows"]:
            continue
        ref = stopcr[acc_key(acc)]
        for c in CODONS:
            lc = c.lower()
            meaning, term = ref[f"fig4c_{c}"], _num(ref[f"term_pct_{c}"])
            p = [_num(r[f"p_{lc}_stop"]) for r in g["rows"]]
            called = [c in r["stop_set"].split(",") for r in g["rows"]]
            per.append({
                "accession": acc, "species": g["species"], "kind": kind, "manifest_stop_set": g["label"],
                "codon": c, "stopcr_meaning": meaning, "term_pct": ref[f"term_pct_{c}"],
                "cds_n": ref[f"cds_n_{c}"] or "NA", "class": codon_class(meaning, term),
                "n_samples": len(g["rows"]),
                "scribonia_stop_frac": round(sum(called) / len(called), 4),
                "p_stop_median": _median(p),
                "d_median": _median([_num(r[f"d_{lc}"]) for r in g["rows"]]),
                "codetta": codetta.get((acc, kind), {}).get(c, "NA"),
            })

    summary = []
    for kind in ("genome", "bin"):
        rows = [r for r in per if r["kind"] == kind]
        if not rows:
            continue
        for cl in CLASSES:
            sub = [r for r in rows if r["class"] == cl]
            summary.append({
                "kind": kind, "class": cl, "n_codons": len(sub),
                "scribonia_stop_frac": round(statistics.mean(r["scribonia_stop_frac"] for r in sub), 4)
                if sub else "NA",
                "p_stop_median": _median([_num(r["p_stop_median"]) for r in sub]),
                "d_median": _median([_num(r["d_median"]) for r in sub]),
                "spearman_p_stop_vs_term": ""})
        summary.append({
            "kind": kind, "class": "all", "n_codons": len(rows), "scribonia_stop_frac": "",
            "p_stop_median": "", "d_median": "",
            "spearman_p_stop_vs_term": spearman([_num(r["term_pct"]) for r in rows],
                                                [_num(r["p_stop_median"]) for r in rows])})

    os.makedirs(outdir, exist_ok=True)
    _write(os.path.join(outdir, "stopcr_per_codon.tsv"), per)
    _write(os.path.join(outdir, "stopcr_summary.tsv"), summary)
    return per, summary


def _write(path, rows):
    with open(path, "w", newline="") as fh:
        if not rows:
            return
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("samples")
    ap.add_argument("results")
    ap.add_argument("outdir")
    ap.add_argument("--stopcr", default=DEFAULT_STOPCR)
    a = ap.parse_args()
    per, summary = compare(a.samples, a.results, a.outdir, a.stopcr)
    print(f"[compare_stopcr] {len(per)} codon rows")
    for r in summary:
        print("\t".join(str(v) for v in r.values()))


if __name__ == "__main__":
    main()
