"""Stop-set calls of Codetta, PhyloFisher and Scribonia, scored against the manifest.

Expects <results>/<tool>/<sample>/ as written by benchmark/bench.sbatch:
  codetta/     out.genetic_code.out   Codetta inference report
  phylofisher/ counts.tsv             run_phylofisher.py
  scribonia/   out.tsv                scribonia classify
  every tool:  time.tsv               GNU time "%e %U %S %M" (wall, user, sys s, max RSS kB)

Stop-set rules (the tools themselves do not output a stop set):
  * Codetta: a stop codon stays in the stop set unless Codetta infers an
    amino acid for it ('?' = not inferred).  Codetta cannot infer that a
    codon is a stop, so no evidence means standard.
  * PhyloFisher: a stop codon is reassigned when it sits at >= PF_MIN_N
    conserved positions and at >= PF_MIN_FRAC of all conserved positions
    counted in the sample.  These thresholds are ours, not PhyloFisher's
    (it only plots counts); the sweep in pf_sweep.tsv shows how much the
    choice matters, including the best setting found on this very set.
  * Scribonia: the stop_set column.

A run without output is FAIL once it has ended and PENDING before; both
stay out of the rates.  Ended means time.tsv has content (GNU time creates
it empty at the start and writes on exit) and is newer than samples.tsv
(older ones are left over from an earlier submission).  Runs may be scored while
others are still going: summary.tsv gives every tool's rates over all its
finished samples (subset "all") and over the samples finished by every
tool (subset "common"), the fair comparison while one tool lags.

Outputs in <out>: per_sample.tsv, summary.tsv, pf_sweep.tsv.

Usage: python3 -m benchmark.score SAMPLES.tsv RESULTS_DIR OUTDIR
"""

import argparse
import csv
import os
import statistics

CODONS = ("TAA", "TAG", "TGA")
STANDARD = frozenset(CODONS)
TOOLS = ("scribonia", "codetta", "phylofisher")
PF_MIN_N = 5
PF_MIN_FRAC = 1e-3


def fmt(stop_set):
    s = [c for c in CODONS if c in stop_set]
    return ",".join(s) if s else "NA"


def parse_set(text):
    text = (text or "").strip()
    if not text or text.upper() == "NA":
        return frozenset()
    return frozenset(p.strip().upper() for p in text.split(","))


# --------------------------------------------------------------------------- parsers

def read_codetta(path):
    """{codon: (inference, n_aligned, n_used)} from a .genetic_code.out."""
    rows = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) == 6 and len(f[0]) == 3 and set(f[0]) <= set("ACGT"):
                rows[f[0]] = (f[1], int(f[4]), int(f[5]))
    if len(rows) != 64:
        raise ValueError(f"{path}: {len(rows)} codon rows, expected 64")
    return rows


def codetta_call(rows):
    stop = frozenset(c for c in CODONS if not rows[c][0].isalpha())
    evidence = sum(r[2] for r in rows.values())
    detail = ";".join(f"{c}:{rows[c][0]}/{rows[c][2]}" for c in CODONS)
    return stop, evidence, detail


def read_pf_counts(path):
    """{codon: {aa: n}} from run_phylofisher.py's counts.tsv."""
    counts = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            counts.setdefault(r["codon"], {})[r["aa"]] = int(r["n"])
    return counts


def pf_call(counts, min_n=PF_MIN_N, min_frac=PF_MIN_FRAC):
    total = sum(sum(v.values()) for v in counts.values())
    stop = set(CODONS)
    parts = []
    for c in CODONS:
        n = sum(counts.get(c, {}).values())
        if n >= min_n and total and n / total >= min_frac:
            stop.discard(c)
        top = max(counts[c].items(), key=lambda t: t[1])[0] if n else "-"
        parts.append(f"{c}:{top}/{n}")
    return frozenset(stop), total, ";".join(parts)


def read_scribonia(path):
    with open(path) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    return parse_set(row["stop_set"]), row.get("n_orf_eval", ""), row.get("confidence", "")


def read_time(path):
    """(wall s, cpu s, max RSS kB) from GNU time's last line, or Nones."""
    try:
        with open(path) as fh:
            lines = [l.split() for l in fh if l.strip()]
        wall, user, sys_, rss = lines[-1][:4]
        return float(wall), float(user) + float(sys_), int(rss)
    except (OSError, IndexError, ValueError):
        return None, None, None


def _ended(time_tsv, since):
    try:
        return os.path.getsize(time_tsv) > 0 and os.path.getmtime(time_tsv) >= since
    except OSError:
        return False


def call(tool, d):
    """(stop set or None on failure, evidence, detail) for one tool run dir."""
    try:
        if tool == "codetta":
            return codetta_call(read_codetta(os.path.join(d, "out.genetic_code.out")))
        if tool == "phylofisher":
            return pf_call(read_pf_counts(os.path.join(d, "counts.tsv")))
        if tool == "scribonia":
            return read_scribonia(os.path.join(d, "out.tsv"))
    except (OSError, ValueError, StopIteration, KeyError):
        pass
    return None, "", "FAIL"


# --------------------------------------------------------------------------- scoring

def rates(rows):
    """Exact, standard called non-standard, non-standard called standard,
    over the non-ambiguous rows with a call."""
    ok = [r for r in rows if r["ambiguous"] == "no" and r["call"] not in ("FAIL", "PENDING")]
    std = [r for r in ok if r["truth"] == fmt(STANDARD)]
    non = [r for r in ok if r["truth"] != fmt(STANDARD)]

    def frac(sub, pred):
        return round(sum(1 for r in sub if pred(r)) / len(sub), 4) if sub else "NA"
    return {"n": len(ok),
            "exact": frac(ok, lambda r: r["call"] == r["truth"]),
            "std_called_nonstd": frac(std, lambda r: r["call"] != fmt(STANDARD)),
            "nonstd_called_std": frac(non, lambda r: r["call"] == fmt(STANDARD))}


def score(samples_tsv, results, outdir):
    with open(samples_tsv) as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))
    since = os.path.getmtime(samples_tsv)
    per, pf_counts = [], {}
    for s in samples:
        for tool in TOOLS:
            d = os.path.join(results, tool, s["sample"])
            stop, evidence, detail = call(tool, d)
            time_tsv = os.path.join(d, "time.tsv")
            wall, cpu, rss = read_time(time_tsv) if _ended(time_tsv, since) else (None,) * 3
            if stop is not None:
                status = fmt(stop)
            else:
                status = "FAIL" if _ended(os.path.join(d, "time.tsv"), since) else "PENDING"
            per.append({"tool": tool, "sample": s["sample"], "kind": s["kind"],
                        "species": s["species"], "size_bp": s["size_bp"],
                        "contamination": s["contamination"], "ambiguous": s["ambiguous"],
                        "truth": fmt(parse_set(s["stop_set"])),
                        "call": status,
                        "evidence": evidence, "detail": detail,
                        "wall_s": wall, "cpu_s": cpu, "max_rss_kb": rss})
            if tool == "phylofisher" and stop is not None:
                pf_counts[s["sample"]] = read_pf_counts(os.path.join(d, "counts.tsv"))
    os.makedirs(outdir, exist_ok=True)
    _write(os.path.join(outdir, "per_sample.tsv"), per)

    unfinished = {r["sample"] for r in per if r["call"] in ("FAIL", "PENDING")}
    summary = []
    for subset in ("all", "common"):
        for tool in TOOLS:
            for kind in ("genome", "bin"):
                rows = [r for r in per if r["tool"] == tool and r["kind"] == kind
                        and (subset == "all" or r["sample"] not in unfinished)]
                if not rows:
                    continue
                cpu = [r["cpu_s"] for r in rows if r["cpu_s"] is not None]
                wall = [r["wall_s"] for r in rows if r["wall_s"] is not None]
                summary.append(dict(
                    {"subset": subset, "tool": tool, "kind": kind}, **rates(rows),
                    failed=sum(1 for r in rows if r["call"] == "FAIL"),
                    pending=sum(1 for r in rows if r["call"] == "PENDING"),
                    cpu_h_total=round(sum(cpu) / 3600, 2) if cpu else "NA",
                    cpu_s_median=round(statistics.median(cpu), 1) if cpu else "NA",
                    wall_s_median=round(statistics.median(wall), 1) if wall else "NA"))
    _write(os.path.join(outdir, "summary.tsv"), summary)

    # PhyloFisher threshold sweep on the same samples.
    sweep = []
    pf_rows = [r for r in per if r["tool"] == "phylofisher" and r["sample"] in pf_counts]
    for min_n in (1, 2, 5, 10, 20, 50):
        for min_frac in (0, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
            for kind in ("genome", "bin"):
                rows = [dict(r, call=fmt(pf_call(pf_counts[r["sample"]], min_n, min_frac)[0]))
                        for r in pf_rows if r["kind"] == kind]
                if rows:
                    sweep.append(dict({"min_n": min_n, "min_frac": min_frac, "kind": kind},
                                      **rates(rows)))
    _write(os.path.join(outdir, "pf_sweep.tsv"), sweep)
    return per, summary, sweep


def _write(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("samples")
    p.add_argument("results")
    p.add_argument("outdir")
    a = p.parse_args(argv)
    _per, summary, _sweep = score(a.samples, a.results, a.outdir)
    for r in summary:
        print("\t".join(str(v) for v in r.values()))


if __name__ == "__main__":
    main()
