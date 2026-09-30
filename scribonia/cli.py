"""Command line interface: scribonia classify | features | simulate | recode | train | evaluate."""

import argparse
import json
import os
import sys
import time

import numpy as np

from . import __version__
from .features import CODONS, FEATURE_NAMES, extract_features, feature_vector
from .fasta import read_fasta
from .tables import canonical_table, format_stop_set

OUTPUT_COLUMNS = (
    "sample bag size_bp n_contigs n50 gc gc3 n_orf_eval "
    "p_taa_stop p_tag_stop p_tga_stop d_taa d_tag d_tga "
    "stop_set table confidence ambiguous_codons contamination_flag "
    "minority_vote_frac classifier_version"
).split()


def _load_sequences(path, max_bp=None, seed=0):
    seqs = []
    for _name, arr in read_fasta(path):
        seqs.append(arr)
    if max_bp and sum(len(s) for s in seqs) > max_bp:
        rng = np.random.default_rng(seed)
        order = rng.permutation(len(seqs))
        kept, total = [], 0
        for i in order:
            kept.append(seqs[i])
            total += len(seqs[i])
            if total >= max_bp:
                break
        seqs = kept
    return seqs


def _d_under_best(summary, c):
    from .rules import evidence_for
    _lab, hs, _kind = evidence_for(summary, c)
    return hs[c]["D"] if hs is not None else float("nan")


def classify_file(path, sample="", bag="", model=None, use_rules=False,
                  max_bp=None, per_contig=False, seed=0):
    from .rules import classify_summary
    seqs = _load_sequences(path, max_bp=max_bp, seed=seed)
    feats = extract_features(seqs, per_contig=per_contig)
    summary = feats["summary"]
    if model is not None and not use_rules:
        res = model.classify(feats, summary)
    else:
        res = classify_summary(summary)
    row = {
        "sample": sample, "bag": bag or os.path.basename(path),
        "size_bp": summary["size_bp"], "n_contigs": summary["n_contigs"],
        "n50": summary["n50"], "gc": f"{summary['gc']:.4f}",
        "gc3": f"{summary['gc3']:.4f}", "n_orf_eval": summary["n_orf_eval"],
        "p_taa_stop": f"{res['p_stop']['TAA']:.4f}",
        "p_tag_stop": f"{res['p_stop']['TAG']:.4f}",
        "p_tga_stop": f"{res['p_stop']['TGA']:.4f}",
        "d_taa": f"{_d_under_best(summary, 'TAA'):.3f}",
        "d_tag": f"{_d_under_best(summary, 'TAG'):.3f}",
        "d_tga": f"{_d_under_best(summary, 'TGA'):.3f}",
        "stop_set": res["stop_set_str"],
        "table": res["table"] if res["table"] is not None else "NA",
        "confidence": res["confidence"],
        "ambiguous_codons": ",".join(res["ambiguous"]) or "none",
        "contamination_flag": "yes" if res["contamination_flag"] else "no",
        "minority_vote_frac": f"{res['minority_vote_frac']:.3f}",
        "classifier_version": res["classifier_version"],
    }
    return row, feats, res


def _write_rows(rows, out):
    fh = open(out, "w") if out and out != "-" else sys.stdout
    try:
        fh.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in OUTPUT_COLUMNS) + "\n")
    finally:
        if fh is not sys.stdout:
            fh.close()


# ---------------------------------------------------------------------------

def cmd_classify(args):
    from .model import load_model
    if len(args.fasta) > 1:
        for flag in ("bag", "json"):
            if getattr(args, flag):
                sys.exit(f"[scribonia] --{flag} takes a single input file")
    model = None if args.rules else load_model(args.model)
    if model is None and not args.rules:
        print("[scribonia] no usable trained model; using the rule-based "
              "classifier", file=sys.stderr)
    rows = []
    for path in args.fasta:
        t0 = time.time()
        row, feats, res = classify_file(
            path, sample=args.sample, bag=args.bag, model=model,
            use_rules=args.rules, max_bp=args.max_bp,
            per_contig=bool(args.per_contig), seed=args.seed)
        rows.append(row)
        if args.verbose:
            print(f"[scribonia] {path}: {row['stop_set']} (table {row['table']}, "
                  f"{row['confidence']}) in {time.time() - t0:.1f} s; "
                  + "; ".join(f"{c}: {res['notes'][c]}" for c in CODONS),
                  file=sys.stderr)
        if args.per_contig:
            with open(args.per_contig, "a" if len(rows) > 1 else "w") as fh:
                if len(rows) == 1:
                    fh.write("bag\tcontig\tlength\tstop_set\t"
                             + "\t".join(f"p_{c}\tD_{c}\tnlong_{c}" for c in CODONS)
                             + "\n")
                for r in feats["contigs"]:
                    fh.write(f"{row['bag']}\t{r['contig']}\t{r['length']}\t{r['stop_set']}\t"
                             + "\t".join(f"{r[f'p_{c}']:.3f}\t{r[f'D_{c}']:.3f}\t{r[f'nlong_{c}']}"
                                         for c in CODONS) + "\n")
        if args.json:
            with open(args.json, "w") as fh:
                json.dump({"row": row, "summary": feats["summary"],
                           "features": {k: feats[k] for k in FEATURE_NAMES}},
                          fh, indent=1, default=str)
    _write_rows(rows, args.output)


def cmd_features(args):
    """Feature vectors for many FASTA files -> npz (training input)."""
    X, names = [], []
    for path in args.fasta:
        seqs = _load_sequences(path, max_bp=args.max_bp, seed=args.seed)
        feats = extract_features(seqs)
        X.append(feature_vector(feats))
        names.append(os.path.basename(path))
        if args.verbose:
            print(f"[scribonia] {path}: best hypothesis "
                  f"{feats['summary']['best_hypothesis']}", file=sys.stderr)
    np.savez(args.output, X=np.array(X), names=np.array(names),
             feature_names=np.array(FEATURE_NAMES))


def cmd_simulate(args):
    from .simulate import simulate_bins
    simulate_bins(args.manifest, args.genome_dir, args.output,
                  n_bins=args.n_bins, seed=args.seed,
                  contig_dist=args.contig_dist,
                  size_range=(args.min_mb * 1e6, args.max_mb * 1e6),
                  contig_range=(args.contig_min, args.contig_max),
                  class_mix=args.class_mix, synthetic_share=args.synthetic_share,
                  threads=args.threads, verbose=args.verbose)


def cmd_recode(args):
    from .recode import read_cds, read_fasta_bytes, recode_cds, write_fasta_bytes
    from .tables import parse_stop_set
    seqs = read_fasta_bytes(args.fasta)
    stats = recode_cds(seqs, read_cds(args.gff), parse_stop_set(args.stop_set),
                       args.fraction, seed=args.seed)
    write_fasta_bytes(args.output, seqs)
    print(json.dumps(stats), file=sys.stderr)
    if stats["cds_used"] == 0:
        sys.exit("[scribonia] recode: no usable CDS in " + args.gff)


def _load_bins(paths):
    """Concatenate the npz files of `scribonia simulate` (e.g. array chunks)."""
    parts = [np.load(p, allow_pickle=True) for p in paths]
    keys = [k for k in ("X", "Y", "groups", "meta") if all(k in d for d in parts)]
    return {k: np.concatenate([d[k] for d in parts]) for k in keys}


def cmd_train(args):
    from .model import StopSetModel
    d = _load_bins(args.features)
    X, Y = d["X"], d["Y"]
    groups = d.get("groups")
    model = StopSetModel.fit(X, Y, groups=groups, version=args.version,
                             n_bags=args.n_bags, calibrate=args.calibrate)
    model.save(args.output)
    print(f"[scribonia] saved {args.output} ({args.version})")


def cmd_evaluate(args):
    from .model import StopSetModel, evaluation_report
    model = StopSetModel.load(args.model)
    d = _load_bins(args.features)
    meta = list(d["meta"]) if "meta" in d else [
        {"size_bp": 0, "contamination": 0.0} for _ in d["X"]]
    rep = evaluation_report(model, d["X"], d["Y"], meta, out_json=args.output)
    print(json.dumps(rep, indent=1))


def build_parser():
    p = argparse.ArgumentParser(
        prog="scribonia",
        description="Fast estimation of the genetic-code stop set of genomes "
                    "and metagenomic bins.")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("classify", help="classify one or more FASTA files")
    c.add_argument("fasta", nargs="+")
    c.add_argument("-o", "--output", default="-", help="TSV (default stdout)")
    c.add_argument("--sample", default="",
                   help="free-text label for the sample column (default empty)")
    c.add_argument("--bag", default=None,
                   help="free-text label for the bag column, e.g. bin or species "
                        "name (default: file name; single input only)")
    c.add_argument("--model", default=None, help="model file (.joblib)")
    c.add_argument("--rules", action="store_true",
                   help="use the rule-based classifier even if a model exists")
    c.add_argument("--max-bp", type=int, default=None,
                   help="subsample contigs to at most this many bp")
    c.add_argument("--per-contig", default=None, help="write per-contig votes TSV")
    c.add_argument("--json", default=None, help="dump summary + features (single input)")
    c.add_argument("--seed", type=int, default=0)
    c.add_argument("-v", "--verbose", action="store_true")
    c.set_defaults(func=cmd_classify)

    f = sub.add_parser("features", help="feature vectors of FASTA files -> npz")
    f.add_argument("fasta", nargs="+")
    f.add_argument("-o", "--output", required=True)
    f.add_argument("--max-bp", type=int, default=None)
    f.add_argument("--seed", type=int, default=0)
    f.add_argument("-v", "--verbose", action="store_true")
    f.set_defaults(func=cmd_features)

    s = sub.add_parser("simulate", help="bin-like training samples from labelled genomes")
    s.add_argument("manifest", help="data/manifest.tsv (verified rows only)")
    s.add_argument("genome_dir", help="directory with <accession>/*.fna[.gz]")
    s.add_argument("-o", "--output", required=True, help="npz with X, Y, groups, meta")
    s.add_argument("--n-bins", type=int, default=20000)
    s.add_argument("--contig-dist", default=None,
                   help="TSV of contig lengths to draw from (default log-uniform in the contig range)")
    s.add_argument("--contig-min", type=int, default=2_000, help="shortest contig (bp)")
    s.add_argument("--contig-max", type=int, default=20_000, help="longest contig (bp)")
    s.add_argument("--min-mb", type=float, default=0.3)
    s.add_argument("--max-mb", type=float, default=30.0)
    s.add_argument("--class-mix", default="1:0.40,TGA:0.30,TAA-TAG:0.15,TAA-TGA:0.10,ambiguous:0.05")
    s.add_argument("--synthetic-share", type=float, default=0.5,
                   help="share of a class's bins drawn from recoded genomes")
    s.add_argument("--threads", type=int, default=1,
                   help="worker processes for feature extraction")
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("-v", "--verbose", action="store_true")
    s.set_defaults(func=cmd_simulate)

    r = sub.add_parser("recode", help="rewrite a standard-code genome's CDS into another stop set")
    r.add_argument("fasta", help="genome FASTA (.gz ok)")
    r.add_argument("gff", help="GFF3 with CDS features (NCBI)")
    r.add_argument("--stop-set", required=True, help="target stop set, e.g. TGA or TAA,TAG")
    r.add_argument("--fraction", type=float, required=True,
                   help="share of donor codons changed to each reassigned codon")
    r.add_argument("-o", "--output", required=True)
    r.add_argument("--seed", type=int, default=0)
    r.set_defaults(func=cmd_recode)

    t = sub.add_parser("train", help="fit the gradient-boosting model")
    t.add_argument("features", nargs="+",
                   help="npz from `scribonia simulate` (several are concatenated)")
    t.add_argument("-o", "--output", required=True)
    t.add_argument("--version", dest="version", default="scribonia-dev")
    t.add_argument("--n-bags", type=int, default=10,
                   help="models per codon, each on a bootstrap sample of genera")
    t.add_argument("--calibrate", action="store_true",
                   help="isotonic calibration on leave-genus-out predictions")
    t.set_defaults(func=cmd_train)

    e = sub.add_parser("evaluate", help="metrics of a model on a held-out npz")
    e.add_argument("model")
    e.add_argument("features", nargs="+",
                   help="npz from `scribonia simulate` (several are concatenated)")
    e.add_argument("-o", "--output", default=None, help="JSON report")
    e.set_defaults(func=cmd_evaluate)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
