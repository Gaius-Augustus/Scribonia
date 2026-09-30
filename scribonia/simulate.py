"""Bin-like training samples from labelled genomes.

Reads data/manifest.tsv (columns: accession species taxid stop_set table
ambiguous label_source citation verified note role), loads every verified genome from
<genome_dir>/[<label>/]<accession>/*.fna[.gz], and draws bins from the
rows with role "primary":

  * contig lengths from an empirical distribution (TSV of lengths, e.g.
    the contigs of real metagenomic bins) or log-uniform, both within
    contig_range (default 2 .. 20 kb; pass a wider range with an empirical
    distribution, real bins reach several 100 kb);
  * bin size log-uniform in [0.3, 30] Mb;
  * with probability 0.5, contamination of 0.3 * Beta(1, 4) of the bp from
    a random genome of another stop set, primary or role "contaminant"
    (label stays the majority);
  * class mix by stop set as given (`--class-mix`); within a class, real
    and recoded genomes are mixed as `synthetic_share` says.

Only features are stored (X), never sequences, together with the labels
(Y, three booleans), the genus (groups) and metadata (size, contamination).
Recoded standard-code genomes for thin classes (recode.recode_cds, `scribonia
recode`) are made by scripts/train.sbatch and enter here as manifest
rows with label_source synthetic_recode.
"""

import csv
import glob
import multiprocessing as mp
import os

import numpy as np

from .fasta import read_fasta
from .features import CODONS, extract_features, feature_vector
from .tables import parse_stop_set

CLASS_KEYS = {"1": frozenset(CODONS), "TGA": frozenset({"TGA"}),
              "TAA-TAG": frozenset({"TAA", "TAG"}),
              "TAA-TGA": frozenset({"TAA", "TGA"}), "TAG": frozenset({"TAG"}),
              "ambiguous": "ambiguous"}


def read_manifest(path, verified_only=True):
    """Rows of the manifest.  The header line may start with '#' (as in
    data/manifest.tsv); later lines starting with '#' are comments."""
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        body = [l for l in fh if l.strip() and not l.startswith("#")]
    rows = []
    for r in csv.DictReader(body, fieldnames=header, delimiter="\t"):
        if not r.get("accession"):
            continue
        if verified_only and (r.get("verified") or "").strip().lower() not in ("yes", "true", "1"):
            continue
        r["stop_set_fs"] = parse_stop_set(r["stop_set"])
        r["ambiguous"] = (r.get("ambiguous") or "").strip().lower() in ("yes", "true", "1")
        r["role"] = (r.get("role") or "primary").strip().lower() or "primary"
        rows.append(r)
    return rows


def _genome_files(genome_dir, accession):
    pats = ("*.fna", "*.fna.gz", "*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz")
    # <genome_dir>/<accession>/ or, as download_genomes.sh writes it,
    # <genome_dir>/<stop set label>/<accession>/
    files = set()
    for base in (os.path.join(genome_dir, accession),
                 os.path.join(genome_dir, "*", accession)):
        for pat in pats:
            files.update(glob.glob(os.path.join(base, "**", pat), recursive=True))
    return sorted(files)


def load_genome(genome_dir, accession):
    seqs = []
    for f in _genome_files(genome_dir, accession):
        for _n, arr in read_fasta(f):
            seqs.append(arr)
    if not seqs:
        raise FileNotFoundError(f"no FASTA under {genome_dir}/{accession}")
    return seqs


def _contig_sampler(contig_dist, rng, contig_range=(2_000, 20_000)):
    """Contig lengths: drawn from an empirical list (restricted to
    contig_range) or log-uniform in contig_range."""
    lo, hi = contig_range
    if contig_dist:
        lens = np.loadtxt(contig_dist, dtype=np.int64, ndmin=1)
        lens = lens[(lens >= lo) & (lens <= hi)]
        if not len(lens):
            raise ValueError(f"no contig length in {contig_dist} within {lo}..{hi}")
        return lambda: int(rng.choice(lens))
    return lambda: int(np.exp(rng.uniform(np.log(lo), np.log(hi))))


def draw_bin(seqs, target_bp, sample_len, rng, min_len=2_000):
    """Random pieces of random contigs until target_bp is reached; source
    contigs shorter than min_len are not used."""
    pieces, total = [], 0
    lens = np.array([len(s) for s in seqs])
    p = lens / lens.sum()
    guard = 0
    while total < target_bp and guard < 100_000:
        guard += 1
        i = int(rng.choice(len(seqs), p=p))
        s = seqs[i]
        L = min(sample_len(), len(s))
        if L < min_len:
            continue
        start = int(rng.integers(0, len(s) - L + 1))
        pieces.append(s[start:start + L])
        total += L
    return pieces


def _bin_features(pieces):
    feats = extract_features(pieces)
    return feature_vector(feats), feats["summary"]["size_bp"]


def simulate_bins(manifest, genome_dir, output, n_bins=20000, seed=0,
                  contig_dist=None, size_range=(3e5, 3e7), contig_range=(2_000, 20_000),
                  class_mix="1:0.40,TGA:0.30,TAA-TAG:0.15,TAA-TGA:0.10,ambiguous:0.05",
                  contamination_prob=0.5, synthetic_share=0.5, threads=1,
                  verbose=False):
    """synthetic_share: probability that a bin of a class with both real and
    recoded genomes (label_source synthetic_recode) is drawn from a recoded
    one; a class with only one kind always uses that kind.  threads: worker
    processes for feature extraction; the output does not depend on it."""
    rng = np.random.default_rng(seed)
    rows = read_manifest(manifest)
    by_class = {}   # class -> (real rows, synthetic_recode rows)
    for r in rows:
        if r["role"] == "contaminant":
            continue
        key = "ambiguous" if r["ambiguous"] else None
        if key is None:
            for k, fs in CLASS_KEYS.items():
                if fs == r["stop_set_fs"]:
                    key = k
        synth = r.get("label_source") == "synthetic_recode"
        by_class.setdefault(key, ([], []))[int(synth)].append(r)
    mix = {}
    for part in class_mix.split(","):
        k, v = part.split(":")
        if k in by_class:
            mix[k] = float(v)
    if not mix:
        raise ValueError("no verified genome of any requested class")
    keys = list(mix)
    probs = np.array([mix[k] for k in keys])
    probs = probs / probs.sum()

    genomes = {}
    def get(acc):
        if acc not in genomes:
            genomes[acc] = load_genome(genome_dir, acc)
        return genomes[acc]

    sample_len = _contig_sampler(contig_dist, rng, contig_range)

    def draw():
        """One bin: its sequence pieces and label record (main process, so
        the draws do not depend on the number of worker processes)."""
        key = keys[int(rng.choice(len(keys), p=probs))]
        real, synth = by_class[key]
        pool = synth if (not real or (synth and rng.random() < synthetic_share)) else real
        r = pool[int(rng.integers(len(pool)))]
        size = float(np.exp(rng.uniform(np.log(size_range[0]), np.log(size_range[1]))))
        contam = 0.0
        pieces = draw_bin(get(r["accession"]), size, sample_len, rng, contig_range[0])
        if rng.random() < contamination_prob:
            # Real genomes only: a recoded copy of a relative is no contaminant.
            others = [o for o in rows if o["stop_set_fs"] != r["stop_set_fs"]
                      and o.get("label_source") != "synthetic_recode"]
            if others:
                o = others[int(rng.integers(len(others)))]
                contam = float(0.3 * rng.beta(1, 4))
                pieces += draw_bin(get(o["accession"]), size * contam, sample_len, rng, contig_range[0])
        return pieces, r, contam

    X, Y, groups, meta = [], [], [], []
    workers = mp.get_context("fork").Pool(threads) if threads > 1 else None
    try:
        while len(X) < n_bins:
            # Bounded batches: bins are drawn only as fast as they are used.
            batch = [draw() for _ in range(min(4 * max(threads, 1), n_bins - len(X)))]
            pieces = [p for p, _, _ in batch]
            res = (workers.map(_bin_features, pieces, chunksize=1) if workers
                   else [_bin_features(p) for p in pieces])
            for (_, r, contam), (vec, size_bp) in zip(batch, res):
                X.append(vec)
                Y.append([c in r["stop_set_fs"] for c in CODONS])
                groups.append(r["species"].split()[0])
                meta.append({"accession": r["accession"], "species": r["species"],
                             "stop_set": r["stop_set"], "size_bp": size_bp,
                             "contamination": contam, "ambiguous": r["ambiguous"]})
            if verbose:
                print(f"[scribonia] {len(X)}/{n_bins} bins", flush=True)
    finally:
        if workers:
            workers.close()
            workers.join()
    np.savez(output, X=np.array(X), Y=np.array(Y, dtype=bool),
             groups=np.array(groups), meta=np.array(meta, dtype=object))
