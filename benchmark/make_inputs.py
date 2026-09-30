"""FASTA inputs for the Codetta / PhyloFisher / Scribonia benchmark.

For every verified primary genome whose species matches --species (the
held-out species of scripts/train.sbatch, so Scribonia never saw
them) and not --exclude:

  * kind "genome": the whole assembly as one plain FASTA;
  * kind "bin": --n-bins bins drawn as in scribonia.simulate.simulate_bins
    (real contig lengths of metagenomic bins from --contig-dist, size log-uniform in
    --min-mb..--max-mb, with probability 0.5 a contamination of
    0.3 * Beta(1, 4) of the bp from a real genome of another stop set).

Writes <outdir>/inputs/<sample>.fa and <outdir>/samples.tsv (one row per
sample: sample kind accession species stop_set ambiguous size_bp
contamination contaminant fasta).  Sample names only use [A-Za-z0-9_.-]:
Codetta rejects other characters in paths.

Usage (compute node, repo root):
  python3 -m benchmark.make_inputs data/manifest.tsv GENOME_DIR OUTDIR \
      --species 'Ichthyophthirius|Stentor|...' --exclude Danio \
      --contig-dist contig_length_dist.tsv --n-bins 5
"""

import argparse
import csv
import os
import re

import numpy as np

from scribonia.simulate import _contig_sampler, draw_bin, load_genome, read_manifest

_DECODE = np.frombuffer(b"ACGTN", dtype=np.uint8)

COLUMNS = ("sample kind accession species stop_set ambiguous size_bp "
           "contamination contaminant fasta").split()


def write_pieces(path, pieces, prefix="c", width=80):
    """uint8-encoded pieces -> FASTA with contigs <prefix>1, <prefix>2, ..."""
    with open(path, "w") as fh:
        for i, arr in enumerate(pieces, 1):
            s = _DECODE[arr].tobytes().decode()
            fh.write(f">{prefix}{i}\n")
            fh.write("\n".join(s[j:j + width] for j in range(0, len(s), width)))
            fh.write("\n")


def safe(text):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


def make_inputs(manifest, genome_dir, outdir, species, exclude=None, n_bins=5,
                min_mb=1.0, max_mb=30.0, contig_dist=None, contig_range=(2_000, 1_000_000),
                contamination_prob=0.5, seed=0):
    rng = np.random.default_rng(seed)
    rows = read_manifest(manifest)
    targets = [r for r in rows if r["role"] == "primary"
               and r.get("label_source") != "synthetic_recode"
               and re.search(species, r["species"])
               and not (exclude and re.search(exclude, r["species"]))]
    if not targets:
        raise ValueError(f"no verified primary genome matches {species!r}")
    os.makedirs(os.path.join(outdir, "inputs"), exist_ok=True)
    sample_len = _contig_sampler(contig_dist, rng, contig_range)
    genomes = {}

    def get(acc):
        if acc not in genomes:
            genomes[acc] = load_genome(genome_dir, acc)
        return genomes[acc]

    out = []
    for r in targets:
        acc = r["accession"]
        base = {"accession": acc, "species": r["species"], "stop_set": r["stop_set"],
                "ambiguous": "yes" if r["ambiguous"] else "no"}
        seqs = get(acc)
        name = f"G_{safe(acc)}"
        fa = os.path.join(outdir, "inputs", name + ".fa")
        write_pieces(fa, seqs)
        out.append(dict(base, sample=name, kind="genome", size_bp=sum(len(s) for s in seqs),
                        contamination=0.0, contaminant="", fasta=fa))
        for k in range(n_bins):
            size = float(np.exp(rng.uniform(np.log(min_mb * 1e6), np.log(max_mb * 1e6))))
            pieces = draw_bin(seqs, size, sample_len, rng, contig_range[0])
            contam, other = 0.0, ""
            if rng.random() < contamination_prob:
                others = [o for o in rows if o["stop_set_fs"] != r["stop_set_fs"]
                          and o.get("label_source") != "synthetic_recode"]
                if others:
                    o = others[int(rng.integers(len(others)))]
                    contam = float(0.3 * rng.beta(1, 4))
                    other = o["accession"]
                    pieces = pieces + draw_bin(get(other), size * contam, sample_len,
                                               rng, contig_range[0])
            name = f"B_{safe(acc)}_{k}"
            fa = os.path.join(outdir, "inputs", name + ".fa")
            write_pieces(fa, pieces)
            out.append(dict(base, sample=name, kind="bin", size_bp=sum(len(p) for p in pieces),
                            contamination=round(contam, 4), contaminant=other, fasta=fa))
        genomes.pop(acc, None)   # primaries are not reused; contaminants stay cached
    with open(os.path.join(outdir, "samples.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("manifest")
    p.add_argument("genome_dir")
    p.add_argument("outdir")
    p.add_argument("--species", required=True, help="regex on the species column")
    p.add_argument("--exclude", default=None, help="regex of species to leave out")
    p.add_argument("--n-bins", type=int, default=5, help="bins per genome")
    p.add_argument("--min-mb", type=float, default=1.0)
    p.add_argument("--max-mb", type=float, default=30.0)
    p.add_argument("--contig-dist", default=None, help="contig lengths, one per line")
    p.add_argument("--contig-min", type=int, default=2_000)
    p.add_argument("--contig-max", type=int, default=1_000_000)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args(argv)
    rows = make_inputs(a.manifest, a.genome_dir, a.outdir, a.species, a.exclude, a.n_bins,
                       a.min_mb, a.max_mb, a.contig_dist, (a.contig_min, a.contig_max),
                       seed=a.seed)
    print(f"[benchmark] {len(rows)} samples in {a.outdir}/samples.tsv")


if __name__ == "__main__":
    main()
