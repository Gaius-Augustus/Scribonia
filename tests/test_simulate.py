import os

import numpy as np

from scribonia.fasta import write_fasta
from scribonia.simulate import read_manifest, simulate_bins
from scribonia.synthetic import make_genome

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

HEADER = ("#accession\tspecies\ttaxid\tstop_set\ttable\tambiguous\tlabel_source"
          "\tcitation\tverified\tnote\trole\n")


def test_real_manifest_parses():
    rows = read_manifest(os.path.join(REPO, "data", "manifest.tsv"))
    assert len(rows) > 30
    assert all(r["verified"] == "yes" for r in rows)
    assert {r["role"] for r in rows} == {"primary", "contaminant"}
    assert not any(r["accession"].startswith("#") for r in rows)


def test_simulate_download_layout_and_contaminants(tmp_path):
    # Layout written by scripts/download_genomes.sh: <dir>/<label>/<accession>/
    genomes = {
        "GCA_1": ("TGA", {"TGA"}, "primary", "Tetrahymena thermophila"),
        "GCA_2": ("TAA,TAG,TGA", {"TAA", "TAG", "TGA"}, "primary", "Stentor coeruleus"),
        "GCA_3": ("TAA,TAG", {"TAA", "TAG"}, "contaminant", "Mycoplasma genitalium"),
    }
    lines = [HEADER, "# a comment line\n"]
    for i, (acc, (label, stops, role, species)) in enumerate(genomes.items()):
        d = tmp_path / "genomes" / label.replace(",", "-") / acc
        d.mkdir(parents=True)
        write_fasta(d / "genomic.fna", make_genome(stops, size_bp=300_000, seed=i))
        lines.append(f"{acc}\t{species}\t1\t{label}\t1\tno\tx\t\tyes\t\t{role}\n")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("".join(lines))

    out = tmp_path / "bins.npz"
    simulate_bins(str(manifest), str(tmp_path / "genomes"), str(out), n_bins=12,
                  seed=3, size_range=(1e5, 2e5), contamination_prob=1.0)
    d = np.load(out, allow_pickle=True)
    assert d["X"].shape[0] == 12
    # The contaminant is never the labelled genome of a bin ...
    assert "Mycoplasma" not in set(d["groups"])
    assert not any(tuple(y) == (True, True, False) for y in d["Y"])
    # ... but is mixed into bins of the other stop sets.
    assert all(m["contamination"] > 0 for m in d["meta"])


def test_synthetic_share(tmp_path):
    lines = [HEADER]
    for i, (acc, source) in enumerate((("GCA_real", "literature"),
                                       ("GCA_real_recode", "synthetic_recode"))):
        d = tmp_path / "genomes" / "TGA" / acc
        d.mkdir(parents=True)
        write_fasta(d / "genomic.fna", make_genome({"TGA"}, size_bp=200_000, seed=i))
        lines.append(f"{acc}\tGenus sp\t1\tTGA\t6\tno\t{source}\t\tyes\t\tprimary\n")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("".join(lines))
    for share, want in ((0.0, "GCA_real"), (1.0, "GCA_real_recode")):
        out = tmp_path / f"bins_{share}.npz"
        simulate_bins(str(manifest), str(tmp_path / "genomes"), str(out), n_bins=4,
                      size_range=(1e5, 1.5e5), contamination_prob=0.0,
                      synthetic_share=share)
        d = np.load(out, allow_pickle=True)
        assert {m["accession"] for m in d["meta"]} == {want}


def test_threads_do_not_change_output(tmp_path):
    d = tmp_path / "genomes" / "TGA" / "GCA_1"
    d.mkdir(parents=True)
    write_fasta(d / "genomic.fna", make_genome({"TGA"}, size_bp=200_000, seed=5))
    d = tmp_path / "genomes" / "TAA-TAG-TGA" / "GCA_2"
    d.mkdir(parents=True)
    write_fasta(d / "genomic.fna", make_genome({"TAA", "TAG", "TGA"}, size_bp=200_000, seed=6))
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(HEADER + "GCA_1\tA a\t1\tTGA\t6\tno\tx\t\tyes\t\tprimary\n"
                        "GCA_2\tB b\t1\tTAA,TAG,TGA\t1\tno\tx\t\tyes\t\tprimary\n")
    outs = []
    for threads in (1, 3):
        out = tmp_path / f"b{threads}.npz"
        simulate_bins(str(manifest), str(tmp_path / "genomes"), str(out), n_bins=7,
                      seed=11, size_range=(1e5, 1.5e5), threads=threads)
        outs.append(np.load(out, allow_pickle=True))
    assert np.array_equal(outs[0]["X"], outs[1]["X"])
    assert np.array_equal(outs[0]["Y"], outs[1]["Y"])
