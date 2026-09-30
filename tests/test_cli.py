import csv
import json

from scribonia.cli import OUTPUT_COLUMNS, _load_bins, build_parser, main
from scribonia.fasta import write_fasta
from scribonia.synthetic import make_genome


def test_classify_writes_table(tmp_path):
    fa = tmp_path / "bag_Tetrahymena.fa"
    write_fasta(fa, make_genome({"TGA"}, size_bp=300_000, seed=2, gc=0.30))
    out = tmp_path / "out.tsv"
    js = tmp_path / "out.json"
    main(["classify", str(fa), "--sample", "01", "--bag", "Tetrahymena",
          "--rules", "-o", str(out), "--json", str(js)])
    with open(out) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert list(rows[0]) == OUTPUT_COLUMNS
    assert rows[0]["bag"] == "Tetrahymena"
    assert rows[0]["stop_set"] == "TGA"
    assert rows[0]["table"] == "6"
    assert rows[0]["sample"] == "01"
    d = json.load(open(js))
    assert "features" in d and "summary" in d


def test_per_contig_output(tmp_path):
    fa = tmp_path / "g.fa"
    write_fasta(fa, make_genome({"TAA", "TAG", "TGA"}, size_bp=300_000,
                                n_contigs=6, seed=4))
    out = tmp_path / "o.tsv"
    pc = tmp_path / "pc.tsv"
    main(["classify", str(fa), "--rules", "-o", str(out), "--per-contig", str(pc)])
    lines = open(pc).read().strip().splitlines()
    assert lines[0].startswith("bag\tcontig\tlength")
    assert len(lines) == 7


def test_bag_and_json_need_a_single_input(tmp_path):
    import pytest
    fas = []
    for i in range(2):
        fa = tmp_path / f"g{i}.fa"
        write_fasta(fa, make_genome({"TAA", "TAG", "TGA"}, size_bp=50_000, seed=i))
        fas.append(str(fa))
    for flag in (["--bag", "x"], ["--json", str(tmp_path / "x.json")]):
        with pytest.raises(SystemExit):
            main(["classify", *fas, "--rules", "-o", str(tmp_path / "o.tsv"), *flag])


def test_explicit_missing_model_is_an_error(tmp_path):
    import pytest
    from scribonia.model import load_model
    with pytest.raises(FileNotFoundError):
        load_model(str(tmp_path / "missing.joblib"))


def test_features_npz(tmp_path):
    import numpy as np
    fa = tmp_path / "g.fa"
    write_fasta(fa, make_genome({"TAA", "TAG", "TGA"}, size_bp=200_000, seed=8))
    out = tmp_path / "f.npz"
    main(["features", str(fa), "-o", str(out)])
    d = np.load(out)
    assert d["X"].shape[0] == 1
    assert len(d["feature_names"]) == d["X"].shape[1]


def test_load_bins_concatenates_chunks(tmp_path):
    import numpy as np
    rng = np.random.default_rng(0)
    paths = []
    for i, n in enumerate((4, 6)):
        p = tmp_path / f"train_{i}.npz"
        np.savez(p, X=rng.random((n, 3)), Y=rng.random((n, 3)) > 0.5,
                 groups=np.array([f"g{i}"] * n),
                 meta=np.array([{"size_bp": i}] * n, dtype=object))
        paths.append(str(p))
    d = _load_bins(paths)
    assert d["X"].shape == (10, 3) and d["Y"].shape == (10, 3)
    assert list(d["groups"][:4]) == ["g0"] * 4 and d["meta"][-1]["size_bp"] == 1


def test_train_and_evaluate_take_several_npz():
    p = build_parser()
    assert p.parse_args(["train", "a.npz", "b.npz", "-o", "m"]).features == ["a.npz", "b.npz"]
    assert p.parse_args(["evaluate", "m", "a.npz"]).features == ["a.npz"]
