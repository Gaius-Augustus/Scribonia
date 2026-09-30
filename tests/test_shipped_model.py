"""The shipped model (DEFAULT_MODEL in scribonia/model.py): it is in the
package, loads with the current feature set, and `classify` uses it.  Needs
scikit-learn and joblib and is skipped without them; CI runs it inside the
release image, where a skip is an error."""

import csv
import os

import numpy as np
import pytest

from scribonia.model import DEFAULT_MODEL, MODEL_DIR

VERSION = os.path.splitext(DEFAULT_MODEL)[0]


def test_shipped_model_is_in_the_package():
    assert os.path.getsize(os.path.join(MODEL_DIR, DEFAULT_MODEL)) > 0


def test_shipped_model_loads():
    pytest.importorskip("sklearn")
    pytest.importorskip("joblib")
    from scribonia.features import FEATURE_NAMES
    from scribonia.model import load_model
    m = load_model()
    assert m is not None
    assert m.version == VERSION
    assert m.feature_names == tuple(FEATURE_NAMES)
    assert set(m.heads) == {"TAA", "TAG", "TGA"}


@pytest.mark.parametrize("stop_set,gc,expected", [
    ({"TAA", "TAG", "TGA"}, 0.45, "TAA,TAG,TGA"),
    ({"TGA"}, 0.30, "TGA"),
    ({"TAA", "TAG"}, 0.45, "TAA,TAG"),
])
def test_classify_uses_the_shipped_model(tmp_path, stop_set, gc, expected):
    pytest.importorskip("sklearn")
    pytest.importorskip("joblib")
    from scribonia.cli import main
    from scribonia.fasta import write_fasta
    from scribonia.synthetic import fragment, make_genome
    contigs = fragment(make_genome(stop_set, size_bp=1_500_000, seed=5, gc=gc),
                       np.random.default_rng(5))
    fa = tmp_path / "g.fa"
    write_fasta(fa, contigs)
    out = tmp_path / "out.tsv"
    main(["classify", str(fa), "-o", str(out)])
    with open(out) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    assert row["classifier_version"] == VERSION
    assert row["stop_set"] == expected
