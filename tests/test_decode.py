import numpy as np

from scribonia.model import CALLABLE_SETS, EXCLUDED_FEATURES, decode_joint
from scribonia.tables import CODONS, canonical_table


def test_joint_decoding_calls_only_sets_with_a_table():
    rng = np.random.default_rng(0)
    for row in decode_joint(rng.uniform(size=(500, 3))):
        s = frozenset(c for c, b in zip(CODONS, row) if b)
        assert canonical_table(s) is not None
    assert len(CALLABLE_SETS) == 5


def test_joint_decoding_agrees_with_thresholds_when_allowed():
    P = np.array([[0.9, 0.9, 0.9], [0.1, 0.1, 0.9], [0.9, 0.8, 0.1],
                  [0.2, 0.6, 0.9]])    # last: {TAG,TGA} has no table
    got = [",".join(c for c, b in zip(CODONS, r) if b) for r in decode_joint(P)]
    assert got[:3] == ["TAA,TAG,TGA", "TGA", "TAA,TAG"]
    assert got[3] == "TGA"     # 0.8*0.4*0.9 beats 0.2*0.6*0.9 (standard)


def test_composition_features_excluded():
    assert set(EXCLUDED_FEATURES) == {"gc", "gc3", "p_TAA", "p_TAG", "p_TGA"}
