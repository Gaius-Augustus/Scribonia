import pytest

from scribonia.tables import (CODONS, STOP_SET_TO_TABLE, TABLE_STOP_SETS,
                              canonical_table, format_stop_set, parse_stop_set,
                              stop_set_label, stop_set_of_table)


def test_canonical_tables_round_trip():
    for stop_set, table in STOP_SET_TO_TABLE.items():
        assert stop_set_of_table(table) == stop_set
        assert canonical_table(stop_set) == table


def test_ciliate_tables_share_stop_sets():
    assert stop_set_of_table(6) == stop_set_of_table(27) == stop_set_of_table(29) \
        == stop_set_of_table(30) == frozenset({"TGA"})
    assert stop_set_of_table(10) == stop_set_of_table(4) == frozenset({"TAA", "TAG"})
    assert stop_set_of_table(12) == stop_set_of_table(1)


def test_no_table_for_unusual_sets():
    assert canonical_table(frozenset({"TAA"})) is None
    assert canonical_table(frozenset()) is None


def test_format_parse():
    assert format_stop_set(frozenset(CODONS)) == "TAA,TAG,TGA"
    assert format_stop_set(frozenset()) == "NA"
    assert parse_stop_set("TGA") == frozenset({"TGA"})
    assert parse_stop_set("tag,taa") == frozenset({"TAA", "TAG"})
    assert parse_stop_set("NA") == frozenset()
    with pytest.raises(ValueError):
        parse_stop_set("TTA")


def test_labels():
    assert stop_set_label(frozenset({"TAG", "TAA"})) == "TAA-TAG"
    assert stop_set_label(frozenset()) == "NONE"


def test_unknown_table():
    with pytest.raises(ValueError):
        stop_set_of_table(99)
    assert 1 in TABLE_STOP_SETS
