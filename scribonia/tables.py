"""NCBI genetic-code tables reduced to what a gene finder needs: the stop set.

Scribonia predicts which of the three canonical stop codons (TAA, TAG, TGA)
act as stops in a genome.  Tables that share a stop set are interchangeable
for gene *structure* and differ only in what the reassigned codons encode,
which matters for protein translation but not for finding genes.  The
mapping from a stop set back to an NCBI table therefore picks one canonical
representative; callers that need the amino-acid identity of a reassigned
codon must bring that information from taxonomy.
"""

CODONS = ("TAA", "TAG", "TGA")

# Stop codons per NCBI translation table (nuclear and organellar), restricted
# to the three canonical stops.  Tables 2 (AGA, AGG), 22 (TCA) and 23 (TTA)
# also use other codons as stops; those never occur as nuclear codes and are
# outside the scope of the classifier, but their canonical-stop subset is
# listed so that a lookup never fails.
TABLE_STOP_SETS = {
    1:  frozenset({"TAA", "TAG", "TGA"}),
    2:  frozenset({"TAA", "TAG"}),          # + AGA, AGG (vertebrate mito)
    3:  frozenset({"TAA", "TAG"}),
    4:  frozenset({"TAA", "TAG"}),          # TGA = Trp (Mycoplasma, mito)
    5:  frozenset({"TAA", "TAG"}),
    6:  frozenset({"TGA"}),                 # TAA/TAG = Gln (ciliates)
    9:  frozenset({"TAA", "TAG"}),
    10: frozenset({"TAA", "TAG"}),          # TGA = Cys (Euplotes)
    11: frozenset({"TAA", "TAG", "TGA"}),
    12: frozenset({"TAA", "TAG", "TGA"}),   # CTG = Ser (sense change only)
    13: frozenset({"TAA", "TAG"}),
    14: frozenset({"TAG"}),
    15: frozenset({"TAA", "TGA"}),          # TAG = Gln (Blepharisma, NCBI)
    16: frozenset({"TAA", "TGA"}),
    21: frozenset({"TAA", "TAG"}),
    22: frozenset({"TAA", "TGA"}),          # + TCA
    23: frozenset({"TAA", "TAG", "TGA"}),   # + TTA
    24: frozenset({"TAA", "TAG"}),
    25: frozenset({"TAA", "TAG"}),
    26: frozenset({"TAA", "TAG", "TGA"}),   # CTG = Ala (sense change only)
    27: frozenset({"TGA"}),                 # TAA/TAG = Gln, TGA ambiguous
    28: frozenset({"TAA", "TAG", "TGA"}),   # all three context-dependent
    29: frozenset({"TGA"}),                 # TAA/TAG = Tyr (Mesodinium)
    30: frozenset({"TGA"}),                 # TAA/TAG = Glu (peritrichs)
    31: frozenset({"TAA", "TAG"}),          # TAA/TAG = Glu and stop; TGA = Trp
    33: frozenset({"TAG"}),                 # TAA = Tyr, TGA = Trp
}

# Canonical NCBI table for every stop set.  Chosen as a table that downstream
# gene finders accept for that set.
STOP_SET_TO_TABLE = {
    frozenset({"TAA", "TAG", "TGA"}): 1,
    frozenset({"TGA"}): 6,
    frozenset({"TAA", "TAG"}): 10,
    frozenset({"TAA", "TGA"}): 15,
    frozenset({"TAG"}): 33,
}

# Stop sets with no NCBI table: {TAA}, {TGA, TAG}, {} (Condylostoma-like).


def stop_set_of_table(table):
    """Stop set (frozenset of codons) of an NCBI translation table."""
    try:
        return TABLE_STOP_SETS[int(table)]
    except (KeyError, ValueError, TypeError):
        raise ValueError(f"unknown NCBI translation table: {table!r}")


def canonical_table(stop_set):
    """Canonical NCBI table for a stop set, or None if no table has it."""
    return STOP_SET_TO_TABLE.get(frozenset(stop_set))


def format_stop_set(stop_set):
    """'TAA,TAG,TGA' in fixed order, or 'NA' for the empty set."""
    s = [c for c in CODONS if c in stop_set]
    return ",".join(s) if s else "NA"


def parse_stop_set(text):
    """Inverse of format_stop_set; accepts 'NA', '', 'TGA', 'TAA,TAG'."""
    if text is None:
        return frozenset()
    text = text.strip()
    if not text or text.upper() == "NA":
        return frozenset()
    parts = {p.strip().upper() for p in text.split(",") if p.strip()}
    bad = parts - set(CODONS)
    if bad:
        raise ValueError(f"not a stop codon: {sorted(bad)}")
    return frozenset(parts)


def stop_set_label(stop_set):
    """File-system safe label, e.g. for directory names: 'TAA-TAG', 'TGA'."""
    s = [c for c in CODONS if c in stop_set]
    return "-".join(s) if s else "NONE"
