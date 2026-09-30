"""Scribonia: fast estimation of the genetic-code stop set of genomes and bins."""

__version__ = "0.1.1"

from .tables import (CODONS, STOP_SET_TO_TABLE, TABLE_STOP_SETS,  # noqa: F401
                     canonical_table, format_stop_set, parse_stop_set,
                     stop_set_of_table)
