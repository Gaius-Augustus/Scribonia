#!/usr/bin/env python3
"""PhyloFisher genetic_code_examiner on one FASTA, with the counts kept.

genetic_code_examiner.py (PhyloFisher 1.2.14) writes only a PDF of bar
plots and prints "Suspicious codon" for every codon whose most frequent
conserved amino acid differs from the standard code; for a stop codon a
single hit is enough.  It also breaks on input paths with a directory or
more than one dot.  This wrapper calls its own functions (collect_queries,
blastDB_tblastn, collect_counts: tblastn of database orthologs, e-value
1e-30, positions conserved in > 70 % of the database alignment) and writes
the counts instead:

  <outdir>/counts.tsv   codon  aa  n   (one row per codon and conserved aa)

The stop-set call from these counts is made in benchmark/score.py.

Needs config-free access to the database: --database is the folder holding
orthologs/, alignments/ (shipped with database v1.0; otherwise made once
with --prepare-alignments, which runs MAFFT and writes into the database
folder) and metadata.tsv (the taxa; --list-taxa prints them).

Usage: python3 run_phylofisher.py --database DB --queries A,B,... --threads N FASTA OUTDIR
"""

import argparse
import csv
import os
import sys


def read_taxa(database):
    """Rows of the database's metadata.tsv (database v1.0 keeps the taxa
    there; newer PhyloFisher versions use a phylofisher.db instead)."""
    with open(os.path.join(database, "metadata.tsv")) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("fasta", nargs="?")
    p.add_argument("outdir", nargs="?")
    p.add_argument("--database", required=True, help="PhyloFisherDatabase_v1.0/database")
    p.add_argument("--queries", default="",
                   help="comma-separated short names; the first present per gene is used")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--conserved", type=float, default=0.7)
    p.add_argument("--evalue", default="1e-30")
    p.add_argument("--prepare-alignments", action="store_true",
                   help="build the ortholog alignments (once, setup) and exit")
    p.add_argument("--list-taxa", action="store_true",
                   help="print short name, long name, higher taxonomy of every taxon and exit")
    a = p.parse_args(argv)

    if a.list_taxa:
        for r in read_taxa(a.database):
            print(f"{r['Unique ID']}\t{r['Long Name']}\t{r['Higher Taxonomy']}")
        return

    from phylofisher.utilities import genetic_code_examiner as gce
    database = os.path.abspath(a.database)
    gce.dfo = database   # module global, set from config.ini in the original
    if a.prepare_alignments:
        gce.prepare_alignments(a.threads)
        return
    if not (a.fasta and a.outdir):
        p.error("FASTA and OUTDIR are required")
    if not gce.check_alignments():
        sys.exit("[phylofisher] no alignments; run with --prepare-alignments first")

    queries = [q for q in a.queries.split(",") if q]
    known = {r["Unique ID"] for r in read_taxa(database)}
    missing = [q for q in queries if q not in known]
    queries = [q for q in queries if q in known]
    if missing:
        print(f"[phylofisher] not in the database, ignored: {','.join(missing)}", file=sys.stderr)
    if not queries:
        sys.exit("[phylofisher] no valid query short name")

    fasta = os.path.abspath(a.fasta)
    os.makedirs(a.outdir, exist_ok=True)
    os.chdir(a.outdir)
    # makeblastdb writes its index next to the input: use a local link.
    if os.path.lexists("genome.fa"):
        os.remove("genome.fa")
    os.symlink(fasta, "genome.fa")
    gce.collect_queries(queries)
    gce.blastDB_tblastn("genome.fa", a.threads, a.evalue)
    stats = gce.collect_counts("genome.fa", a.conserved)
    with open("counts.tsv", "w") as fh:
        fh.write("codon\taa\tn\n")
        for codon in sorted(stats):
            aas = stats[codon]
            for aa in sorted(set(aas)):
                fh.write(f"{codon}\t{aa}\t{aas.count(aa)}\n")
    for f in os.listdir("."):
        if f.startswith("genome.fa.n") or f == "blast_out.xml":
            os.remove(f)


if __name__ == "__main__":
    main()
