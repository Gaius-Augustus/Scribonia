#!/bin/bash
#SBATCH --job-name=gcbench_setup
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=gcbench_setup_%j.out
#SBATCH --error=gcbench_setup_%j.err
# =============================================================================
# Benchmark setup, part 2 of 2: unpack, compile, align, smoke test
# =============================================================================
# Run after benchmark/setup.sh, from the Scribonia repo root, with the same
# GCB_ENV, CODETTA_DIR and PF_DIR (partition: sbatch -p or SBATCH_PARTITION):
#   sbatch benchmark/setup.sbatch
#
#   * unpacks Pfam-A_enone and compiles Codetta's bundled HMMER 3.3.2
#     (codetta/setup.sh); presses the database if the tarball did not;
#   * runs Codetta on its bundled example, which also writes the Pfam
#     emissions cache (resources/*.emissions_dict.p) once, before array
#     tasks would race to write it;
#   * unpacks the PhyloFisher database and builds its ortholog alignments
#     (MAFFT; genetic_code_examiner --prepare_alignments);
#   * lists the database taxa and checks the PhyloFisher query names
#     (PF_QUERIES, as in bench.sbatch).
# =============================================================================

set -euo pipefail
source ~/.bashrc 2>/dev/null || true

GCB_ENV=${GCB_ENV:?set GCB_ENV to the conda env prefix for the benchmark tools}
CODETTA_DIR=${CODETTA_DIR:?set CODETTA_DIR to the Codetta checkout}
PF_DIR=${PF_DIR:?set PF_DIR to the PhyloFisher database directory}
PF_DB=$PF_DIR/PhyloFisherDatabase_v1.0/database
PF_DB_TAR=Tice_etal.PhyloFisherDatabase_v1.0_Apr.11.2021.tar.gz
THREADS=${SLURM_CPUS_PER_TASK:-8}
REPO=$PWD
export PATH=$GCB_ENV/bin:$PATH

echo "[$(date)] Codetta: HMMER and Pfam"
cd "$CODETTA_DIR"
[ -x hmmer-3.3.2/bin/hmmscan ] || bash setup.sh
[ -x hmmer-3.3.2/bin/hmmscan ] || { echo "[ERROR] HMMER build failed" >&2; exit 1; }
cd resources
if [ -s Pfam-A_enone.tar.gz ] && [ ! -s Pfam-A_enone.hmm ]; then
    tar xf Pfam-A_enone.tar.gz && rm Pfam-A_enone.tar.gz
fi
[ -s Pfam-A_enone.hmm.h3m ] || ../hmmer-3.3.2/bin/hmmpress Pfam-A_enone.hmm
cd "$REPO"

# Codetta's README example: Nasuia deltocephalinicola, UGA = Trp.  The
# stop-to-sense call on TGA (15th codon) is what the benchmark relies on.
echo "[$(date)] Codetta smoke test (writes the emissions cache)"
SMOKE=$(mktemp -d "${TMPDIR:-/tmp}/codetta_smoke.XXXXXX")
cp "$CODETTA_DIR/examples/GCA_014211875.1.fna" "$SMOKE/"
python3 benchmark/run_codetta.py --codetta-dir "$CODETTA_DIR" --threads "$THREADS" \
    "$SMOKE/GCA_014211875.1.fna" "$SMOKE/GCA_014211875.1"
code=$(tail -1 "$SMOKE/GCA_014211875.1.genetic_code.out")
echo "inferred $code"
echo "README   FFLLSSSSYY??CCWWL?LLPPPPHHQQ????I?IMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
[ "${code:14:1}" = W ] || { echo "[ERROR] Codetta smoke test: TGA not inferred as W" >&2; exit 1; }

echo "[$(date)] PhyloFisher database"
if [ ! -d "$PF_DB" ]; then
    tar xzf "$PF_DIR/$PF_DB_TAR" -C "$PF_DIR"
fi
[ -d "$PF_DB/orthologs" ] || { echo "[ERROR] no $PF_DB/orthologs after unpacking" >&2; exit 1; }
n_orth=$(ls "$PF_DB"/orthologs/*.fas | wc -l)
n_aln=$(ls "$PF_DB"/alignments/*.fas.aln 2>/dev/null | wc -l)
echo "orthologs $n_orth, alignments $n_aln"
if [ "$n_aln" -ne "$n_orth" ]; then
    python3 benchmark/run_phylofisher.py --database "$PF_DB" --threads "$THREADS" --prepare-alignments
fi

echo "[$(date)] PhyloFisher taxa (short name, long name, higher taxonomy): $PF_DIR/taxa.tsv"
python3 benchmark/run_phylofisher.py --database "$PF_DB" --list-taxa > "$PF_DIR/taxa.tsv"
wc -l < "$PF_DIR/taxa.tsv"
PF_QUERIES=${PF_QUERIES:-$(sed -n "s/^PF_QUERIES=\${PF_QUERIES:-\(.*\)}$/\1/p" benchmark/bench.sbatch)}
for q in ${PF_QUERIES//,/ }; do
    if cut -f1 "$PF_DIR/taxa.tsv" | grep -qx "$q"; then echo "query ok      $q"; else echo "query MISSING $q"; missing=1; fi
done
[ -z "${missing:-}" ] || { echo "[ERROR] PhyloFisher query names missing from the database" >&2; exit 1; }
echo "[$(date)] setup done"
