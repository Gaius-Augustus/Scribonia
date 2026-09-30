#!/bin/bash
# =============================================================================
# Submit the benchmark: inputs -> three tool arrays -> score
# =============================================================================
# Usage (login node, Scribonia repo root; only submits jobs):
#   export BENCH=/path/to/bench_dir GCB_ENV=... CODETTA_DIR=... PF_DIR=...
#   export SCRIBONIA_TRAIN_DIR=...   (or GENOMES and CONTIG_DIST)
#   bash benchmark/submit_bench.sh
# Set the partition with SBATCH_PARTITION (read by sbatch itself).
#
# Every variable benchmark/bench.sbatch reads (BENCH, N_BINS, SPECIES,
# EXCLUDE, PF_QUERIES, ...) is passed through.  The array size is the
# number of samples make_inputs.py will write: matching genomes x
# (1 + N_BINS).  MAX_PARALLEL caps running tasks per tool (default 20).
# =============================================================================

set -euo pipefail
cd "$(dirname "$0")/.."

export BENCH=${BENCH:?set BENCH to the benchmark working directory}
for v in GCB_ENV CODETTA_DIR; do
    [ -n "${!v:-}" ] || { echo "[ERROR] set $v (see benchmark/setup.sh)" >&2; exit 1; }
done
[ -n "${PF_DB:-}${PF_DIR:-}" ] || { echo "[ERROR] set PF_DIR or PF_DB" >&2; exit 1; }
[ -n "${GENOMES:-}${SCRIBONIA_TRAIN_DIR:-}" ] || { echo "[ERROR] set SCRIBONIA_TRAIN_DIR or GENOMES" >&2; exit 1; }
export N_BINS=${N_BINS:-5}
export SPECIES=${SPECIES:-$(sed -n "s/^HELDOUT='\(.*\)'$/\1/p" scripts/train.sbatch)}
export EXCLUDE=${EXCLUDE-Danio}
S=benchmark/bench.sbatch
LOG=$BENCH/logs
mkdir -p "$LOG"

# Same filter as make_inputs.py: verified primary rows, no recodes.
n_genomes=$(awk -F'\t' -v re="$SPECIES" -v ex="$EXCLUDE" '
    /^#/ || $1 == "" { next }
    $9 == "yes" && ($11 == "primary" || $11 == "") && $7 != "synthetic_recode" \
        && $2 ~ re && !(ex != "" && $2 ~ ex)' data/manifest.tsv | wc -l)
n=$(( n_genomes * (1 + N_BINS) ))
[ "$n" -gt 0 ] || { echo "[ERROR] no genome matches SPECIES" >&2; exit 1; }

inputs=$(STAGE=inputs sbatch --parsable --export=ALL --cpus-per-task=2 --mem=64G \
         --time=04:00:00 -o "$LOG/inputs_%j.out" -e "$LOG/inputs_%j.err" "$S")
tools=()
for t in scribonia codetta phylofisher; do
    cpus=16; [ "$t" = scribonia ] && cpus=2
    tools+=("$(STAGE=$t sbatch --parsable --export=ALL --dependency=afterok:"$inputs" \
               --array=0-$((n - 1))%"${MAX_PARALLEL:-20}" --cpus-per-task=$cpus \
               -o "$LOG/${t}_%A_%a.out" -e "$LOG/${t}_%A_%a.err" "$S")")
done
score=$(STAGE=score sbatch --parsable --export=ALL \
        --dependency=afterany:"$(IFS=:; echo "${tools[*]}")" --cpus-per-task=1 --mem=4G \
        --time=00:30:00 -o "$LOG/score_%j.out" -e "$LOG/score_%j.err" "$S")
echo "$n samples ($n_genomes genomes); inputs $inputs  scribonia/codetta/phylofisher ${tools[*]}  score $score"
