#!/bin/bash
# =============================================================================
# Split training run: prep -> simulation array -> train/evaluate
# =============================================================================
# Usage (login node, Scribonia repo root; only submits jobs):
#   export SCRIBONIA_TRAIN_DIR=/path/to/train_dir
#   VERSION=scribonia_v3 bash scripts/submit_train.sh
# Set the partition with SBATCH_PARTITION (read by sbatch itself).
#
# Same outputs as `sbatch scripts/train.sbatch`, but the bins are
# simulated in N_CHUNKS_TRAIN + N_CHUNKS_HELDOUT array tasks (one node each)
# instead of on one node: ~30 min instead of ~3.5 h for 20000 + 3000 bins.
# Every variable train.sbatch reads (VERSION, N_TRAIN, CONTIG_DIST, ...)
# is passed through.
# =============================================================================

set -euo pipefail
cd "$(dirname "$0")/.."

export VERSION=${VERSION:-scribonia_v3}
export N_CHUNKS_TRAIN=${N_CHUNKS_TRAIN:-10}
export N_CHUNKS_HELDOUT=${N_CHUNKS_HELDOUT:-3}
S=scripts/train.sbatch
: "${SCRIBONIA_TRAIN_DIR:?set SCRIBONIA_TRAIN_DIR to the training working directory}"
export SCRIBONIA_TRAIN_DIR
LOG=${LOG_PREFIX:-scribonia_train_$VERSION}

prep=$(STAGE=prep sbatch --parsable --export=ALL --cpus-per-task=48 --mem=64G \
       --time=02:00:00 -o "${LOG}_prep_%j.out" -e "${LOG}_prep_%j.err" "$S")
sim=$(STAGE=sim sbatch --parsable --export=ALL --dependency=afterok:"$prep" \
      --array=0-$((N_CHUNKS_TRAIN + N_CHUNKS_HELDOUT - 1)) --cpus-per-task=48 \
      --mem=100G --time=04:00:00 -o "${LOG}_sim_%A_%a.out" -e "${LOG}_sim_%A_%a.err" "$S")
final=$(STAGE=final sbatch --parsable --export=ALL --dependency=afterok:"$sim" \
        --cpus-per-task=32 --mem=64G --time=03:00:00 \
        -o "${LOG}_final_%j.out" -e "${LOG}_final_%j.err" "$S")
echo "prep $prep  sim $sim  final $final  ($VERSION)"
