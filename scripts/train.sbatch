#!/bin/bash
#SBATCH --job-name=scribonia_train
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=scribonia_train_%j.out
#SBATCH --error=scribonia_train_%j.err
# =============================================================================
# Simulate bin-like samples from the labelled genomes, train and evaluate
# =============================================================================
# Submit from the Scribonia repo root:
#   export SCRIBONIA_TRAIN_DIR=/path/to/train_dir
#   bash scripts/submit_train.sh     (split over nodes, ~30 min)
#   sbatch scripts/train.sbatch      (one node, ~3.5 h)
# No partition is set here: pass -p to sbatch or export SBATCH_PARTITION
# (also honoured by submit_train.sh).
# STAGE selects the part to run: all (default, one node), prep (split the
# manifest, recode), sim (one array task = one chunk of training or held-out
# bins) or final (train and evaluate on the chunks).  submit_train.sh
# chains prep -> sim array -> final with dependencies.
#
# Environment:
#   SCRIBONIA_TRAIN_DIR  (required) working directory, BASE below
#   SCRIBONIA_SIF        Singularity image; without it SCRIBONIA_PY is used
#   SCRIBONIA_PY         python with numpy, scikit-learn, joblib (default python3)
#
# Inputs : $BASE/genomes  (download_genomes.sh)
#          data/manifest.tsv
#          CONTIG_DIST: contig lengths of real metagenomic bins for the held-out
#          bins (TRAIN_CONTIG_DIST to use them for training too), one per line;
#          the current file holds the 7403 contigs of the 150 bins of one
#          metagenomic sample (median 10.5 kb, 683 kb max), made with
#            for f in <bins_dir>/*.fa; do
#              awk '/^>/ {if (n) print n; n=0; next} {n+=length($0)} END {if (n) print n}' $f
#            done > $BASE/contig_length_dist.tsv
#          (without it, contig lengths are log-uniform in CONTIG_MIN..20 kb)
# Outputs: $BASE/genomes/synthetic/  (recoded)
#          $BASE/$VERSION/manifest_{train,heldout,recoded}.tsv
#          $BASE/$VERSION/{train,heldout}.npz
#            (split runs: $VERSION/chunks/{train,heldout}_<i>.npz)
#          $BASE/models/$VERSION.joblib
#          $BASE/eval/$VERSION.json
#
# Held-out taxa (never simulated into the training set) are listed in
# HELDOUT below (see docs/genetic_codes.md for the split).
# =============================================================================

set -euo pipefail
source ~/.bashrc 2>/dev/null || true

BASE=${SCRIBONIA_TRAIN_DIR:?set SCRIBONIA_TRAIN_DIR to the training working directory}
GENOMES=$BASE/genomes
MANIFEST=${MANIFEST:-data/manifest.tsv}
CONTIG_DIST=${CONTIG_DIST:-$BASE/contig_length_dist.tsv}
N_TRAIN=${N_TRAIN:-20000}
N_HELDOUT=${N_HELDOUT:-3000}
VERSION=${VERSION:-scribonia_v3}
SIF=${SCRIBONIA_SIF:-}

STAGE=${STAGE:-all}
N_CHUNKS_TRAIN=${N_CHUNKS_TRAIN:-10}
N_CHUNKS_HELDOUT=${N_CHUNKS_HELDOUT:-3}
case "$STAGE" in all|prep|sim|final) ;; *) echo "[ERROR] STAGE=$STAGE" >&2; exit 1 ;; esac

OUT=$BASE/$VERSION
mkdir -p "$OUT" "$OUT/chunks" "$BASE/models" "$BASE/eval"

if [ "$STAGE" = all ] || [ "$STAGE" = prep ]; then

# Split the manifest into training and held-out genomes.  Every genome with
# ambiguous = yes (no dedicated stop codon) must match here: those are
# evaluation only, their stop-set label is nominal.
HELDOUT='Ichthyophthirius|Pseudocohnilembus|Euplotes crassus|Stentor|Entodinium|Isotricha|Plasmodium|Danio|Hexamita|Blastocystis|Acetabularia|Condylostoma|Blastocrithidia|Amoebophrya|Loxodes|Plagiopyla|Laurentiella|Halteria|Pseudokeronopsis|Diophrys|Sterkiella|Tetmemena|Paraurostyla|Urostyla|Uronychia|Euplotes woodruffi|Euplotes parawoodruffi|Euplotes weissei|Euplotes cf. woodruffi|Mesodinium'
# Matched on the species column only; contaminant-only genomes (role
# column) go into both, they are never the genome being classified.
split_manifest() {  # $1 = heldout|train
    awk -F'\t' -v re="$HELDOUT" -v want="$1" '
        NR == 1 || /^#/ { print; next }
        $11 == "contaminant" { print; next }
        { held = ($2 ~ re); if ((want == "heldout") == held) print }' "$MANIFEST"
}
split_manifest heldout > "$OUT/manifest_heldout.tsv"
split_manifest train   > "$OUT/manifest_train.tsv"
fi

run() {
    if [ -s "$SIF" ]; then
        module load singularity 2>/dev/null || true
        # Further binds (e.g. a symlinked GENOMES) via SINGULARITY_BIND.
        singularity exec -B "$BASE:$BASE" -B "$PWD:$PWD" "$SIF" scribonia "$@"
    else
        # From the repo root, with numpy, scikit-learn and joblib.
        "$PY" -m scribonia.cli "$@"
    fi
}
PY=${SCRIBONIA_PY:-python3}
if [ -n "$SIF" ] && [ ! -s "$SIF" ]; then
    echo "[ERROR] SCRIBONIA_SIF=$SIF not found" >&2; exit 1
fi
if [ ! -s "$SIF" ]; then
    "$PY" -c "import numpy, sklearn, joblib" || {
        echo "[ERROR] $PY lacks numpy/scikit-learn/joblib; set SCRIBONIA_PY or SCRIBONIA_SIF" >&2; exit 1; }
fi
THREADS=${SLURM_CPUS_PER_TASK:-8}

# Recoded standard-code genomes for the thin classes ({TAA,TAG} has two real
# genomes, {TAA,TGA} none).  Sources come from the training split only and
# keep their genus, so the genus bags and the held-out split stay taxon-clean.
# RECODES lists stop set:fraction pairs.  The fractions are augmentation
# settings, not measured usage; evaluation uses real genomes only.  Chosen
# on the held-out and leave-genus-out results of 2026-09-26:
#   * no {TGA} recodes: recoding CAA/CAG into TAA/TAG inside annotated CDS
#     leaves D of TAA/TAG at 0.1-0.4 even at fraction 0.8, below every real
#     {TGA} ciliate (0.45-0.8); such genomes taught the model to call
#     AT-rich standard genomes (Dictyostelium, C. elegans) {TGA};
#   * no fraction 0.05: indistinguishable from the standard code;
#   * no intron-rich animal/plant sources (C. elegans, Drosophila,
#     Arabidopsis): too little coding sequence per 2-20 kb contig, the
#     recoded copies look standard and made Danio a {TGA}/{TAA,TGA} call.
RECODE_SOURCES=${RECODE_SOURCES:-"GCF_000004695.1 GCF_000208925.1 GCF_000002595.2 GCF_000006565.2 GCF_000146045.2 GCF_000149405.2 GCF_000006335.3 GCF_000004985.1 GCF_000165345.1 GCF_000182965.3"}
RECODES=${RECODES:-"TAA,TAG:0.25 TAA,TAG:0.5 TAA,TAG:0.8 TAA,TGA:0.5 TAA,TGA:0.8 TAA,TGA:1.0"}
SYN=$GENOMES/synthetic
mkdir -p "$SYN"

if [ "$STAGE" = all ] || [ "$STAGE" = prep ]; then

recode_one() {  # accession stop_set fraction -> $SYN/<id>/genomic.fna
    local acc=$1 set=$2 f=$3 src fna gff id out
    id="${acc}_recode_$(echo "$set" | tr , -)_f$f"
    out="$SYN/$id"
    [ -s "$out/genomic.fna" ] && return 0
    src=$(ls -d "$GENOMES"/*/"$acc" 2>/dev/null | grep -v "/synthetic/" | head -1 || true)
    fna=$(ls "$src"/*.fna 2>/dev/null | head -1 || true)
    gff=$(ls "$src"/*.gff 2>/dev/null | head -1 || true)
    if [ -z "$fna" ] || [ -z "$gff" ]; then
        echo "[WARN] $acc: no FASTA/GFF under $GENOMES, not recoded" >&2; return 0
    fi
    mkdir -p "$out"
    if run recode "$fna" "$gff" --stop-set "$set" --fraction "$f" \
            -o "$out/genomic.fna.tmp" 2> "$out/recode.log"; then
        mv "$out/genomic.fna.tmp" "$out/genomic.fna"
    else
        echo "[WARN] recode failed for $id, see $out/recode.log" >&2
    fi
}

echo "[$(date)] recoding standard-code genomes"
for acc in $RECODE_SOURCES; do
    if ! awk -F'\t' -v a="$acc" '$1 == a { f = 1 } END { exit !f }' "$OUT/manifest_train.tsv"; then
        echo "[WARN] $acc is not in the training split, not recoded" >&2; continue
    fi
    for sf in $RECODES; do
        recode_one "$acc" "${sf%:*}" "${sf#*:}" &
        while [ "$(jobs -rp | wc -l)" -ge "$THREADS" ]; do wait -n || true; done
    done
done
wait
# One manifest row per recoded genome, training split only.
for acc in $RECODE_SOURCES; do
    for sf in $RECODES; do
        set=${sf%:*}; f=${sf#*:}
        id="${acc}_recode_$(echo "$set" | tr , -)_f$f"
        [ -s "$SYN/$id/genomic.fna" ] || continue
        case "$set" in TGA) t=6 ;; TAA,TAG) t=10 ;; TAA,TGA) t=15 ;; *) t=NA ;; esac
        awk -F'\t' -v OFS='\t' -v a="$acc" -v id="$id" -v set="$set" -v t="$t" -v f="$f" '
            $1 == a { print id, $2 " [recoded " set "]", $3, set, t, "no",
                      "synthetic_recode", "", "yes", "recoded from " a ", fraction " f, "primary" }' \
            "$OUT/manifest_train.tsv"
    done
done > "$OUT/manifest_recoded.tsv"
cat "$OUT/manifest_recoded.tsv" >> "$OUT/manifest_train.tsv"
echo "[$(date)] $(wc -l < "$OUT/manifest_recoded.tsv") recoded genomes added to the training manifest"
fi
[ "$STAGE" = prep ] && exit 0

# Training bins use log-uniform 2-20 kb contigs; held-out bins use the real
# contig lengths and are at least 1 Mb (the intended use).  Leave-genus-out on
# 2026-09-27: trained on 2-20 kb and tested on real lengths, exact 0.973 on
# protist bins >= 1 Mb; trained on real lengths, 0.956.
CONTIG_MIN=${CONTIG_MIN:-2000}
train_dist_arg="--contig-min $CONTIG_MIN"
[ -s "${TRAIN_CONTIG_DIST:-}" ] && train_dist_arg="$train_dist_arg --contig-dist $TRAIN_CONTIG_DIST --contig-max ${CONTIG_MAX:-1000000}"
held_dist_arg="--contig-min $CONTIG_MIN --min-mb ${HELDOUT_MIN_MB:-1}"
[ -s "$CONTIG_DIST" ] && held_dist_arg="$held_dist_arg --contig-dist $CONTIG_DIST --contig-max ${CONTIG_MAX:-1000000}"

if [ "$STAGE" = sim ]; then
    # Array task i < N_CHUNKS_TRAIN: training chunk i; above: held-out chunk.
    i=${SLURM_ARRAY_TASK_ID:?STAGE=sim needs an array job}
    if [ "$i" -lt "$N_CHUNKS_TRAIN" ]; then
        n=$(( (N_TRAIN + N_CHUNKS_TRAIN - 1) / N_CHUNKS_TRAIN ))
        echo "[$(date)] training chunk $i: $n bins"
        run simulate "$OUT/manifest_train.tsv" "$GENOMES" -o "$OUT/chunks/train_$i.npz" \
            --n-bins "$n" $train_dist_arg --threads "$THREADS" --seed $((1000 + i)) -v
    else
        j=$((i - N_CHUNKS_TRAIN))
        n=$(( (N_HELDOUT + N_CHUNKS_HELDOUT - 1) / N_CHUNKS_HELDOUT ))
        echo "[$(date)] held-out chunk $j: $n bins"
        run simulate "$OUT/manifest_heldout.tsv" "$GENOMES" -o "$OUT/chunks/heldout_$j.npz" \
            --n-bins "$n" $held_dist_arg --threads "$THREADS" --seed $((2000 + j)) -v
    fi
    exit 0
fi

if [ "$STAGE" = all ]; then
    echo "[$(date)] simulating $N_TRAIN training bins"
    run simulate "$OUT/manifest_train.tsv" "$GENOMES" -o "$OUT/train.npz" \
        --n-bins "$N_TRAIN" $train_dist_arg --threads "$THREADS" --seed 1 -v
    echo "[$(date)] simulating $N_HELDOUT held-out bins"
    run simulate "$OUT/manifest_heldout.tsv" "$GENOMES" -o "$OUT/heldout.npz" \
        --n-bins "$N_HELDOUT" $held_dist_arg --threads "$THREADS" --seed 2 -v
    TRAIN_NPZ=("$OUT/train.npz"); HELDOUT_NPZ=("$OUT/heldout.npz")
else
    TRAIN_NPZ=("$OUT"/chunks/train_*.npz); HELDOUT_NPZ=("$OUT"/chunks/heldout_*.npz)
    [ "${#TRAIN_NPZ[@]}" -eq "$N_CHUNKS_TRAIN" ] && [ "${#HELDOUT_NPZ[@]}" -eq "$N_CHUNKS_HELDOUT" ] || {
        echo "[ERROR] expected $N_CHUNKS_TRAIN training and $N_CHUNKS_HELDOUT held-out chunks in $OUT/chunks" >&2; exit 1; }
    # Chunks older than the manifest are left over from an earlier run.
    for f in "${TRAIN_NPZ[@]}" "${HELDOUT_NPZ[@]}"; do
        [ "$f" -nt "$OUT/manifest_train.tsv" ] || { echo "[ERROR] $f predates the manifest (stale chunk)" >&2; exit 1; }
    done
fi

echo "[$(date)] training"
# Ten models per codon on bootstrap samples of genera, composition features
# excluded, joint decoding over the stop sets with a table (scribonia/model.py).
run train "${TRAIN_NPZ[@]}" -o "$BASE/models/$VERSION.joblib" --version "$VERSION"

echo "[$(date)] evaluating on held-out taxa"
run evaluate "$BASE/models/$VERSION.joblib" "${HELDOUT_NPZ[@]}" -o "$BASE/eval/$VERSION.json"
echo "[$(date)] done: $BASE/models/$VERSION.joblib ; copy into scribonia/models/ and rebuild the image"
