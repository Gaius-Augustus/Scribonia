#!/bin/bash
#SBATCH --job-name=gcbench
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=gcbench_%j.out
#SBATCH --error=gcbench_%j.err
# =============================================================================
# Scribonia vs Codetta vs PhyloFisher on eukaryotic genomes and bins
# =============================================================================
# Submit from the Scribonia repo root, after benchmark/setup.sh and
# benchmark/setup.sbatch (environment: see submit_bench.sh):
#   bash benchmark/submit_bench.sh      (chains all stages below)
#
# STAGE selects the part:
#   inputs       whole genomes + N_BINS simulated bins per genome as FASTA
#                (benchmark/make_inputs.py) -> $BENCH/samples.tsv, inputs/
#   scribonia    one array task per samples.tsv row (task i = row i), each
#   codetta      tool writes $BENCH/results/<tool>/<sample>/ with its output
#   phylofisher  and time.tsv (GNU time: wall, user, sys s, max RSS kB);
#                a finished sample (file "done") is skipped on resubmission
#   score        benchmark/score.py -> $BENCH/report/{summary,per_sample,pf_sweep}.tsv
#
# Samples: the held-out species of scripts/train.sbatch (never simulated
# into Scribonia's training set) minus EXCLUDE (default Danio: a 1.4 Gb
# animal genome, outside the protist scope, ~150 CPU h for Codetta).
# Bins follow the held-out evaluation: real contig lengths of metagenomic
# bins, 1-30 Mb, contamination in half of them.
#
# The score stage records what was run in $BENCH/versions.tsv.  If the
# working tree differs from HEAD, the difference is kept next to it
# (versions_dirty.txt, versions_dirty.diff); run from a clean checkout of a
# tagged commit for numbers that go into the docs.
#
# All tools get the same cores for wall time; compare CPU time for cost.
# Scribonia is single-threaded.  Codetta runs with its defaults through
# benchmark/run_codetta.py (hmmscan chunks in parallel, inference unchanged);
# PhyloFisher genetic_code_examiner through benchmark/run_phylofisher.py
# (its defaults: e-value 1e-30, conservation 0.7; counts kept).
# =============================================================================

set -euo pipefail
source ~/.bashrc 2>/dev/null || true

BENCH=${BENCH:?set BENCH to the benchmark working directory}
GENOMES=${GENOMES:-${SCRIBONIA_TRAIN_DIR:?set SCRIBONIA_TRAIN_DIR or GENOMES}/genomes}
MANIFEST=${MANIFEST:-data/manifest.tsv}
CONTIG_DIST=${CONTIG_DIST:-${SCRIBONIA_TRAIN_DIR:-}/contig_length_dist.tsv}
N_BINS=${N_BINS:-5}
SPECIES=${SPECIES:-$(sed -n "s/^HELDOUT='\(.*\)'$/\1/p" scripts/train.sbatch)}
EXCLUDE=${EXCLUDE-Danio}
GCB_ENV=${GCB_ENV:?set GCB_ENV to the conda env prefix for the benchmark tools}
CODETTA_DIR=${CODETTA_DIR:?set CODETTA_DIR to the Codetta checkout}
PF_DB=${PF_DB:-${PF_DIR:?set PF_DIR or PF_DB}/PhyloFisherDatabase_v1.0/database}
# genetic_code_examiner uses, per ortholog, the first of these taxa present
# in the database alignment as tblastn query.  A fixed, taxonomically
# broad list: metagenomic bins come without a reliable taxonomy.  Checked
# against the database by setup.sbatch.
PF_QUERIES=${PF_QUERIES:-Homosapi,Arabthal,Dictdisc,Naeggrub,Plasfalc,Tetrther}
SCRIBONIA_PY=${SCRIBONIA_PY:-python3}
THREADS=${SLURM_CPUS_PER_TASK:-8}
STAGE=${STAGE:?set STAGE to inputs, scribonia, codetta, phylofisher or score}

[ -n "$SPECIES" ] || { echo "[ERROR] no HELDOUT line in scripts/train.sbatch; set SPECIES" >&2; exit 1; }

case "$STAGE" in
inputs)
    # Bins must use the real contig lengths of the held-out evaluation;
    # CONTIG_DIST=none opts into make_inputs' log-uniform lengths instead.
    dist=()
    if [ "$CONTIG_DIST" != none ]; then
        [ -s "$CONTIG_DIST" ] || { echo "[ERROR] no contig lengths at '$CONTIG_DIST'; set CONTIG_DIST (or none)" >&2; exit 1; }
        dist=(--contig-dist "$CONTIG_DIST")
    fi
    "$SCRIBONIA_PY" -m benchmark.make_inputs "$MANIFEST" "$GENOMES" "$BENCH" \
        --species "$SPECIES" --exclude "$EXCLUDE" --n-bins "$N_BINS" ${dist[@]+"${dist[@]}"}
    ;;

scribonia|codetta|phylofisher)
    [ -x /usr/bin/time ] || { echo "[ERROR] GNU time (/usr/bin/time) not found" >&2; exit 1; }
    i=${SLURM_ARRAY_TASK_ID:?STAGE=$STAGE needs an array job}
    row=$(awk -F'\t' -v n=$((i + 2)) 'NR == n' "$BENCH/samples.tsv")
    [ -n "$row" ] || { echo "no sample for task $i"; exit 0; }
    sample=$(cut -f1 <<< "$row")
    fasta=$(cut -f10 <<< "$row")
    OUT=$BENCH/results/$STAGE/$sample
    [ -e "$OUT/done" ] && { echo "$STAGE $sample: done before"; exit 0; }
    mkdir -p "$OUT"
    export PATH=$GCB_ENV/bin:$PATH
    T=(/usr/bin/time -f "%e %U %S %M" -o "$OUT/time.tsv")
    echo "[$(date)] $STAGE $sample ($fasta)"
    case "$STAGE" in
    scribonia)
        "${T[@]}" "$SCRIBONIA_PY" -m scribonia.cli classify "$fasta" --bag "$sample" \
            -o "$OUT/out.tsv" -v ;;
    codetta)
        "${T[@]}" "$GCB_ENV/bin/python3" benchmark/run_codetta.py \
            --codetta-dir "$CODETTA_DIR" --threads "$THREADS" "$fasta" "$OUT/out" ;;
    phylofisher)
        "${T[@]}" "$GCB_ENV/bin/python3" benchmark/run_phylofisher.py --database "$PF_DB" \
            --queries "$PF_QUERIES" --threads "$THREADS" "$fasta" "$OUT" ;;
    esac
    touch "$OUT/done"
    echo "[$(date)] $STAGE $sample finished: $(cat "$OUT/time.tsv")"
    ;;

score)
    # Job logs (untracked) do not make the tree dirty; changed or new
    # sources do, and are saved so that the run can be reconstructed.
    rm -f "$BENCH/versions_dirty.txt" "$BENCH/versions_dirty.diff"
    dirty=$(git status --porcelain | grep -vE '^\?\? (logs/|.*\.(out|err)$)' || true)
    if [ -n "$dirty" ]; then
        echo "$dirty" > "$BENCH/versions_dirty.txt"
        git diff HEAD > "$BENCH/versions_dirty.diff"
        echo "[WARN] working tree differs from HEAD, see $BENCH/versions_dirty.txt" >&2
    fi
    # The conda package version: PhyloFisher 1.2.14 reports itself as 1.2.13.
    pf_version=$(ls "$GCB_ENV"/conda-meta/phylofisher-*.json 2>/dev/null | sed -n '1s/.*phylofisher-\([^-]*\)-.*/\1/p')
    {
        echo "scribonia	$(git rev-parse HEAD)${dirty:+ dirty}"
        echo "scribonia_version	$("$SCRIBONIA_PY" -c 'import scribonia; print(scribonia.__version__)')"
        echo "codetta	$(git -C "$CODETTA_DIR" rev-parse HEAD)"
        echo "phylofisher	${pf_version:-$("$GCB_ENV/bin/python3" -c 'import pkg_resources as p; print(p.get_distribution("phylofisher").version)' 2>/dev/null || echo unknown)}"
        echo "species	$SPECIES"
        echo "exclude	$EXCLUDE"
        echo "pf_queries	$PF_QUERIES"
    } > "$BENCH/versions.tsv"
    "$SCRIBONIA_PY" -m benchmark.score "$BENCH/samples.tsv" "$BENCH/results" "$BENCH/report"
    ;;
*)
    echo "[ERROR] STAGE=$STAGE" >&2; exit 1 ;;
esac
