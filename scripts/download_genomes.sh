#!/bin/bash
# =============================================================================
# Download the labelled training genomes listed in data/manifest.tsv
# =============================================================================
# Usage (on a login node -- network I/O only, no compute):
#   bash scripts/download_genomes.sh [manifest.tsv] [genome_dir]
#
# Defaults: data/manifest.tsv and $SCRIBONIA_TRAIN_DIR/genomes
#
# Each verified row is fetched with the NCBI Datasets CLI into
#   <genome_dir>/<stop_set label>/<accession>/   (genome FASTA + GFF3 when any)
# where the label is TAA-TAG-TGA, TGA, TAA-TAG, TAA-TGA or TAG.  Rows whose
# `verified` column is not "yes" are skipped unless ALLOW_UNVERIFIED=1 (for
# evaluation-only genomes).  Where an NCBI annotation exists, transl_table=
# in the GFF3 is compared with the manifest's table and a mismatch is
# reported in the log (the manifest itself is not changed).
#
# Genomes are hundreds of MB to GB each: needs a fast, unmetered connection.
# =============================================================================

set -euo pipefail

MANIFEST="${1:-$(dirname "$0")/../data/manifest.tsv}"
GENOME_DIR="${2:-${SCRIBONIA_TRAIN_DIR:?give genome_dir or set SCRIBONIA_TRAIN_DIR}/genomes}"
ALLOW_UNVERIFIED="${ALLOW_UNVERIFIED:-0}"

command -v datasets >/dev/null 2>&1 || {
    echo "[ERROR] NCBI 'datasets' CLI not found (conda install -c conda-forge ncbi-datasets-cli)" >&2
    exit 1
}
command -v unzip >/dev/null 2>&1 || { echo "[ERROR] unzip not found" >&2; exit 1; }

mkdir -p "$GENOME_DIR"
LOG="$GENOME_DIR/download_$(date +%Y%m%d_%H%M%S).log"
echo "manifest: $MANIFEST" | tee "$LOG"

label_of() {
    # "TAA,TAG,TGA" -> "TAA-TAG-TGA"
    echo "$1" | tr ',' '-'
}

# Tabs become \037 first: a tab in IFS is whitespace, so read would merge
# empty fields (e.g. an empty citation) and shift every later column.
grep -v '^#' "$MANIFEST" | tr '\t' '\037' | while IFS=$'\037' read -r accession species taxid stop_set table ambiguous label_source citation verified note role; do
    [ -z "$accession" ] && continue
    if [ "$verified" != "yes" ] && [ "$ALLOW_UNVERIFIED" != "1" ]; then
        echo "SKIP (unverified) $accession $species" | tee -a "$LOG"
        continue
    fi
    label=$(label_of "$stop_set")
    target="$GENOME_DIR/$label/$accession"
    if compgen -G "$target/*.fna" >/dev/null || compgen -G "$target/*.fna.gz" >/dev/null; then
        echo "have $accession ($species)" | tee -a "$LOG"
        continue
    fi
    mkdir -p "$target"
    case "$accession" in
        GCA_*|GCF_*)
            zip="$target/ncbi_dataset.zip"
            if datasets download genome accession "$accession" --include genome,gff3 --filename "$zip" </dev/null >>"$LOG" 2>&1; then
                (cd "$target" && unzip -o -q ncbi_dataset.zip && rm -f ncbi_dataset.zip)
                find "$target" -name '*.fna' -exec mv -t "$target" {} + 2>/dev/null || true
                find "$target" -name '*.gff' -exec mv -t "$target" {} + 2>/dev/null || true
                gff=$(ls "$target"/*.gff 2>/dev/null | head -1 || true)
                if [ -n "$gff" ]; then
                    # Majority table over CDS features; a CDS without a
                    # transl_table tag is table 1 (NCBI tags organelle CDS
                    # only, so counting tags alone reports the mito code).
                    seen=$(awk -F'\t' '$3 == "CDS" { t = 1
                        if (match($9, /transl_table=[0-9]+/)) t = substr($9, RSTART + 13, RLENGTH - 13)
                        n[t]++ } END { for (t in n) if (n[t] > m) { m = n[t]; b = t }; print b }' "$gff")
                    if [ -n "$seen" ]; then
                        if [ "$seen" = "$table" ]; then
                            echo "OK   $accession ($species): annotation transl_table=$seen matches manifest" | tee -a "$LOG"
                        else
                            echo "WARN $accession ($species): annotation transl_table=$seen, manifest says $table" | tee -a "$LOG"
                        fi
                    fi
                fi
                echo "$species	$taxid	$stop_set	$table	$citation" > "$target/LABEL.tsv"
            else
                echo "FAIL $accession ($species): datasets download failed" | tee -a "$LOG"
                rmdir "$target" 2>/dev/null || true
            fi
            ;;
        *)
            echo "MANUAL $accession ($species): not a GCA/GCF accession; fetch by hand (WGS/BioProject) into $target" | tee -a "$LOG"
            rmdir "$target" 2>/dev/null || true
            ;;
    esac
done

echo "done; log: $LOG"
