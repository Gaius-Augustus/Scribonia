#!/bin/bash
# =============================================================================
# Benchmark setup, part 1 of 2: downloads only (login node)
# =============================================================================
# Usage (Scribonia repo root):
#   export GCB_ENV=... CODETTA_DIR=... PF_DIR=...
#   bash benchmark/setup.sh
#   sbatch benchmark/setup.sbatch      (part 2: unpack, compile, align)
#
# Fetches (~5 GB):
#   * conda env $GCB_ENV: PhyloFisher $PF_VERSION (bioconda; brings BLAST+,
#     MAFFT, Biopython, Python 3.7) plus numpy/scipy for Codetta;
#   * Codetta at commit $CODETTA_COMMIT into $CODETTA_DIR, and the Pfam 35
#     profile database Codetta ships (Pfam-A_enone.tar.gz, eddylab.org);
#   * the PhyloFisher database v1.0 (Tice et al., figshare file 29093409).
# Nothing is unpacked or compiled here.
# =============================================================================

set -euo pipefail

GCB_ENV=${GCB_ENV:?set GCB_ENV to the conda env prefix for the benchmark tools}
CODETTA_DIR=${CODETTA_DIR:?set CODETTA_DIR to the Codetta checkout}
PF_DIR=${PF_DIR:?set PF_DIR to the PhyloFisher database directory}
CODETTA_COMMIT=863359ed326276602d44e48227b6003ac6ffd266   # master, 2025-05-05 (v2.0)
PF_VERSION=1.2.14
PFAM_URL=http://eddylab.org/publications/Shulgina21/Pfam-A_enone.tar.gz
PF_DB_URL=https://ndownloader.figshare.com/files/29093409
PF_DB_TAR=Tice_etal.PhyloFisherDatabase_v1.0_Apr.11.2021.tar.gz

if [ ! -x "$GCB_ENV/bin/python3" ]; then
    CONDA=$(command -v mamba || command -v conda) || { echo "[ERROR] no mamba/conda" >&2; exit 1; }
    "$CONDA" create -y -p "$GCB_ENV" -c conda-forge -c bioconda \
        "phylofisher=$PF_VERSION" numpy scipy wget
fi

if [ ! -d "$CODETTA_DIR/.git" ]; then
    git clone https://github.com/kshulgina/codetta.git "$CODETTA_DIR"
fi
git -C "$CODETTA_DIR" checkout -q "$CODETTA_COMMIT"

if [ ! -s "$CODETTA_DIR/resources/Pfam-A_enone.hmm" ] && [ ! -s "$CODETTA_DIR/resources/Pfam-A_enone.tar.gz" ]; then
    wget -q -O "$CODETTA_DIR/resources/Pfam-A_enone.tar.gz.part" "$PFAM_URL"
    mv "$CODETTA_DIR/resources/Pfam-A_enone.tar.gz.part" "$CODETTA_DIR/resources/Pfam-A_enone.tar.gz"
fi

mkdir -p "$PF_DIR"
if [ ! -d "$PF_DIR/PhyloFisherDatabase_v1.0" ] && [ ! -s "$PF_DIR/$PF_DB_TAR" ]; then
    wget -q -O "$PF_DIR/$PF_DB_TAR.part" "$PF_DB_URL"
    mv "$PF_DIR/$PF_DB_TAR.part" "$PF_DIR/$PF_DB_TAR"
fi

echo "downloads done; next: sbatch benchmark/setup.sbatch"
