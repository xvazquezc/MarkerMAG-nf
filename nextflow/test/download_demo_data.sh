#!/usr/bin/env bash
# Download the MarkerMAG demo dataset from Zenodo and prepare a small
# subsampled test set suitable for a quick local pipeline run.
#
# Usage (from repo root or nextflow/):
#   bash nextflow/test/download_demo_data.sh
#
# Requires: curl or wget, unzip, seqtk (in the markermag-nf Conda environment)
# Output:   nextflow/test/data/{demo_16S.fasta, mags/, test_R1.fasta, test_R2.fasta}
#           nextflow/test/samplesheet.csv (overwritten)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
ZIP_URL="https://zenodo.org/records/6466784/files/MarkerMAG_demo_data.zip"
ZIP_PATH="${DATA_DIR}/MarkerMAG_demo_data.zip"
DEMO_DIR="${DATA_DIR}/MarkerMAG_demo_data"

# Number of read PAIRS to subsample (keeps the test fast)
N_READS=5000

mkdir -p "${DATA_DIR}"

# ── Download ───────────────────────────────────────────────────────────────────
if [[ ! -f "${ZIP_PATH}" ]]; then
    echo "[1/4] Downloading demo dataset from Zenodo (~78 MB)..."
    if command -v curl &>/dev/null; then
        curl -L --progress-bar -o "${ZIP_PATH}" "${ZIP_URL}"
    else
        wget -q --show-progress -O "${ZIP_PATH}" "${ZIP_URL}"
    fi
else
    echo "[1/4] Archive already downloaded — skipping."
fi

# ── Extract ────────────────────────────────────────────────────────────────────
if [[ ! -d "${DEMO_DIR}" ]]; then
    echo "[2/4] Extracting..."
    unzip -q "${ZIP_PATH}" -d "${DATA_DIR}"
fi

# ── Copy 16S and MAGs ──────────────────────────────────────────────────────────
echo "[3/4] Staging 16S and MAG files..."
mkdir -p "${DATA_DIR}/mags"
cp "${DEMO_DIR}/demo_16S.fasta" "${DATA_DIR}/16s/demo_16S.fasta" 2>/dev/null || \
    { mkdir -p "${DATA_DIR}/16s" && cp "${DEMO_DIR}/demo_16S.fasta" "${DATA_DIR}/16s/demo_16S.fasta"; }
cp "${DEMO_DIR}/demo_MAGs/"*.fa "${DATA_DIR}/mags/"

# ── Subsample reads ────────────────────────────────────────────────────────────
echo "[4/4] Subsampling to ${N_READS} read pairs with seqtk..."
seqtk sample -s 42 "${DEMO_DIR}/demo_R1.fasta" "${N_READS}" > "${DATA_DIR}/test_R1.fasta"
seqtk sample -s 42 "${DEMO_DIR}/demo_R2.fasta" "${N_READS}" > "${DATA_DIR}/test_R2.fasta"

# ── Write samplesheet ──────────────────────────────────────────────────────────
cat > "${SCRIPT_DIR}/samplesheet.csv" << EOF
sample,r1,r2,mag_dir,mag_ext,16s_reads,16s_fasta
demo,${DATA_DIR}/test_R1.fasta,${DATA_DIR}/test_R2.fasta,${DATA_DIR}/mags,fa,,${DATA_DIR}/16s/demo_16S.fasta
EOF

echo ""
echo "Done. Test data ready in ${DATA_DIR}/"
echo "  16S:  ${DATA_DIR}/16s/demo_16S.fasta"
echo "  MAGs: ${DATA_DIR}/mags/"
echo "  R1:   ${DATA_DIR}/test_R1.fasta  (${N_READS} reads)"
echo "  R2:   ${DATA_DIR}/test_R2.fasta  (${N_READS} reads)"
echo ""
echo "Run the test with:"
echo "  cd nextflow"
echo "  PATH=/path/to/conda/envs/markermag-nf/bin:\$PATH nextflow run main.nf \\"
echo "    -profile test --input test/samplesheet.csv --outdir test/results"
