#!/usr/bin/env bash
# ==========================
# Template: single-gene strip plot
# ==========================
# Copy this file, edit the variables below, then run:
#   bash plot_single_gene.run.sh
#
# RESULTS_DIR is the --outdir of a finished da.proteome.R run.
# By default every <comparison>_limma.csv in RESULTS_DIR/data/de_data is used;
# set LIMMA_FILES to a comma-separated list to restrict the comparisons shown.

set -euo pipefail

# ---- Edit these ----
PIPELINE_DIR="/path/to/proteomics_core_pipelines"
RESULTS_DIR="/path/to/results"
SAMPLESHEET="/path/to/samplesheet.csv"
GENES=("CD59")                # one or more gene symbols or UniProt accessions
OUT_DIR="${RESULTS_DIR}/figures/single_gene"
LIMMA_FILES=""                # leave empty to use all comparisons
WIDTH=7                       # inches
# --------------------

IMPUTED="${RESULTS_DIR}/data/protein_imputed_matrix.csv"
RAW="${RESULTS_DIR}/data/protein_raw_matrix.csv"

# Mark imputed values when the pre-imputation matrix exists (older runs lack it)
RAW_ARGS=()
if [[ -f "${RAW}" ]]; then
  RAW_ARGS=(--raw "${RAW}")
fi

if [[ -z "${LIMMA_FILES}" ]]; then
  LIMMA_FILES=$(ls "${RESULTS_DIR}"/data/de_data/*_limma.csv | paste -sd, -)
fi

mkdir -p "${OUT_DIR}"

for GENE in "${GENES[@]}"; do
  Rscript --vanilla "${PIPELINE_DIR}/inst/scripts/plot_single_gene.R" \
    --gene "${GENE}" \
    --imputed "${IMPUTED}" \
    --limma "${LIMMA_FILES}" \
    --samplesheet "${SAMPLESHEET}" \
    --out "${OUT_DIR}/${GENE}_strip_plot.png" \
    --width "${WIDTH}" \
    "${RAW_ARGS[@]}"
done
