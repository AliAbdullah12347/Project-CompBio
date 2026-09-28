#!/usr/bin/env bash
# Reproduce the public-data recreation of Boltz et al. 2024 end to end. Run from this directory:
#   ./download_data.sh && ./run_all.sh
# Runtime on a 4-core laptop: ~5 min for steps 00-03, ~3 h for 04 and 06 (bMIND MCMC). Steps 04 and
# 06 checkpoint every 200 genes in data/processed/bmind_chunks/, so an interrupted run resumes.
set -euo pipefail
cd "$(dirname "$0")"
source ../env.sh

echo "=== 00  parse GEO, TPM, Boltz's expression filter ==="
python src/00_prepare_data.py

echo; echo "=== 01  published tables -> data/reference/ ==="
python src/01_extract_reference.py

echo; echo "=== 02  cell-type proportions: Table 1, logistic regressions, Fig S6 ==="
Rscript src/02_celltype_proportions.R

echo; echo "=== 03  bulk limma DE: lithium, case/control, cell-adjusted; Krebs overlap ==="
Rscript src/03_bulk_de.R

echo; echo "=== 04  bMIND cell-type expression and bmind_de ==="
Rscript src/04_bmind.R

echo; echo "=== 05  bMIND evaluation: Fig 1B, Fig S2, Table S3 vs DICE, bmind_de vs Boltz ==="
Rscript src/05_bmind_evaluation.R

echo; echo "=== 06  preprint prediction analysis (single 70/30 split) ==="
Rscript src/06_preprint_prediction.R

echo; echo "=== 07  validation scorecard ==="
python src/07_validate.py

echo; echo "Done. Tables and figures in results/"
