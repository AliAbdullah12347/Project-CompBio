#!/usr/bin/env bash
# Reproduce Krebs et al. 2020 end to end. Run from this directory.
#   ./run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
source ../env.sh

echo "=== 00  parse GEO, QC, gene filter ==="
python src/00_prepare_data.py

echo; echo "=== 01  extract published tables as validation ground truth ==="
python src/01_extract_ground_truth.py

echo; echo "=== 02  differential expression (BD and lithium) ==="
Rscript src/02_differential_expression.R

echo; echo "=== 03  WGCNA co-expression network ==="
Rscript src/03_wgcna.R

echo; echo "=== 04  cell-type composition ==="
Rscript src/04_celltype.R

echo; echo "Done. Tables in results/"
