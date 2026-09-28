#!/usr/bin/env bash
# Fetch every input this project needs (about 29 MB total). Idempotent.
set -euo pipefail
cd "$(dirname "$0")"
GEO="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124326"
CAM="https://static.cambridge.org/content/id/urn%3Acambridge.org%3Aid%3Aarticle%3AS0033291719002745/resource/name"

mkdir -p data/raw data/supplementary
[ -f data/raw/GSE124326_count_matrix.txt.gz ]  || curl -fL -o data/raw/GSE124326_count_matrix.txt.gz  "$GEO/suppl/GSE124326_count_matrix.txt.gz"
[ -f data/raw/GSE124326_series_matrix.txt.gz ] || curl -fL -o data/raw/GSE124326_series_matrix.txt.gz "$GEO/matrix/GSE124326_series_matrix.txt.gz"
[ -f "data/supplementary/File_S2_WGCNA_modules (sup001).xlsx" ]            || curl -fL -o "data/supplementary/File_S2_WGCNA_modules (sup001).xlsx"            "$CAM/S0033291719002745sup001.xlsx"
[ -f "data/supplementary/Supplementary_Methods_Tables_Figures (sup002).docx" ] || curl -fL -o "data/supplementary/Supplementary_Methods_Tables_Figures (sup002).docx" "$CAM/S0033291719002745sup002.docx"
[ -f "data/supplementary/File_S1_DEG_lists (sup003).xlsx" ]                || curl -fL -o "data/supplementary/File_S1_DEG_lists (sup003).xlsx"                "$CAM/S0033291719002745sup003.xlsx"

echo "Verifying checksums ..."
( cd data/raw && sha256sum -c CHECKSUMS.sha256 ) || echo "  (raw checksums differ - GEO may have re-deposited)"
( cd data/supplementary && sha256sum -c CHECKSUMS.sha256 ) || echo "  (supplementary checksums differ)"
du -sh data/raw data/supplementary
