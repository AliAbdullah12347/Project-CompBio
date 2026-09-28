#!/usr/bin/env bash
# Fetch every public input this recreation needs (about 125 MB). Idempotent: existing files are kept.
#
#   data/raw/            GSE124326 counts + series matrix (GEO)                          17 MB
#   data/supplementary/  Boltz et al. 2024 Document S1-S2 and Tables S1-S15 (Elsevier)   82 MB
#   data/external/       GENCODE v19 GTF (gene lengths for TPM)                          37 MB
#                        DICE / Schmiedel 2018 per-sample TPM, 3 cell types (Table S3)   52 MB, gzipped
#
# Boltz's own cohort (dbGaP phs002856.v1) is controlled-access and is NOT downloaded.
set -euo pipefail
cd "$(dirname "$0")"
GEO="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124326"
ELS="https://ars.els-cdn.com/content/image/1-s2.0-S0002929723004500"   # PMC serves a JS challenge; this CDN does not
GENC="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19"
DICE="https://dice-database.org/download"

get() { [ -s "$2" ] || { echo "  fetching $(basename "$2")"; curl -fsSL -A Mozilla/5.0 -o "$2" "$1"; }; }

mkdir -p data/raw data/supplementary data/external/gencode_v19 data/external/DICE_Schmiedel2018

get "$GEO/suppl/GSE124326_count_matrix.txt.gz"  data/raw/GSE124326_count_matrix.txt.gz
get "$GEO/matrix/GSE124326_series_matrix.txt.gz" data/raw/GSE124326_series_matrix.txt.gz

# Elsevier file ids (mmcN) -> descriptive names. Order and captions from the article's
# Supplemental information section.
while IFS='|' read -r id name; do
  get "$ELS-$id" "data/supplementary/$name"
done <<'EOF'
mmc1.pdf|Document_S1_Figures_S1-S8 (mmc1).pdf
mmc2.xlsx|Table_S1_DeconeQTL_eGenes (mmc2).xlsx
mmc3.xlsx|Table_S2_replication_bulk_vs_celltype (mmc3).xlsx
mmc4.xlsx|Table_S3_scRNAseq_expression_correlation (mmc4).xlsx
mmc5.xlsx|Table_S4_TWAS_coloc_Bulk (mmc5).xlsx
mmc6.xlsx|Table_S5_TWAS_coloc_MemoryB (mmc6).xlsx
mmc7.xlsx|Table_S6_TWAS_coloc_NaiveB (mmc7).xlsx
mmc8.xlsx|Table_S7_TWAS_coloc_Monocytes (mmc8).xlsx
mmc9.xlsx|Table_S8_TWAS_coloc_Neutrophils (mmc9).xlsx
mmc10.xlsx|Table_S9_TWAS_coloc_RestingNK (mmc10).xlsx
mmc11.xlsx|Table_S10_TWAS_coloc_MemoryCD4T (mmc11).xlsx
mmc12.xlsx|Table_S11_TWAS_coloc_NaiveCD4T (mmc12).xlsx
mmc13.xlsx|Table_S12_TWAS_coloc_CD8T (mmc13).xlsx
mmc14.xlsx|Table_S13_lithium_interaction_eQTL (mmc14).xlsx
mmc15.xlsx|Table_S14_Li-eGene_summary (mmc15).xlsx
mmc16.xlsx|Table_S15_limma_DE_results (mmc16).xlsx
mmc17.pdf|Document_S2_article_plus_supplement (mmc17).pdf
EOF

get "$GENC/gencode.v19.annotation.gtf.gz" data/external/gencode_v19/gencode.v19.annotation.gtf.gz
for f in B_CELL_NAIVE_TPM MONOCYTES_TPM CD4_NAIVE_TPM; do
  out="data/external/DICE_Schmiedel2018/$f.csv.gz"
  [ -s "$out" ] || { echo "  fetching $f.csv"; curl -fsSL -o "${out%.gz}" "$DICE/$f.csv"; gzip -9n "${out%.gz}"; }
done

echo "Verifying checksums ..."
( cd data/raw && sha256sum -c --quiet CHECKSUMS.sha256 ) && echo "  raw OK" || echo "  raw differ (GEO re-deposit?)"
( cd data/supplementary && sha256sum -c --quiet CHECKSUMS.sha256 ) && echo "  supplementary OK" || echo "  supplementary differ"
# gzip -n omits the file name and timestamp, so a re-download reproduces the same checksum.
( cd data/external && sha256sum -c --quiet CHECKSUMS.sha256 ) && echo "  external OK" || echo "  external differ"
du -sh data/raw data/supplementary data/external
