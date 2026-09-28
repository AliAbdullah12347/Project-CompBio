# Boltz et al. 2024 — reference recreation

A from-scratch reimplementation of

> Boltz T, Schwarz T, Bot M, Hou K, Caggiano C, Lapinska S, Duan C, Boks MP, Kahn RS, Zaitlen N,
> Pasaniuc B, Ophoff R (2024). **Cell-type deconvolution of bulk-blood RNA-seq reveals biological
> insights into neuropsychiatric disorders.** *American Journal of Human Genetics* 111(2), 323–337.
> [doi:10.1016/j.ajhg.2023.12.018](https://doi.org/10.1016/j.ajhg.2023.12.018)

and of the lithium-prediction analysis that appears only in its bioRxiv preprint
([10.1101/2023.05.24.542156](https://doi.org/10.1101/2023.05.24.542156), v1, 25 May 2023).

The authors released no code. **Their data is controlled-access** (dbGaP phs002856.v1: 1,730
samples plus genotypes), so this folder runs their methods on the public **GSE124326** cohort —
444 post-QC samples from the same Utrecht cohort and the same lab (Krebs et al. 2020) — and checks
each result against their published tables. Everything that needs genotypes or medical records is
out of reach and is listed, with the reason, in [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).
Read that file before relying on any number here.

---

## Layout

```
Boltz Paper/
├── README.md               ← this file
├── REPRODUCIBILITY.md      ← what is and is not reproducible, every deviation, and why
├── STATUS.md               ← scorecard and run log
├── download_data.sh        ← fetch all public inputs (~190 MB)
├── run_all.sh              ← run the whole pipeline
├── paper/
│   ├── Boltz2024_AJHG.pdf
│   └── Boltz2023_bioRxiv_preprint_v1.xml     ← preprint full text (JATS)
├── docs/                   ← extracted text of paper, supplement and preprint, for grepping
├── data/
│   ├── raw/                ← GEO downloads, never edited (+ SHA-256)
│   ├── supplementary/      ← all 17 AJHG supplementary files, renamed by content (+ SHA-256)
│   ├── external/           ← GENCODE v19 GTF; DICE (Schmiedel 2018) TPM for Table S3 (+ SHA-256)
│   ├── reference/          ← the authors' published results as CSV — the validation target
│   └── processed/          ← derived; regenerable from src/
├── src/
│   ├── 00_prepare_data.py          ← parse GEO, TPM from counts, Boltz's expression filter
│   ├── 01_extract_reference.py     ← supplementary xlsx + main-text numbers → CSV
│   ├── 02_celltype_proportions.R   ← Table 1, cell-type selection, logistic regressions, Fig S6
│   ├── 03_bulk_de.R                ← limma-trend DE (Fig 4C, Fig S8), cell-adjusted DE, Krebs overlap
│   ├── 04_bmind.R                  ← bMIND cell-type expression; bmind_de
│   ├── 05_bmind_evaluation.R       ← Fig 1B, Fig S2, Table S3 vs DICE, bmind_de vs Boltz
│   ├── 06_preprint_prediction.R    ← the preprint's 70/30-split lithium prediction
│   └── 07_validate.py              ← scorecard against the paper
└── results/                ← all output tables and figures
```

## Running it

Uses the shared R 4.6.1 library at `C:\Users\hp\R\win-library\4.6` (limma, edgeR, plus MIND 0.3.3
from github.com/randel/MIND and its dependencies MCMCglmm, nnls, matrixcalc, Biobase, BisqueRNA).
Python needs only pandas and numpy.

```bash
./download_data.sh   # idempotent, verifies checksums
./run_all.sh         # ~2.5 h on a 4-core laptop, almost all of it bMIND MCMC; resumable
```

## What reproduces, in one table

Full detail in `results/validation_scorecard.csv` (49 claims) and `REPRODUCIBILITY.md`.

| Boltz claim | Boltz (n = 1,730) | Here (n = 444) | Verdict |
|:--|:--|:--|:--|
| Neutrophils higher in lithium users | p = 1.5e-9 | p = 0.012 | replicated |
| Neutrophils higher in cases; NK higher in controls | p = 2.3e-8; 1.2e-7 | p = 2.8e-6; 3.2e-4 | replicated |
| Lithium non-users vs controls: no cell type differs | none | none of their 8 | replicated |
| Bulk lithium DEGs | 100 | 46 | consistent (power) |
| Their 100 lithium DEGs, same sign here | — | 85%, logFC r = 0.80 | replicated |
| 33 of their lithium DEGs in Krebs | OR 6.43 | 32, OR 5.92, using Krebs' *text* list | replicated |
| bmind_de lithium: no cell-type DEG | 0 | 0 | replicated |
| bmind_de case/control: 21 naive B, 4 neutrophil DEGs | 21 / 4 | 0 / 0 | differs |
| bMIND vs scRNA (DICE) R², naive CD4 | 0.84 | 0.52 (bulk baseline 0.31) | consistent, lower |
| Preprint: expression adds nothing over proportions | — | lithium AUC 0.59 vs ≤ 0.60 | replicated |

## Four things the paper does not say

1. **The prediction result is preprint-only.** The single-70/30-split lithium prediction ("gene
   expression does not provide additional predictive value over the cell type proportions") is in
   the bioRxiv v1 but was removed from the AJHG paper. Over 100 random splits here, the same
   proportions-only model's test AUC ranges 0.51–0.72: one split cannot support that claim.
2. **Their Krebs comparison used a list Krebs never deposited** — the documented-model lithium DEGs
   from Krebs' text, not either deposited sheet.
3. **The case/control logFC range in the text is wrong** (−0.126 stated; −0.217 deposited).
4. **Their bMIND validation lacks a baseline.** Raw bulk expression correlates with DICE monocytes
   almost as well as "bMIND monocytes" do (R² 0.42 vs 0.45). Deconvolution adds real cell-type
   signal for naive CD4 T cells, little for monocytes, none for naive B cells.

## Relevance to our project

- **Arm 1** exists because of finding 1. The preprint recreation here is the baseline to beat, and
  it shows each flaw concretely: gene selection and residualisation before the split, one split,
  and more genes than training samples (AUC stops changing beyond ~200 genes because `glm` drops
  the aliased columns).
- **Illness retest:** in GSE124326, case/control DE with only PCs as covariates is mostly
  assessment-group batch (236 DEGs → 80 with group adjusted). Any case/control model must adjust
  for it.
- **Deconvolution:** with the deposited LM22 fractions, memory B and CD8 are effectively absent, so
  any mediator or feature built on them is empty. Monocytes are ~4× overestimated.
- **bmind_de has an exact closed form** (`src/04_bmind.R`, `REPRODUCIBILITY.md` §3) that agrees with
  the MCMC at r ≥ 0.989 and runs in seconds instead of hours — usable inside CV folds.
