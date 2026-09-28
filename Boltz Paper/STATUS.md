# Status — complete, 2026-09-28

All scripts written and run clean. A fresh end-to-end run from an empty `data/processed/` and
`results/` reproduced all 24 result files byte-for-byte (bMIND's MCMC is seeded per gene).

## Pipeline

| Script | Output | Runtime (4-core laptop) |
|:--|:--|--:|
| `00_prepare_data.py` | `data/processed/` — 13/13 fact checks | 1 min |
| `01_extract_reference.py` | `data/reference/` — Tables S1–S3, S13–S15 + 51 main-text numbers | <1 min |
| `02_celltype_proportions.R` | Table 1, 3 logistic-regression contrasts, Fig S6 | <1 min |
| `03_bulk_de.R` | bulk DE, cell-adjusted DE, batch sensitivity, Krebs overlap, Fig 4C, Fig S8 | 1 min |
| `04_bmind.R` | exact bmind_de (all genes) + real bmind_de check (400 genes) + real bMIND (2,000 genes) | 85 min |
| `05_bmind_evaluation.R` | Fig 1B, Fig S2, Table S3 vs DICE, bmind_de vs Boltz | <1 min |
| `06_preprint_prediction.R` | preprint 70/30 prediction + 100-split check | 35 min |
| `07_validate.py` | `results/validation_scorecard.csv` | <1 min |

Steps 04 and 06 checkpoint every 200 genes, keyed by a hash of inputs, code and MIND version.

## Scorecard (49 claims)

| Verdict | n | Meaning |
|:--|--:|:--|
| EXACT | 10 | paper text vs its own deposit; real vs exact bmind_de |
| REPLICATED | 14 | same direction and significant here |
| CONSISTENT | 11 | same direction or conclusion, weaker (n 444 vs 1,730) |
| DIFFERS | 12 | each with a stated cause, mostly the deposited LM22 fractions |
| DISCREPANCY | 2 | the paper disagrees with its own deposited tables |

Headline numbers are in `README.md`; every deviation and skip is in `REPRODUCIBILITY.md`.

## Decisions made by the user (2026-09-28)

No dbGaP access → GSE124326 substitute. Deposited CIBERSORT fractions instead of CIBERSORTx.
Installed MIND 0.3.3 and dependencies (~14 MB). Downloaded GENCODE v19 GTF and DICE TPM (3 types).

## Not done

- Everything needing genotypes or medical records (see `REPRODUCIBILITY.md` §1).
- Table S3 BLUEPRINT rows (large download, not requested).
- The preprint's Supplementary Table 2 values: an image behind bioRxiv's bot wall; the Chrome
  extension was not connected. Comparison with it is qualitative only.
