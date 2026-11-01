# ann-041: module overlap with lithium / bipolar DEG lists

**Date:** 2026-11-01
**Analysis:** `network/runs/ann-041` — Fisher exact test of 71 base-001 modules against three DEG lists (de_final results), BH across modules within each list
**Evidence:** 16/25 (E=3 R=4 B=2 D=3 C=4) · **Interest:** 5/10 · **Δ:** −1.4

## Claim
Are the modules enriched for genes called differentially expressed for lithium or bipolar disorder?

## Results (background N=12,368)
| DEG list | genes in background | modules BH<0.05 |
|:--|--:|:--|
| WB_BPD | 4 | none (too few genes to test; uninformative) |
| CT_LI | 120 | yellow (20 overlap vs 5.3 expected, OR 4.5), darkred, lightyellow, mediumpurple2 (4 vs 0.34, OR 13.6) |
| WB_LI | 1426 | 15 modules; strongest tan (125 vs 36.3, OR 5.4), darkred (84 vs 20.3, OR 7.4), lightyellow, white, yellow |

Total 19 module-list overlaps significant (4 + 15). The CT_LI hits (yellow, darkred, lightyellow, mediumpurple2) are all also WB_LI hits.
Full table: `results_deg_overlap.csv`.

## Interpretation and limits
- The lithium DEG list is concentrated in a subset of modules, as expected given that lithium effects are
  largely cell-composition shifts (finding #2, trait-033/035). I did not check whether these particular modules are the
  ILR-associated ones, so this is a consistent-with, not a test.
- BH is applied across modules within each list. The 19 hits are not corrected across the three lists, so this is a description, not a pooled significance claim (README §6.6/6.8).
- The CT_LI and WB_LI lists are related and not independent evidence.
- WB_BPD has 4 genes in the background: no statement about a bipolar signature is possible here.
- Gene lists were taken as supplied; I did not re-derive them.

<!-- SCORES: E=3 R=4 B=2 D=3 C=4 -->
