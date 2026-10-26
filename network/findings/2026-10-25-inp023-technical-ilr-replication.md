# inp-023: technical + ILR residualisation replicates inp-022's rescue

**Date:** 2026-10-25 (completed 2026-10-26T12:05Z)
**Analysis:** `network/runs/inp-023` (age+sex+rin + ILR b1–b4 residualised, 42-module reference on 234 controls, 20 lith draws)
**Evidence:** 24/25 · **Interest:** 5/10 · **Δ:** −4.6

## Design

Residualise all expression data for both technical covariates (age, sex, RIN) AND cell
composition (ILR b1–b4) before building the reference and test networks. This is the most
conservative form of residualisation: any preserved structure must survive both sources of
confounding. `inp-022` tested ILR-only; `inp-023` adds the three technical covariates
used in the de_v2 primary model.

Pre-specified covariates: `~ age + sex + rin + b1 + b2 + b3 + b4`. Plate and seqPCs
omitted to avoid over-fitting (stated in `config/inp-023.json` before computation).

| arm | expression | modules | instrument |
|:--|:--|--:|:--|
| base-001 (raw) | offset-removed, raw | 71 | 50 perms/draw, 20 draws |
| inp-022 (ILR) | offset-removed, ILR residualised | 44 | 50 perms/draw, 20 draws |
| inp-023 (this) | offset-removed, technical+ILR residualised | 42 | 50 perms/draw, 20 draws |

## Reference network

42 modules (deepSplit=4, mergeCutHeight=0, power=12, bicor, signed, n=234 controls after
doubly-residualising). Two fewer than inp-022 (44), consistent with residualising for
three additional technical covariates removing a small amount of additional co-expression
structure.

SFT curve (recorded in sft_r2_curve.csv):

| arm | first power R²≥0.80 | max R² | powerEstimate |
|:--|--:|--:|--:|
| base-001 (raw) | never | 0.770 | NA |
| inp-022 (ILR) | power 9 | 0.880 | 16 |
| inp-023 (this) | power 13 | 0.836 | NA |

Adding technical covariates on top of ILR degrades the SFT criterion slightly (0.880 →
0.836, no powerEstimate returned) — consistent with technical covariates capturing some
genuine regulatory variation whose removal reduces the scale-free structure.

## Final results (20 of 20 lith draws)

**bp_nolith** (500 perms, n=74 fixed, complete):

| metric | inp-023 | inp-022 (for reference) |
|:--|--:|--:|
| gold Zsummary | 28.08 | 30.46 |
| gold medianRank | 43/43 (last) | 45/45 (last) |
| modules above gold (Zsummary) | 25 of 42 (60%) | 22 of 44 (50%) |

**bp_lith** (50 perms per draw, 20 draws, complete):

| metric | inp-023 final | inp-022 final |
|:--|--:|--:|
| gold mean Z (across 20 draws) | **16.98** (sd 10.79) | 15.55 (sd 8.13) |
| gold medianRank | **43/43 in all 20 draws** | 45/45 in all 20 draws |
| modules with mean Z > gold mean Z | **40 of 42 (95%)** | 41 of 44 (93%) |
| real modules mean Z | 35.08 (sd 13.74) | — |
| gold Z range across draws | 6.75 – 38.21 | 7.6 – 33.1 |

Gold ranks **last** in every single draw (43/43 = rank 43 out of 43 entries, 20/20 draws),
fully matching the inp-022 pattern. 40/42 real modules beat gold by mean Zsummary; all 42
beat it by medianRank in every draw.

## Interpretation

With technical covariates (age, sex, RIN) residualised on top of ILR composition removal,
the network behaves nearly identically to inp-022:

- The 42-module reference is smaller by only 2 modules vs inp-022 (44 → 42; −5%). Technical
  covariates contributed almost nothing to the composition-driven co-expression architecture.
- Gold always ranks last: all real modules beat the random gene set in every draw.
- 95% of modules beat gold's mean Zsummary — essentially the same as inp-022 (93%).
- Draw-to-draw variance in gold Z remains high (sd 10.79 vs 8.13), confirming genuine
  heterogeneity in bp_lith subsamples that persists after full residualisation.

**This confirms the composition finding is robust to technical covariate residualisation.**
The within-arm evidence is complete; the within-arm structure is real. Cross-arm Zsummary
comparisons (base-001 vs inp-022 vs inp-023) remain invalid because the reference networks
differ (71 / 44 / 42 modules, different gene sets).

## Caveat

Cross-arm Zsummary comparisons are not on a common scale (different reference networks,
different module counts, different gene membership). The evidence base for this finding
is the within-arm comparison only: gold vs. real modules, within the same network. That
comparison is valid.

## Score

E=5 R=5 B=5 D=5 C=4

- E=5 (Evidence): 20/20 draws complete for both lith and nolith; question directly answered with
  full pre-specified data.
- R=5 (Relevance): directly addresses the arm's key question — does composition + technical
  residualisation leave any co-expression structure?
- B=5 (Robustness): gold always last across all 20 draws; replicates inp-022 pattern exactly;
  two independent residualisation specs give the same answer.
- D=5 (Decisiveness): gold last in all 20 draws is an unambiguous within-arm result.
- C=4 (Completeness): within-arm evidence is complete; −1 because cross-arm comparison
  (to base-001 / inp-022) is structurally not possible (different reference networks).
