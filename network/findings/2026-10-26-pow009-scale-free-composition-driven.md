# Scale-free topology fails on all groups even after offset removal; composition removal rescues it

**Date:** 2026-10-26
**Analysis:** `network/runs/pow-009` (pickSoftThreshold, all 3 groups, offset-removed expression)
**Evidence:** 22/25 · **Interest:** 7/10 · **Δ:** -1.8

## Design

`pickSoftThreshold(powerVector=1:20, networkType="signed", corFnc="bicor")` on offset-removed
expression from `expr_for()`, all 12,368 genes. Three groups: control (n=234), bp_nolith
(n=74), bp_lith (n=74 subsampled at seed 20261009). The SFT R² curves are themselves the
deliverable (per spec). Config pre-specified in `network/config/pow-009.json`.

## Results

| group | n | max R² | at power | powerEstimate | threshold (0.80) |
|:--|--:|--:|--:|--:|:--|
| control | 234 | **0.770** | 20 | NA | **FAILS** |
| bp_nolith | 74 | **0.755** | 20 | NA | **FAILS** |
| bp_lith | 74 | **0.701** | 20 | NA | **FAILS** |

All groups fail the scale-free criterion (R²<0.80) at every power in 1–20. `powerEstimate`
returns NA (function reports no valid power), which is the "honest" failure mode — compare
Phase 0 on raw expression, where `powerEstimate` returned 1 (a silently wrong answer).

Full curves in `network/runs/pow-009/sft_r2_all_groups.csv`.

## Context: raw vs offset-removed vs composition-removed

| expression | control max R² | bp_nolith max R² | bp_lith max R² | fits (≥0.80)? |
|:--|--:|--:|--:|:--|
| Phase 0 (raw, pre-offset removal) | 0.66 | 0.41 | 0.54 | No |
| pow-009 (offset-removed, this run) | 0.770 | 0.755 | 0.701 | No |
| inp-022 SFT (ILR b1–b4 residualised) | — | — | — | Yes (power 9, R²=0.802) |
| inp-023 SFT (technical + ILR) | — | — | — | Yes (power 13, R²=0.836) |

Offset removal substantially improves the R² curves (control: 0.66 → 0.770; bp_nolith:
0.41 → 0.755; bp_lith: 0.54 → 0.701) and eliminates the non-monotone behaviour seen in
raw expression. But all three groups still fall below the 0.80 threshold.

Composition removal (inp-022, inp-023) crosses the threshold. This establishes that the
residual scale-free failure on offset-removed expression is **composition-driven**: the
co-expression structure violates the WGCNA assumption because cell-type proportions
create a dense, composition-correlated backbone that is neither scale-free nor modular.

## Interpretation

Scale-free failure is not fully an artifact of the library offset. It is substantially
driven by cell composition variance in the expression matrix. Once composition is
accounted for (ILR residualisation in inp-022/inp-023), the network satisfies the
scale-free criterion. This corroborates the arm's central finding (#2, E=24/25):
73% of the lithium network dimension is cell composition.

Using power=12 (WGCNA FAQ guideline for signed networks, n>40, when no power satisfies
scale-free) is therefore appropriate: the criterion fails not because of small samples or
unusual biology, but because cell composition variance is a specific known confound that
can be removed.

## Relationship to retracted finding #27

Finding #27 (`2026-10-09-scale-free-fails.md`) was retracted because it was measured on
raw expression and the retraction note claimed "scale-free does not fail on offset-removed
expression." pow-009 shows this was an overstatement: scale-free still fails (max R²=0.770
vs 0.80 threshold), but the failure is milder and the mechanism is now identifiable (composition).

The original finding's claim (scale-free fails) stands for all currently tested expressions
except the residualised networks. The retraction should be read as "the raw-expression failure
included an artifact component" rather than "scale-free fits on corrected expression."

## Score

E=5 R=4 B=4 D=4 C=5

- E=5: Three groups, all 20 powers, offset-removed expression, full CSV output pre-specified.
- R=4: Consistent with Phase 0 (direction confirmed) and inp-022 SFT (composition rescue
  confirmed). Robustness reduced from 5 because the bp_lith subsample is a single draw.
- B=4: Composition as the driver is strongly supported by inp-022/inp-023 SFT; −1 because
  no direct test of the composition→scale-free mechanism in this run.
- D=4: All three groups fail — consistent, directional. −1 because the near-miss at 0.770
  (just below 0.80) leaves slight ambiguity about the exact threshold.
- C=5: Complete — all three groups, full grid, output CSV, config pre-specified.
