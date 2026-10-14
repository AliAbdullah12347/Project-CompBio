# Darkturquoise: control-lith contrast survives composition adjustment; nolith-lith does not

**Date:** 2026-10-14
**Analysis:** `network/analyses/007-darkturquoise-across-groups`
**Evidence:** 17/25 · **Interest:** 8/10 · **Δ:** +0.8

## Claim

Extending to all 460 samples (control n=234, bp_nolith n=74, bp_lith n=152) reveals that
the darkturquoise module eigengene (the ribosomal module identified in analyses 005–006)
differs between groups, and that the pattern of significance after composition adjustment
is asymmetric:

| Contrast | p_unadjusted | p_bonf_unadjusted | p_bonf_adjusted | verdict |
|:--|--:|--:|--:|:--|
| control vs bp_nolith | 0.549 | 1.000 | 0.831 | NS both ways |
| control vs bp_lith | < 0.0001 | < 0.001 | **0.001** | **survives adjustment** |
| bp_nolith vs bp_lith | 0.005 | 0.015 | 0.339 | attenuates |

Group means (darkturquoise ME, unadjusted):
control = +0.0083, bp_nolith = +0.0048, bp_lith = −0.0129.

## What this means

Analysis 002 (BP-only, n=226) found the lithium~PC2 association disappears after
composition adjustment. That result stands — within BP patients, the nolith-lith
darkturquoise difference is entirely compositional (p_bonf_adjusted = 0.339).

But when controls are included, the control-lith contrast survives composition adjustment
(p=0.001). Two possible explanations:

1. **Lithium has a real effect on ribosomal gene expression** that is independent of
   cell composition as measured by ILR b1–b4. The direction (lith lower than control)
   is consistent with a lithium-mediated downregulation of ribosomal transcription.

2. **Unmeasured composition**: Lithium-treated patients may differ from controls in cell
   subtypes or states not resolved by LM22's 22-cell-type palette. The ILR adjustment
   may over-correct within BP (where the composition-darkturquoise relationship is
   strong: r_b2 = −0.444) but under-correct for control-vs-lith (r_b2 in control = −0.130).

The composition-darkturquoise correlation is much weaker in controls (r_b2 = −0.130) than
in BP patients (bp_nolith r = −0.390, bp_lith r = −0.444). This asymmetry could reflect:
- A disease-state change in the ribosomal module's cell-type dependency
- Different B-cell activation states in controls vs BP

## Overall composition R²

Across all 460 samples, ILR b1–b4 explain 32.6% of darkturquoise variance (R²=0.326).
The composition-adjusted ANOVA has higher F (15.3) than unadjusted (10.6), because
composition partially confounds the group comparison in opposite directions.

## Scoring

E=4 (pre-specified tests; clean linear model; N=460; but three-group comparison with
4-ILR composition model; can't rule out unmeasured composition as driver of adjusted result).
R=4 (extends darkturquoise beyond BP-only; reveals control-lith asymmetry not visible in
analysis 002; changes the picture).
B=2 (single cohort, cross-sectional; no replicate).
D=3 (the asymmetric pattern is specific and interpretable, but the causal claim is uncertain).
C=4 (all pairwise reported; unadjusted and adjusted; data written).

Interest 8: genuinely surprising — analysis 002 found "all composition," but extending to
controls reveals a residual that survives. Δ +0.8 — ahead of the data given the unmeasured-
composition alternative explanation.

## Limits

The ILR balances were estimated from the same RNA (potential circularity, addressed in
analysis 003 for the BP-only context but not repeated here for all 474 samples).

The adjusted group comparison uses composition as a linear covariate, but the relationship
between composition and darkturquoise ME is not linear (r varies across groups), so a
linear adjustment may be inadequate. A group × composition interaction test would be the
appropriate follow-up.
