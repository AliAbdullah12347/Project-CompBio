# Darkturquoise group×composition interaction: the control-lith contrast survives proper slope correction

> **CORRECTION 2026-10-15 — `b2` was mislabelled.** The sequential binary partition in
> `mediation/scripts/00b_ilr_trace.R:96-99` defines **b2 = gran vs mono**, not a B-cell or
> lymphocyte balance (b1 = myeloid vs lymphoid, b3 = T vs NK+B, b4 = NK vs B). Labels below
> have been corrected in place. Two directional claims built on the wrong label were also
> wrong and are flagged inline: b2 by group is bp_lith 0.3079 > control 0.2122 > bp_nolith
> 0.1816 (lithium is **highest**, not lowest), and r(darkturquoise, b2) is **negative** in
> every group (-0.130 ctrl, -0.390 nolith, -0.444 lith). The corrected reading — lithium
> raises granulocytes, granulocytes dilute the lymphocyte ribosomal program, darkturquoise
> falls — fits the arm's own r(darkturquoise, gran) = -0.510 (p = 2.3e-16).

**Date:** 2026-10-14
**Analysis:** `network/analyses/008-darkturquoise-interaction`
**Evidence:** 18/25 · **Interest:** 7/10 · **Δ:** −0.2

## Claim

Analysis 007 used an additive composition adjustment, but the b2 (granulocyte/monocyte)
slope on darkturquoise ME differs significantly between groups. This analysis tests
whether the surviving control-lith contrast (analysis 007) is an artefact of incorrect
slope assumptions.

## Per-group b2 slopes

| Group | n | slope_b2 | 95% CI | r_b2 |
|:--|--:|--:|:--|--:|
| control | 234 | −0.027 | [−0.053, −0.001] | −0.130 |
| bp_nolith | 74 | −0.063 | [−0.097, −0.028] | −0.390 |
| bp_lith | 152 | −0.071 | [−0.094, −0.048] | −0.444 |

The b2 slope in bp_lith is ~2.6× steeper than in controls. The group×b2 interaction
is significant (LRT F(2,451)=3.53, p=0.030); the individual interaction term for
`groupbp_lith:b2` = −0.041 (p=0.0083).

## Pairwise contrasts in interaction model

Computed via the variance-covariance matrix of the interaction model
`dt_ME ~ group * b2 + b1 + b3 + b4`; Bonferroni correction n=3.

| Contrast | Δ (marginal mean) | SE | t | p_bonf |
|:--|--:|--:|--:|--:|
| control vs bp_nolith | −0.006 | 0.0051 | −1.22 | 0.670 (NS) |
| **control vs bp_lith** | **−0.014** | **0.0041** | **−3.38** | **0.002** |
| bp_nolith vs bp_lith | −0.008 | 0.0055 | −1.38 | 0.509 (NS) |

Contrasts evaluated at mean b1, b2, b3, b4 across all 460 samples.

## Conclusion

The control-lith difference in the darkturquoise module survives proper, group-specific
composition adjustment (p_bonf=0.002). It is NOT an artefact of assuming equal b2 slopes.
The nolith-lith contrast remains NS under all adjustment models.

The remaining open question is whether this control-lith difference reflects:
(a) a genuine lithium effect on ribosomal gene co-expression independent of cell composition, or
(b) compositional differences not captured by LM22's 22-cell palette.

## Scoring

E=4 (pre-specified interaction test; vcov-based contrast SE; N=460; but cross-sectional,
single cohort, same RNA composition estimate).
R=4 (confirms analysis 007 was not an artefact; refines the mechanistic story).
B=2 (single cohort, cross-sectional; no replication).
D=4 (interaction significant; control-lith contrast robust across model specifications).
C=4 (all contrasts reported; proper SE computation documented).

Interest 7: the robustness check resolves a legitimate methodological concern from 007.
Δ −0.2 — evidence and interest well-aligned.
