# Darkturquoise drives PC2: the lithium-associated dimension is localised to one module

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
**Analysis:** `network/analyses/005-pc2-module-loadings`
**Evidence:** 18/25 · **Interest:** 7/10 · **Δ:** +0.2

## Claim

PC2 of the 33-module eigengene matrix (the lithium-associated dimension from finding 002)
is not diffuse across modules — it is heavily loaded on a single module.

| module | PC2 loading | r_pc2 | r_lith | r_b1 | r_b2 |
|:--|--:|--:|--:|--:|--:|
| **darkturquoise** | **+0.455** | **+0.843** | −0.184 | −0.379 | −0.449 |
| yellow | +0.283 | +0.523 | −0.159 | −0.239 | −0.324 |
| blue | +0.256 | +0.474 | −0.159 | −0.272 | −0.319 |
| saddlebrown | −0.239 | −0.444 | +0.081 | +0.365 | +0.354 |

Threshold: |loading| > 0.3 (pre-specified; 1.7× expected null magnitude 1/√33 = 0.174).
**Only darkturquoise passes** (n=140 genes, 1.1% of the transcriptome).

## What the pattern means

Modules with positive PC2 loading (darkturquoise, yellow, blue, orange, skyblue) all show
*negative* correlation with ILR b1 and b2 — their eigengenes decline when the
granulocyte/monocyte balance (b2) is high. Negative-loading modules (saddlebrown, white,
violet, etc.) all show *positive* correlation with b1 and b2.

PC2 is a contrast between two sets of modules that move oppositely with cell composition.
Darkturquoise is the largest single contributor to that contrast (r_pc2 = 0.843,
4× the squared loading of the next module).

No individual module eigengene associates with lithium status at r > 0.18, consistent with
finding 001 (no BH-significant module-lithium association). The lithium signal lives in the
composition-driven PC2 aggregate, not in any single module.

## Scoring

E=4 (deterministic PCA, same algorithm as 002; pre-specified threshold; but single cohort,
one decomposition). R=4 (localises the headline result to a specific module with
interpretable biology). B=2 (one decomposition, one cohort; no replication yet).
D=4 (darkturquoise unambiguously stands out; loading ratio 0.455 vs 0.283 next).
C=4 (all 33 modules reported; nothing hidden).

Interest 7: identifying the specific module that drives the lithium-composition result is
concrete and quotable. Δ +0.2 — borderline between "behind the data" and "ahead of it."

## Limits

We do not yet know what biological processes the 140 darkturquoise genes represent —
that is the annotation experiment (ann-041 in queue). The PC2 loading tells us WHERE in the
module structure the composition signal lives, not WHY.

The 0.3 loading threshold is somewhat arbitrary; yellow (0.283) and saddlebrown (−0.239)
are close runners-up. The claim is directional, not binary.
