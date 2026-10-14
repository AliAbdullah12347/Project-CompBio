# The lithium dimension replicates across decompositions, and never survives composition

**Date:** 2026-10-14
**Analysis:** `network/analyses/004-decomposition-robustness`
**Evidence:** 22/25 · **Interest:** 7/10 · **Δ:** -1.8

## Claim
All ten cells of the deepSplit x mergeCutHeight grid were rebuilt and asked the same
question under the same a-priori rule.

- **6 of 10** decompositions show a lithium-associated dimension at BH < 0.05.
- **4 of 4 usable** decompositions (>=10 modules, largest <25%) are among them.
- It is **PC2** in every case but one.
- **0 of 6 survive composition adjustment** (BH across the six: all 0.166).

## The unexpected part
**correlation(n_modules, composition share) = +0.879.** The better the decomposition, the
*more* of the lithium dimension is composition — 73.3% at 33 modules, 27.1% at 3.

If the compositional explanation were a parameter artifact it would peak somewhere
arbitrary. Instead it strengthens monotonically as the decomposition improves, which is
what a real effect being resolved more sharply looks like. It also means the degenerate
structure `base-001` used would have *understated* the composition story, not invented it.

## Scoring
E=5 (ten decompositions, one a-priori rule, deterministic, BH applied consistently).
R=4 (tests the headline finding's generality rather than extending it). B=4 (this IS the
robustness analysis; still one cohort and one composition estimator). D=5 (replicates in
every usable cell; never survives adjustment in any). C=4.

Interest 7: the monotone relationship between decomposition quality and composition share
was not anticipated. Δ −1.8 — well supported, deliberately undramatic.

## Limits
Parameter robustness, not external validity: one dataset, one cohort, one estimator. Does
not separate a compositional mediator from a confounder (Arm 2's question).
