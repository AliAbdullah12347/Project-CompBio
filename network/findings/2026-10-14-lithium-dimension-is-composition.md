# A lithium-associated co-expression dimension exists — and it is 73% cell composition

**Date:** 2026-10-14
**Analysis:** `network/analyses/002-effective-dimensionality`
**Evidence:** 22/25 · **Interest:** 9/10 · **Δ:** +0.2

## Claim
The 33 module eigengenes span ~**6** effective dimensions (Li & Ji); two PCs carry 80% of
the variance. Testing lithium against those 6 dimensions, **PC2 is associated with lithium**
(β = −0.908, p = 4.8e-4, BH = 0.0029; permutation null over 2,000 label shuffles
**p = 0.0015**).

**That dimension is predominantly cell composition.** ILR b1–b4 explain **73.3%** of PC2's
variance; lithium alone explains 5.3%. Adjusting for composition attenuates the lithium
association by **68.1%** (73.7% with technical covariates too), leaving p = 0.081.

Technical covariates account for 9%; age for 5.7%. Composition does essentially all of it.

## Why it matters
This is the arm's question answered on the network side, and it supports the project's
founding hypothesis: a blood signature attributed to bipolar disorder may be a lithium
signature acting through cell proportion rather than through any cell changing behaviour.

## Scoring
E=5 (permutation null, nested models, dimension count fixed a priori by Li & Ji, PCA
deterministic and computed without reference to lithium so there is no selection
circularity). R=5 (the arm's question, directly). B=3 (**raised from 2**: survives the LM22 circularity check in analysis 003 with the
association and composition share both intact; still one decomposition). D=5 (clear answer with a quantified
mechanism). C=4 (config, CSVs, README, permutation null; no cross-decomposition replicate).

Interest 9: the central hypothesis of the project, confirmed crisply with the mechanism
quantified. Δ +0.2 — the evidence has caught up with the story after the circularity check.

## Caveats that belong with any quotation of this
1. **Attenuated, not abolished.** p = 0.081 after full adjustment; a residual effect
   cannot be excluded at n = 226.
2. **Mediator vs confounder is not resolved here.** Lithium raises neutrophil proportion,
   so adjusting for composition removes part of lithium's own causal path. This analysis
   cannot distinguish a compositional mediator from a compositional confounder — that is
   Arm 2's question, and this should be read against `mediation/results/`.
3. One decomposition at one parameter set.
