# Removing the library offset rebuilds the network into a well-behaved one

**Date:** 2026-10-18
**Analysis:** `network/analyses/016-normalisation-decision`
**Evidence:** 22/25 · **Interest:** 9/10 · **Δ:** new

## Claim

The decision left open on 2026-10-15 — whether the 49 queued experiments are worth
running on the un-normalised matrix — resolves clearly: **no.** Rebuilding the reference
network with the per-sample `log2(TPM+1)` offset removed changes every structural
property that the arm has been reporting as a finding.

Control samples (n=234), top 4,000 most-variable genes, identical WGCNA parameters
(power 12, signed, bicor, deepSplit 4, minModuleSize 30, mergeCutHeight 0):

| | RAW (what the arm used) | offset removed |
|:--|--:|--:|
| mean gene–gene r | +0.425 | **+0.042** |
| % random pairs r > 0.5 | 41.9% | **8.8%** |
| best scale-free R² | 0.670 (power 20) | **0.769 (power 18)** |
| fit slope at power 12 | −0.603 | **−1.19** |
| modules detected | 25 | **35** |
| % genes unassigned (grey) | 0.9% | 5.3% |
| median \|r\| between eigengenes | 0.616 | **0.286** |
| effective dimensions | **2.33** of 25 | **5.06** of 35 |

## Why each row matters

**The scale-free curve changes shape, not just height.** The raw data gives a
non-monotone R² curve — 0.751 at power 4, collapsing to 0.007 at power 8, climbing back
to 0.670 at power 20. That is not a soft-threshold curve; it is what a degenerate
correlation matrix produces. Offset-removed gives a clean monotone rise, 0.214 → 0.769.

**The slope is the diagnostic.** Scale-free topology wants a fit slope near −1 to −2.
At the arm's chosen power 12 the raw network has slope **−0.603**; offset-removed has
**−1.19**. Finding #16 ("scale-free topology fails at every power") was reporting the
artifact, not a property of blood co-expression.

**Effective dimensionality more than doubles** (2.33 → 5.06) and median inter-eigengene
|r| halves (0.616 → 0.286). Analysis 002's headline observation that 33 modules span only
~6 effective dimensions is substantially a statement about the offset.

**Only 0.9% of genes were unassigned in the raw network.** In a network where 42% of
random gene pairs exceed r = 0.5, essentially everything joins something. The
offset-removed network leaves 5.3% grey, which is ordinary.

## Consequence

The 49 queued preservation experiments all read `expr_for()`, which returns raw
`log2(TPM+1)`. Running them produces preservation statistics on a network whose module
structure is substantially an artifact — and the offset **differs by group**
(control 3.4637, bp_nolith 3.2688, bp_lith 3.1811; Welch p = 1.14e-09), so between-group
preservation is contaminated at the source, not merely noisy.

**Recommendation: fix `expr_for()` before spending the queue's ~100 h of compute.**
That means removing the per-sample offset (or applying the TMM factors already loaded
and unused in `d$tmm`), invalidating `network/cache/`, and rebuilding the reference.
It is a rebuild, not a patch, and it is the user's call — but the evidence that it
changes the answer is now direct rather than inferred.

## Caveat

Bounded by design: 4,000 most-variable genes rather than all 12,368, run at 1 thread so
it would not starve the concurrent `base-001` job on this 4-core machine. The gene-gene
density figures therefore differ from the all-gene measurement in
`2026-10-15-library-offset-artifact.md` (+0.535 / 61.9%) because that used all genes;
the direction and magnitude of the change are the same. A full-gene rerun should
accompany any actual rebuild.

Note also that an initial version of `run.R` computed `cor()` on a genes × samples
matrix without transposing, which measures sample–sample rather than gene–gene
correlation and returned 100% for both arms. Caught, fixed, and the CSV corrected
before these numbers were recorded.

E=5 R=5 B=4 D=5 C=3
