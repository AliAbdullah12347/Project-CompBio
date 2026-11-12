# The BPD module enrichment survives an honest null — but the null it was originally tested against was wrong by sixteen modules

**Date:** 2026-11-12
**Analysis:** `network/experiments/perm-054.R` -> `network/runs/perm-054` (K=500)
**Evidence:** 22/25 · **Interest:** 8/10 · **Δ:** +0.6

## The result

**27 of 69** modules are rank-enriched for BPD-associated genes. Against a
sample-label permutation null that recomputes the entire DE fit 500 times:

| | |
|:--|--:|
| observed | **27** |
| null mean | **16.5** (sd 4.9) |
| null median | 17 |
| null 95th percentile | 24 |
| null maximum | 30 |
| draws reaching 27 | 8 of 500 |
| **empirical p (add-one)** | **0.018** |
| p floor at K=500 | 0.002 |

The p is nine times the floor, so unlike `perm-031`, `perm-032`, `hub-037`
and `cent-056` this is a measurement with real resolution rather than a
budget artifact. **The enrichment is real.**

## And the null it used to be tested against was badly wrong

`analyses/019`'s original null drew random size-matched gene sets and gave
**median 0, 95th percentile 1**. The honest null's median is **17**.

The reason is the one the test was designed around: a module's genes are
co-expressed by construction, so their DE statistics move together, and a
Wilcoxon that assumes independence is anti-conservative on them. Random gene
sets are uncorrelated, so that null cannot see the inflation at all. It ruled
out module size, which was never the threat.

Permuting the diagnosis label leaves the expression matrix untouched -- every
gene-gene correlation survives intact -- and destroys only the group signal.
Sixteen of the twenty-seven modules are what that correlation structure
produces from noise.

## What the claim becomes

Before: *27 of 69 modules are BPD-enriched, against a null of essentially
zero* -- a major distributed signal.

After: *27 against an expected 16.5, p = 0.018* -- about **eleven modules
above chance**, nominally significant and nothing more.

**The paradox that made this interesting largely dissolves.** The puzzle was
that whole blood yields only **4** BPD DEGs at BH<0.05 against **1,426** for
lithium, yet BPD lit up *more* modules (27) than lithium (21). Most of that
gap is now accounted for: a contrast with almost no real signal still lands
~16 modules by correlation alone, so 27 is a modest excess rather than
evidence of a large hidden signal.

## How to state it

> Module-level enrichment for bipolar diagnosis exceeds a sample-label
> permutation null (27 observed against 16.5 expected, p = 0.018, K = 500),
> but the margin is small, and the previously reported contrast with lithium
> reflects an earlier null that assumed gene independence.

The bare sentence "27 of 69 modules are enriched for BPD" should not appear
anywhere. On its own it implies a comparison against zero.

## What is still missing

`perm-059` is running the identical test on the **lithium** contrast, which
`019` scored against the same invalid null. Until it lands the BPD-versus-
lithium comparison is not symmetric -- the number that made BPD look
surprising has been corrected while the number it was surprising *relative
to* has not. With 1,426 DEGs lithium should clear its own null comfortably,
but that is an expectation, not a result.

## Caveats

- **The permutation is mildly optimistic.** Covariates stay with their
  subjects while the label moves, so each permuted design also loses the real
  diagnosis-tobacco association (19% of controls, 32.7% of BP1) and is
  slightly better conditioned than the truth. That makes the null easier to
  beat, so p = 0.018 is if anything generous.
- **p = 0.018 is one marginal result** in a dataset that has been examined
  hard. It is a sentence in a discussion, not a section.
- The contrast runs against the QC-filtered control group, which is not
  random with respect to diagnosis (Fisher p = 3.04e-10).

See also [[2026-11-09-perm031-at-157-draws-is-a-clean-negative]].
