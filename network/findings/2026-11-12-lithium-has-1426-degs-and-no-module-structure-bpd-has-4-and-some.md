# Lithium has 1,426 DEGs and no module structure; bipolar has 4 and some

**Date:** 2026-11-12
**Analyses:** `network/experiments/perm-054.R` (BPD) and `perm-059.R` (lithium), K=500 each
**Evidence:** 22/25 · **Interest:** 9/10 · **Δ:** +1.0

## The pair

The same test, the same partition, the same 500-draw sample-label permutation
null, run on both contrasts:

| contrast | DEGs at BH<0.05 | modules enriched | null mean | draws >= obs | **p** |
|:--|--:|--:|--:|--:|--:|
| **WB_LI** (lithium) | **1,426** | 21 of 69 | 16.4 | 96 of 500 | **0.194** |
| **WB_BPD** (bipolar) | **4** | 27 of 69 | 16.5 | 8 of 500 | **0.018** |

**The contrast with 1,426 differentially expressed genes shows no module-level
concentration beyond chance. The contrast with four does.**

This is the reverse of what I expected before running it, and it is the more
interesting way round.

## Why it is not a contradiction

A module rank-enrichment test is *relative*: it asks whether a module's genes
sit higher in the DE ranking than every other gene. A contrast that moves
almost everything a little raises all modules together, so none stands out —
large, diffuse and unlocalised looks exactly like nothing here.

That is precisely what this arm has concluded lithium does. `analyses/002`
put 73% of the lithium-associated co-expression dimension down to cell
composition; `trait-053` found 15 of 71 marginal module associations
collapsing to 1 under the composition balances. A shift in the cell mixture
moves the whole transcriptome at once. **Lithium's 1,426 DEGs are real and
they are global, which is why they do not concentrate in modules.**

Bipolar diagnosis is the opposite shape: almost no individual gene clears
correction, but what little moves appears to be structured — concentrated in
particular co-expression modules rather than spread evenly.

## What each number is worth

**Lithium: withdrawn.** `019`'s "21 of 69 modules enriched for lithium" is
within what co-expression alone produces. It should not be reported.

**Bipolar: kept, with a smaller headline.** 27 against an expected 16.5,
p = 0.018, about eleven modules above chance. The old null — random
size-matched gene sets, median 0 — was wrong by sixteen modules because
uncorrelated gene sets cannot see the inflation that correlated ones cause.
"27 of 69 modules are BPD-enriched" implies a comparison against zero and
must not be written.

## Why this matters to the project's argument

The project's thesis is that a bipolar blood signature is really a lithium
signature. This result sharpens that rather than contradicting it, by
separating two things that were being conflated:

- **In magnitude**, lithium dominates completely: 1,426 genes against 4.
- **In structure**, lithium contributes nothing at module level and bipolar
  contributes a little.

So the right statement is not "bipolar has no signal" but **"bipolar's signal
is too small to see gene by gene and is only visible as coordinated structure,
while lithium's is large, diffuse, and compositional."** Those are different
kinds of effect, and a study that looks only at DEG counts sees the second
and misses the first entirely.

## Caveats

- **p = 0.018 is one marginal test.** It is a sentence in a discussion, not a
  section, and this dataset has been examined hard.
- **The two rows are not formally compared.** 27 against 21 on the same null
  is suggestive, but they come from different sample sets (n=308 and n=226)
  and no test of the difference has been run.
- **Both permutations are mildly optimistic.** Covariates stay with their
  subjects while the label moves, so each permuted design also loses the real
  label-covariate associations and is slightly better conditioned than the
  truth. Both p-values are therefore generous, which matters more for the one
  that passed.
- The BPD contrast runs against the QC-filtered control group, which is not
  random with respect to diagnosis (Fisher p = 3.04e-10).
- These are `grp-051`'s 69 modules. Per `ann-058`, only 6 of 69 reconstruct
  across decompositions above Jaccard 0.5, so module-level results here are
  partly decomposition-specific.

See also [[2026-11-12-bpd-module-enrichment-survives-but-the-old-null-was-wrong-by-16-modules]]
and [[2026-10-14-lithium-dimension-is-composition]].
