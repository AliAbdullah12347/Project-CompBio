# How to read the 28 preservation rows: the perm split is by design, the comparison was not

## The pattern

28 of the 31 pending rows that call `modulePreservation` use two permutation budgets:
**50 per draw** for the 20 subsampled bp_lith draws, **500** for the single bp_nolith run.

This looked like a defect when I first audited it. It is not. `METHODS.md:470-479`
pre-specifies it, with the reasoning stated:

> 20 independent subsamples of bp_lith (each n=74 from 152) gives a usable distribution
> of Zsummary values **for descriptive purposes**; 50 permutations per draw is sufficient
> for an individual null. bp_nolith, which cannot be subsampled, is run once with 500
> permutations (ten times the per-draw count) **for a stable point estimate**.

bp_nolith is the whole group at n=74 — there is nothing to subsample — so one precise
estimate is the right spend. bp_lith is 152, so the subsampling variance is itself a
quantity worth measuring, and 20 draws buys that more cheaply than 20 precise ones.

## What went wrong was the reading, not the design

The design describes the lith Zsummary distribution as **descriptive**. It never claimed
parity with the nolith point estimate. But downstream work compared them anyway —
`RESULTS.md`'s headline (now withdrawn) set inp-022's bp_lith gold of 15.55 against
base-001's bp_nolith gold of 14.01 and called it "same reference, matched instrument."
Neither half was true, and the perm split is one of the three reasons.

The scale of the effect, measured on base-001's own data: gold reads **21.74** in the
50-perm lith draws and **14.01** in the 500-perm nolith run, from one reference.

## So how should these 28 rows be read

**Valid, within one arm:**
- each module against **gold from the same run and the same draw**
- **medianRank**, which is ordinal within a run and therefore immune both to the perm
  budget and to the size inflation that lets gold score 21.74 while ranking last of 72
- the spread of the lith distribution across draws — the thing the 20 draws exist to show

**Not valid:**
- any bp_lith Zsummary against any bp_nolith Zsummary, including within a single row
- any Zsummary against Langfelder's published bands

## Where the central question actually gets answered

`base-002-ceil` is the only row that scores **one reference against two groups at one
budget**: the ctrl160 network against 74 held-out controls and against bp_nolith, both at
500 perms, n=74. It is running now, and `analyses/018-matched-ceiling-comparison/`
is written and waiting for it — with the reference-identity and instrument-match checks
enforced as preconditions that abort rather than warn.

No script changes were made. Changing 28 pre-specified permutation budgets would be a
change to the experimental design, not a bug fix, and that is the user's call.
