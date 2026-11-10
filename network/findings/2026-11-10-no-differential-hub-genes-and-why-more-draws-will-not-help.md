# No gene changes its intramodular connectivity under lithium, and raising K further would not settle it

**Date:** 2026-11-10
**Analysis:** `network/experiments/hub-037.R` -> `network/runs/hub-037` (K=500)
**Evidence:** 19/25 · **Interest:** 6/10 · **Δ:** +0.5

## The result

Per-gene change in intramodular connectivity, kME in the lithium network
against kME in the control network, tested against a null that permutes the
group labels. At **K=500 permuted draws**:

- **0 of 12,360 genes** at BH 5%. Smallest BH is **0.149** (`pink`).
- **11 genes** beat all 500 null draws, against **24.7 expected by chance** at
  a floor of 1/501. Fewer than chance, not more.
- 797 genes at unadjusted p<0.05 against 618 expected - a **1.29x** excess.
- Of the 20 smallest p-values, **17 are less connected in lithium**, 3 more.

The modest nominal excess with a one-sided direction is the global density
difference showing up at gene level: the lithium graph is 11-22% sparser
across three independent subsamples
([[2026-11-10-lithium-graph-is-sparser-in-three-independent-draws]]), so kME
drifts down almost everywhere by a little. It is not a set of hub genes.

## Why K=500 is the right place to stop

The usual move when a permutation test returns nothing is to raise K, and
this row had already been cut from its specified 500 to 20 for budget. The
power rule used elsewhere in this arm - `K > 20n/m - 1` for `m` tests tied at
the floor among `n` - says `pink` would need **K = 1488** to clear BH 5%,
about two hours of cluster time. That is affordable, and it should still not
be run.

**`m` is not a constant.** A gene sits at the floor because it beat every
draw so far, so raising K removes genes from the floor unless their effect is
real. The requirement `20n/m` therefore *recedes* as K rises for a null
result and *falls* for a real one. This row shows which case it is in: at
K=20 there were 1028 genes at the floor, which looked like it needed only
K=240; at K=500 there are 11, below the 24.7 that chance alone predicts. The
floor emptied out the way a null does.

So the honest statement is not "underpowered" but "tested at a resolution of
0.002 and found below-chance occupancy of the tail". My own earlier note that
this needed K>=240 was extrapolating from the K=20 floor count and was wrong
for exactly this reason.

## The number this row reported before today was not a measurement

At K=20 it said "0 of 12360 genes". That was true by accident and for two
reasons at once, neither of them the data:

1. **The results column was empty.** `pull_own_kme()` builds its vector with
   `sapply()` over `seq_along()`, which returns an **unnamed** vector, and
   every caller subsets it by gene id - `kme_own_ctrl[g_in_mod]`. Subsetting
   an unnamed vector by character gives NA for every element, so
   `diff_connectivity.csv` was 12,360 rows of NA and the headline came from
   `sum(NA < 0.05, na.rm = TRUE)`, which is 0 whatever the data says.
   `delta_kme.csv` was always correct, because it indexes the kME matrices,
   which carry rownames.
2. **Underneath that, the estimator was invalid.** With names attached the
   same 20 draws give 1028 genes at BH 5%, every one of them at empirical
   p = 0 from `mean(null >= obs)`. A p of 0 is not a probability. Under
   add-one the count is 0 again - the same defect found in `perm-031` and
   `perm-032` the previous day, in a third script.

Both are fixed. The answer did not change; it became an answer.

## Usable output

`delta_kme.csv` - 12,360 genes with kME in control, bp_nolith and bp_lith,
and both deltas. This is the table for anything descriptive, including
picking candidate hubs. `diff_connectivity.csv` now carries real p-values
too, but no gene in it survives correction.

See also [[2026-11-09-perm031-at-157-draws-is-a-clean-negative]].
