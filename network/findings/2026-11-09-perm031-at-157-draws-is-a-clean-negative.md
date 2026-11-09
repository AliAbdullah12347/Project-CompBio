# perm-031 at 157 draws: the K=20 lead was a floor artifact, and there is no module-specific lithium effect

**Date:** 2026-11-09
**Analysis:** `network/experiments/perm-031.R` → `network/runs/perm-031` (157 draws)
**Evidence:** 22/25 · **Interest:** 8/10 · **Δ:** +0.9

## The result

**Zero of 71 modules** are less preserved in the lithium network than the
group-label permutation null at BH 5%. The smallest empirical p is **0.038**,
which is not significant even before correction.

The four modules that sat below all 20 draws in the K=20 run, and were the
reason for raising K at all:

| module | obs Z | null Z | p | BH |
|:--|--:|--:|--:|--:|
| grey60 | 24.41 | 31.10 | 0.0380 | 0.820 |
| blue | 44.93 | 53.35 | 0.0506 | 0.820 |
| purple | 49.80 | 63.86 | 0.0570 | 0.820 |
| pink | 22.04 | 26.99 | 0.0759 | 0.832 |

**No module is below every null draw any more.** At K=20 these four were, which
is what "p = 0.05" meant; with 157 draws the null distribution fills in around
them and they land in the ordinary tail.

## Two separate things were wrong at K=20

**The budget could not produce an answer.** With `m` modules tied at the floor
`1/(K+1)` among `n = 71`, the BH value for the last of them is
`(1/(K+1)) · n/m`, so a positive needs `K > 20n/m − 1`:

| modules at the floor | draws needed | K=20 | K=157 |
|--:|--:|:--|:--|
| 4 | 354 | no | no |
| 8 | 177 | no | no |
| 12 | 118 | no | yes |

At K=20 nothing could clear BH whatever the biology. **This calculation should
have been done before the run, not after** — it decides whether the experiment
is worth doing at all.

**The estimator was invalid.** Both `perm-031` and `perm-032` computed
`mean(nz <= obs)`, which returns exactly 0 when the observed value is below
every null draw. A p of 0 is not a probability; it says the true p is below
1/K, not that it is zero. `perm-032` reported **"8 modules significant at BH
5%"** on the back of eight such zeros. With the add-one estimator
`(1 + sum(nz <= obs)) / (1 + K)` the same draws give **0 significant modules**,
because the smallest attainable BH value at K=20 across 72 modules is 0.429.

Both scripts now use add-one. `perm-032`'s corrected table is in
`results_permutation_test_addone.csv`; **its published count of 8 should not be
cited.**

## What is left over

Every one of the top eight is in the same direction — observed below null —
and they bunch at p ≈ 0.04–0.09 rather than spreading uniformly. `gold`, the
all-genes pseudo-module, is among them (p = 0.095, obs 21.7 vs null 32.3).

That pattern is a weak *global* tendency for the lithium network to be slightly
less preserved, not a set of specific modules. It is what `pres-030` already
reported from a different angle: Zdensity 48.3 against Zconnectivity 15.6, with
70 of 71 modules below the diagonal. A global shift of that size is what
subsampling 152 lithium users down to 74 produces, and the matched ceiling in
`base-002-ceil` is the control for it — ceiling medianRank 34.0 vs bp_nolith
33.0, Wilcoxon p = 0.724.

## How this reads in the paper

This is a negative result that is now worth stating, because the test has the
resolution to detect something and does not. Combined with the 21 parameter
variants that all agree (`mod-013`…`mod-020`, `pow-010`…`pow-012`,
`corr-005/006`, `net-007`, `gene-025`…`gene-029`, `ref-002`: **0 of 21** show
any lithium/non-lithium preservation difference, p = 0.44–0.99), the network
arm's conclusion is stable under every parameter we varied.

## Caveat

157 of the intended 200 draws. The 43 gaps are the slices that hit the cluster's
24 h ceiling, so which draws are missing depends on pod scheduling and not on
any outcome — the draws are exchangeable and a partial set is a valid null. It
raises the floor to 1/158 = 0.0063, which does not change the conclusion: not
one module reached even the unadjusted 0.05.

See also [[2026-11-09-eosinophil-module-survives-composition]].
