# The lithium co-expression graph is 11-22% sparser than the non-lithium one, in three independent draws, and it localises to no module

**Date:** 2026-11-10
**Analyses:** `comm-038`, `comm-039`, `cent-040` -> `network/runs/{comm-038,comm-039,cent-040}`
**Evidence:** 18/25 · **Interest:** 7/10 · **Δ:** +0.7

## The observation

Three rows each drew their own 74-sample subsample from the 152 lithium
users, under three different seeds, and built a co-expression graph from it.
All three came out thinner than the 74-sample non-lithium graph:

| graph | n | seed | edges at \|A\| > 0.1 | vs bp_nolith |
|:--|--:|--:|--:|--:|
| control | 234 | - | 1,036,347 | - |
| bp_nolith | 74 | - | 1,118,242 | - |
| bp_lith (`comm-038`) | 74 | 20299010 | 994,895 | **-11%** |
| bp_lith (`comm-039`) | 74 | 20300010 | 951,485 | **-15%** |
| bp_lith (`cent-040`) | 74 | 20301009 | mean degree -22% | **-22%** |

The `comm-038` and `comm-039` subsamples share 33 of 74 samples against a
chance expectation of 36, so these are effectively independent draws rather
than three views of one.

**This is not a sample-size artifact.** All three lithium graphs and the
non-lithium graph have exactly 74 samples. Note that control, with 234, has
*fewer* edges than bp_nolith with 74 - at a fixed threshold, smaller n gives
noisier and therefore larger absolute correlations, so edge count rises as n
falls. That is why only the 74-vs-74 comparison means anything here, and it
is the one that shows the difference.

## It does not localise to any module

This is the same global difference `pres-030` reported from the other
direction - Zdensity 48.3 against Zconnectivity 15.6, with 70 of 71 modules
below the diagonal. The arm's properly-nulled tests look for it at the module
level and do not find it:

- `perm-031` at K=157, against a group-label permutation null: **0 of 71
  modules** at BH 5%, smallest p 0.038.
- **0 of 25 parameter variants** show a lith/nolith preservation difference,
  p 0.441 to 0.989 (`runs/synth-048/primary_contrast_by_variant.csv`).
- `cent-040`'s 52-of-71 degree hits are this one shift counted once per
  module - 46 of the 52 move in the same direction, and betweenness on a
  thresholded graph, the only density-insensitive metric, gives 0 of 71.

So the finding is a global one: **the lithium transcriptome is slightly less
coherently co-expressed, and no module is responsible.**

## What it is not evidence of

Diagnosis is held constant - both groups are bipolar I - so this is not a
case-control effect. But "lithium thins the correlation structure" and "the
lithium group is more heterogeneous" predict the same thing, and nothing here
separates them. Greater spread in cell composition, illness state or
unrecorded co-medication among users would thin a correlation graph exactly
this way, and the deposit has no BMI, no other medications and no lithium
dose (`CLAUDE.md`), so the heterogeneity account cannot be tested on this
data.

The honest claim is the measurement, not the mechanism.

## Community structure, as a separate question

`comm-038` (Louvain) and `comm-039` (Leiden) at resolution 1:

| group | n communities (Louvain / Leiden) | ARI vs WGCNA | Leiden vs Louvain ARI |
|:--|:--|--:|--:|
| control | 3061 / 3063 | 0.148 / 0.123 | 0.499 |
| bp_nolith | 1618 / 1617 | 0.176 / 0.184 | 0.827 |
| bp_lith | 2048 / 1797 | 0.153 / 0.153 | 0.605 |

**The two algorithms agree with each other far better than either agrees with
WGCNA** - ARI 0.50 to 0.83 between them, 0.12 to 0.18 against WGCNA modules,
with NMI steady at 0.47 to 0.52 throughout. Graph community detection on a
thresholded correlation graph finds 1600-3100 groups where WGCNA's TOM-based
clustering finds 71. These are different objects, and a module is not a
community; the agreement that exists is real but weak, and the two should not
be used interchangeably.

Leiden and Louvain agree least on the lithium graph (0.605) and most on the
non-lithium one (0.827), which is what a sparser, noisier graph produces.
That is the same global difference again, a third time.

## Two bugs had to be fixed before any of this could be read

Both rows exited 0 and wrote plausible numbers.

- **`comm-039` called Leiden with igraph's default objective function.**
  `cluster_leiden()` defaults to the Constant Potts Model, where the
  resolution parameter is compared against internal edge *density*. This
  graph's density is 0.0124, so at gamma = 1 no community can qualify and
  every gene becomes its own: 12,368 communities, ARI 0 against everything.
  That was recorded as a result about community structure and was a unit
  mismatch. `objective_function = "modularity"` is the resolution the row was
  pre-specified in and the one `comm-038`'s `cluster_louvain()` already uses.
- **Both rows computed NMI with the subsetting on the wrong side of the
  division.** `pab[pab > 0] / outer(pa, pb)` divides a vector of the non-zero
  cells by the whole matrix, R recycles, and the `[pab > 0]` outside then
  samples the result. Under CPM it happened to land on 0.639 and look
  plausible; under modularity it returned negative numbers, which NMI cannot
  take. Verified after the fix: 1.000 against itself, 0.005 for independent
  partitions, 0.817 for a coarsening.

See also [[2026-11-10-centrality-difference-is-one-global-shift]] and
[[2026-11-09-perm031-at-157-draws-is-a-clean-negative]].
