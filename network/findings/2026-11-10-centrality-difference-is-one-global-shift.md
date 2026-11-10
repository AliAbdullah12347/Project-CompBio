# cent-040's "52 of 71 modules differ" is one global density shift counted 52 times

**Date:** 2026-11-10
**Analysis:** `network/experiments/cent-040.R` -> `network/runs/cent-040`
**Evidence:** 14/25 · **Interest:** 6/10 · **Δ:** +0.4

## What the row reports, and why it reads wrong

cent-040 compares degree, eigenvector and betweenness centrality between
groups. As it ran on the cluster it tested only `control` against each patient
group, and those two columns look like a result:

| metric | ctrl vs lith | ctrl vs nolith |
|:--|--:|--:|
| degree | 50/71 | 51/71 |
| eigenvector | 70/71 | 66/71 |
| betweenness | 20/71 | 27/71 |

Both columns are nearly the same size, which is the first clue: whatever makes
the lithium graph differ from control also makes the non-lithium graph differ
from control by about as much. Neither column is about lithium. They are also
confounded outright - control has 234 samples and both patient graphs have 74,
and centrality on a correlation graph rises with sample size.

**The matched comparison was missing.** `bp_nolith` and `bp_lith` are both
n=74, so they are the only pair in this row that can be compared without
sample size in the way. It has been added and the row re-run from its
checkpoints:

| metric | nolith vs lith (BH 5%) | direction of the significant ones |
|:--|--:|:--|
| degree | **52/71** | 46 lower in lithium, 6 higher |
| eigenvector | **58/71** | 40 higher, 18 lower |
| betweenness | **0/71** | - |

## Why this is not 52 findings

**The degree result is one number.** Mean weighted degree over all 12,368
genes is 87.6 in `bp_nolith` and 68.7 in `bp_lith` - the whole graph is 22%
less connected. 46 of the 52 significant modules move in that one direction.
A Wilcoxon over the ~53 genes of a module treats them as independent
observations when they all come from a single correlation matrix and all move
together when its overall level moves, so a global shift of this size is
re-detected in every module large enough to reach significance. The test
measures module size as much as it measures biology.

**Eigenvector centrality is normalised to a maximum of 1**, so it reports how
the graph's connectivity is distributed rather than how much there is. When
total degree falls, the distribution re-spreads; mean eigenvector goes the
other way (0.109 nolith, 0.122 lith). That is renormalisation, not a second
independent signal.

**Betweenness is the one metric computed on a thresholded graph** (`A > 0.05`),
which strips the global level out by construction. It shows **0 of 71**. The
metric insensitive to the global shift sees nothing at all.

## The global shift is real; the per-module counting is not

`bp_lith` here is **one** fixed draw of 74 from the 152 users, seed 20301009,
because betweenness is O(n^3) and 20 draws were unaffordable. On its own that
would be too thin to stand on. It does not stand alone: `comm-038` (seed
20299010) and `comm-039` (seed 20300010) drew their own subsamples, sharing
33 of 74 samples with this one against a chance expectation of 36, and both
built a graph 11% and 15% thinner than the 74-sample non-lithium graph. Three
independent draws, same direction, 11 to 22% - see
`findings/2026-11-10-lithium-graph-is-sparser-in-three-independent-draws.md`.

So the global difference is reproducible. **What is not supported is reading
it off as 52 modules.** That number is the one shift re-detected in every
module large enough to reach significance, and the test that would show this
is a control-split null: cut the 234 controls into two disjoint halves of 74
and run the same comparison. If two halves of one group produce a comparable
count, the procedure is measuring the global level and module size, not
biology. That has not been run - two adjacency builds, about an hour - and
until it is, the `padj_nolith_lith` column stays descriptive.

## How this reads in the paper

It does not contradict the arm's negative. `perm-031` at K=157 tests the
module-specific question properly, with a group-label permutation null, and
returns 0 of 71 at BH 5%. `pres-030` already reported the same global tendency
from preservation statistics - Zdensity 48.3 against Zconnectivity 15.6, 70 of
71 modules below the diagonal. cent-040's matched comparison is a third view
of that one global difference, measured with a test that cannot separate it
from module-specific rewiring.

The usable output of this row is `centrality_control/nolith/lith.csv`:
per-gene degree, eigenvector and betweenness in each group, which is a fine
descriptive table for picking hub genes. The `padj_nolith_lith` column should
not be cited as a count of differing modules.

See also [[2026-11-10-lithium-graph-is-sparser-in-three-independent-draws]]
and [[2026-11-09-perm031-at-157-draws-is-a-clean-negative]].
