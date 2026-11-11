# There is no 28-versus-1 discrepancy: both decompositions find exactly one module, and it is the same one

**Date:** 2026-11-11
**Analysis:** `network/experiments/recon-055.R` -> `network/runs/recon-055`
**Evidence:** 23/25 · **Interest:** 7/10 · **Δ:** +0.6

## The question, and why it dissolves

`NEXT_STEPS` item 3 asked why a group-specific decomposition yields **28**
lithium associations surviving composition adjustment where the fixed
verified reference yields **1**. `trait-053` had already ruled out the model
and the sample set, leaving the module set as the only candidate, and that
would have been an interesting finding: a decomposition built on one group
seeing composition-independent signal a fixed reference misses.

It is not what happened. Running `grp-051`'s own `control_A` partition
through `trait-053`'s model, with the expression matrix, sample sets, models
and BH correction all unchanged, gives:

| sample set | model | fixed reference | grp-051 control_A |
|:--|:--|--:|--:|
| BP1 only | `li + b1..b4` | 1 / 71 | **1 / 69** |
| BP1 only | `li + dx + technical + b1..b4` | 1 / 71 | **1 / 69** |
| all | `li + b1..b4` | 26 / 71 | 33 / 69 |
| all | `li + dx + technical + b1..b4` | 1 / 71 | **1 / 69** |

**One, not 28.** And `grp-051b`'s own output file says the same thing. In
`runs/grp-051/module_trait_joint.csv`:

- `li_padj` — lithium adjusted for diagnosis and technical covariates but
  **not** composition: **29 of 69**
- `li_adj_padj` — the same model **plus the four ILR balances**: **1 of 69**

The 28 in `NEXT_STEPS` was the unadjusted column, read as though it were the
adjusted one. The row it came from has been reporting 1 all along.

## What the two partitions actually agree on

The single survivor under the fixed reference is `mediumpurple2`. Under
`grp-051`'s partition it is `orangered3`. Those are the same 35-ish genes:
Jaccard 0.846 against a background of 0.184, established in
[[2026-11-09-eosinophil-module-survives-composition]].

So two decompositions that agree on almost nothing else — built on different
sample sets, yielding 71 and 69 modules, with no shared naming — each find
exactly one module whose lithium association survives adjustment for cell
composition, and it is the same module. That is a stronger statement about
the eosinophil result than either row made on its own, and it arrives as a
by-product of chasing a discrepancy that was not there.

## The partition does matter, just not where it was looked for

Under the weaker model — `li + b1..b4`, no diagnosis and no technical
covariates — the two partitions give 26 and 33 on the full sample set. That
gap is real and is about the module set. It closes completely once `dx`, age,
sex, RIN and plate enter: both drop to 1. So the partition's influence is
entirely absorbed by the covariates the full model already carries, which is
the reassuring direction for the arm's conclusion.

## The lesson worth keeping

A count read from the wrong column of a results file survived in
`NEXT_STEPS` long enough to be scheduled as cluster work, and it was only
caught by running the comparison it asked for. Three other numbers in this
arm failed the same way in the last two days — `hub-037`'s "0 of 12360" came
from an all-NA column, `comm-039`'s "ARI = 0" from a unit mismatch,
`perm-032`'s "8 significant" from an invalid estimator. The common shape is
a plausible number with nothing checking it against the thing it claims to
measure.

See also [[2026-11-09-eosinophil-module-survives-composition]].
