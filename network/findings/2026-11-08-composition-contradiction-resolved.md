# Resolved: lithium associations do not survive composition adjustment — finding 001 was right

**Date:** 2026-11-08
**Analysis:** `network/experiments/trait-053-composition-reconcile.R` → `network/runs/trait-053`
**Evidence:** 23/25 · **Interest:** 8/10 · **Δ:** +0.8

## The contradiction

Finding 001: **0 of 33** lithium–eigengene associations survive ILR b1–b4 adjustment.
`grp-051b` (second machine): **28 of 69** survive. Three things differed at once —
module set, sample set, and model — so neither refuted the other.

## The design

Hold the module set **fixed** at the verified 71-module offset-removed reference
(`expr_mode=offset_removed`, gated by `ref_cache_ok()`), and cross the other two.

| sample set | model | modules at BH < 0.05 |
|:--|:--|--:|
| BP1 only (n=226) | `li` (marginal) | 15 / 71 |
| BP1 only | `li + b1..b4` — **001's model** | **1 / 71** |
| BP1 only | `li + dx + technical + b1..b4` — **051b's model** | **1 / 71** |
| all (n=460) | `li` (marginal) | 44 / 71 |
| all | `li + b1..b4` — 001's model | **26 / 71** |
| all | `li + dx + technical + b1..b4` — 051b's model | **1 / 71** |

## The answer

**Finding 001 is correct. One module of 71 survives.** Both routes to controlling
diagnosis agree exactly: restricting to BP1 (001) and adjusting for `dx` (051b) each give
1 of 71.

**The `26 / 71` cell is the trap.** Using all 460 samples while adjusting only for
composition lets the `li` coefficient absorb **diagnosis**, because every lithium user is
a patient and no control is. "Lithium" then means "bipolar and on lithium" versus
"everyone else, mostly controls". Adding `dx` collapses it from 26 to 1.

This is the same confound `grp-051b`'s own comment warns about —

> *Without dx in the model, "lithium" is partly just "bipolar", because every lithium user
> is a patient.*

— and its joint model does include `dx`. So the residual gap between its 28/69 and the
1/71 here is **not** the model and **not** the sample set: it is the **module set**. Its
69 modules come from a group-specific `control_A` network; these 71 come from the fixed
verified reference. Group-specific decomposition is doing the work.

## What this settles

The arm's headline — a lithium-associated co-expression signal that is substantially cell
composition — **holds on the corrected network.** It was first established on the raw
matrix (finding 002), and the raw-matrix bug could have invalidated it. It does not: on
the verified 71-module reference, 15 marginal associations within BP1 fall to 1 after
composition adjustment.

## Limits

- 1 of 71 is not 0 of 71. One module does survive; it is not identified here and deserves
  a look.
- This is eigengene association, not network structure. The structural instruments
  (preservation, hub overlap, partition agreement) are reported separately and mostly null.
- The group-specific-decomposition result is not explained away, only localised. Why a
  `control_A` network yields 28 survivors where the fixed reference yields 1 is now the
  open question — a decomposition that produces more composition-independent lithium
  signal is either finding something real or overfitting its own group.

E=5 R=5 B=4 D=5 C=4
