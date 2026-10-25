# The one fully-matched comparison behaves incoherently — stop before interpreting it

**Date:** 2026-10-25
**Analysis:** `network/runs/base-002-ceil`
**Evidence:** 19/25 · **Interest:** 8/10 · **Δ:** new

## The comparison finally exists

`base-002-ceil` scores **one** reference — the 160-control, 36-module ctrl160 network —
against three test sets. Two of them are fully matched:

| test set | n | perms | reference |
|:--|--:|--:|:--|
| 74 held-out controls (`ctrl_ceiling`) | 74 | 500 | ctrl160 |
| bp_nolith (`nolith`) | 74 | 500 | ctrl160 |

Same reference, same gene sets, same instrument. This is exactly the design the
2026-10-25 synthesis finding said did not exist — **that claim was wrong and is corrected
here.** No new experiment is needed to make the comparison; what is needed is to trust it.

## And it does not behave

| | ceiling (controls) | bp_nolith |
|:--|--:|--:|
| gold (random genes) Zsummary | 22.58 | **4.64** (−79%) |
| real modules, median Zsummary | 18.69 | **36.81** (+97%) |
| modules lower in bp_nolith | — | **4 of 36** |
| **cor(log10 module size, Zsummary)** | **+0.911** | **−0.036** |
| Spearman(Z_ceiling, Z_nolith) across modules | — | **+0.072** |
| medianRank, paired Wilcoxon | — | **p = 0.76** |

Three things here are not jointly believable:

1. **Real modules preserve *better* in patients than in held-out controls** (18.69 →
   36.81, 32 of 36 rising). A split-half ceiling is supposed to be the upper bound. A
   disease group exceeding it by ~2× is not a biological result, it is a symptom.

2. **Gold moves the opposite way to the modules** (−79% while modules go +97%) in the
   same comparison, on the same reference, with the same permutation budget.

3. **Zsummary is essentially module size in one arm and unrelated to it in the other**
   (+0.911 vs −0.036). This arm already knows Zsummary tracks size — it recorded
   cor(Zsummary, log10 size) = 0.879 earlier. Seeing that dependence vanish between two
   runs of the same statistic points at the computation, not the biology.

The only statistic that behaves consistently is the ordinal one: gold ranks 36 then 37
of 37 (worst in both), and across real modules the paired Wilcoxon on `medianRank` gives
**p = 0.76 — no detectable difference** between held-out controls and bipolar patients.
But Spearman between the two arms' rankings is only +0.312, so even the orderings
disagree more than they should.

## What to do before any of the remaining 45 rows is interpreted

Diagnose the extraction, not the biology. `base-002-ceil.R:184-188` pulls
`mp$preservation$Z[[key_obs]][[key_test]]` by constructed name
(`inColumnsAlsoPresentIn.<test>`). The first thing to check is whether the ceiling and
nolith blocks are being read from the slots they are assumed to be, and whether
`modulePreservation` was given the same `networkType`/`corFnc` for both.

Until that is settled, **no preservation comparison in this arm should be quoted** —
including the within-arm gold-anchored ones, which survive scale problems but not a
mis-indexed Z block.

## Status of the claims this displaces

- The 2026-10-25 synthesis finding's "the experiment does not exist yet" is wrong;
  corrected there.
- The 2026-10-24 inp-022 "rescue" finding is separately corrected: its cross-arm
  comparison used two different reference networks.

E=4 R=5 B=3 D=5 C=2
