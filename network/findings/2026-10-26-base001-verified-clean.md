# base-001 re-run on a verified network: gold ranks last every time, but beats most modules on Zsummary

**Date:** 2026-10-26
**Analysis:** `network/runs/base-001` — first run with reference **and** test both verified offset-removed
**Evidence:** 23/25 · **Interest:** 8/10 · **Δ:** new

## Claim

`base-001` has been re-run after the raw-test-matrix fix. The reference cache was
quarantined, rebuilt with `assert_expr_mode()` checking the matrix, stamped
`expr_mode=offset_removed`, and reused. This is the first preservation result in the arm
where both sides are verified rather than assumed.

Control reference (71 modules), test = 74 bp_lith drawn from 152, 20 draws, 50 perms:

| | Zsummary | medianRank |
|:--|--:|--:|
| gold (random genes), bp_lith, 20 draws | **21.74** (sd 7.05) | **72.0 of 72 — last in every draw** |
| gold, bp_nolith, 500 perms | 14.01 | 72 of 72 |

| 71 real modules | count |
|:--|--:|
| beat gold (Zsummary) in all 20 draws | **17** |
| in ≥18 of 20 | 32 |
| **in fewer than 10 of 20** | **18** |

## The fix changed the answer, and in the unflattering direction

| | contaminated (raw test) | verified |
|:--|--:|--:|
| gold Zsummary | 4.06 (sd 0.93) | **21.74** (sd 7.05) |
| modules beating gold in all 20 draws | 46 of 69 | **17 of 71** |

The raw test matrices had **suppressed gold** and so inflated the apparent separation
between real modules and random genes. The retracted claim "46 of 69 modules beat gold in
all 20 draws" overstated it by roughly threefold. A bug fix that makes your result weaker
is the one worth trusting.

## The two statistics disagree, and the disagreement is the point

Gold has a **high Zsummary** (21.74, above most real modules) and the **worst medianRank**
(72 of 72, in all 20 draws, zero variance). Both are true of the same random gene set.

This is the size effect the arm measured earlier as cor(Zsummary, log10 module size) =
0.879: gold is a large random sample, Zsummary rewards size, and `medianRank` — which
aggregates ranks across the component preservation statistics — does not. Gold's rank of
last in 20 of 20 draws is unambiguous; its Zsummary of 21.74 would, read against
Langfelder's ">10 = strong" band, call random genes strongly preserved.

**Report medianRank, or Zsummary against gold from the same run. Never Zsummary against
the published bands.** That instruction survives from the superseded finding #5 and is
now better supported.

## Candidate disrupted modules

18 modules fall below gold in more than half the draws, several in none:
skyblue2 (Zmean 8.64, beats gold 0/20), coral2 (11.02, 0/20), mediumorchid (11.30, 0/20),
lightcoral (13.96, 0/20), mediumpurple2 (13.08, 1/20), skyblue1 (14.55, 1/20).

These are where the arm's question lives — but they cannot yet be called
lithium-disrupted. Nothing here compares bp_lith against bp_nolith under one reference
and one instrument. `base-002-ceil` is re-running now to provide exactly that.

## Caveat

gold's Zsummary sd is **7.05**, against a median per-module sd of 4.09 — the random
reference is *more* variable across subsamples than the real modules. With 20 draws the
mean is reasonably determined, but any single-draw statement about beating gold is noise.

E=5 R=5 B=4 D=5 C=4
