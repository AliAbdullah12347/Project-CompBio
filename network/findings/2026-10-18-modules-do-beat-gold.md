# Real modules do beat the random gene set — once the decomposition is real

> **CORRECTION 2026-10-26 — numbers updated from base-001 verified re-run.**
>
> The nolith arm figures in this finding came from the first base-001 run (2026-10-18).
> The verified re-run (2026-10-26, committed after raw-test-matrix fix) gives:
>
> - "63 of 69 real modules above gold" → **65 of 71** (same gold Z=14.01, same rank 72/72)
> - "sd 0.93" and "sd 1.86 at 15 draws" are from contaminated lith draws (raw test matrices).
>   Verified lith gold Z=21.74, **sd 7.05** (20 draws). See
>   `findings/2026-10-26-base001-verified-clean.md`.
>
> The main claim — real modules decisively beat random once the decomposition is real —
> stands and is strengthened. The nolith bp arm gives gold rank 72/72 (last) in a 71-module
> network; gold also ranks last in all 20 lith draws.

**Date:** 2026-10-18
**Analysis:** `network/runs/base-001` (offset-removed, 71-module reference)
**Evidence:** 21/25 · **Interest:** 9/10 · **Δ:** new

## Claim

Finding #5 reported that no module reliably beats WGCNA's `gold` module (a random sample
of network genes), and concluded the decomposition was degenerate. On the rebuilt network
the result reverses, decisively.

Control reference vs bp_nolith, n = 74, 500 permutations — the identical comparison:

| | raw matrix (4-module) | offset removed (71-module) |
|:--|--:|--:|
| gold Zsummary | 20.26 | **14.01** |
| gold **medianRank** | 4 of 4 | **72 of 72 (last)** |
| real modules above gold | **1 of 4** | **63 of 69** |
| real-module Zsummary median | — | 31.21 |
| real-module Zsummary range | — | 8.34 – 75.31 |

Top preserved: cyan 75.31, midnightblue 59.78, yellow 58.63, purple 58.60, white 57.59.
Weakest: mediumorchid 8.34, skyblue2 10.53, mediumpurple2 10.61.

## What survives of finding #5, and what does not

**Does not survive: "no module beats random."** That was a property of a decomposition in
which 97.5% of genes sat in one module. Given a real partition, random genes rank *last*
of 72 and are beaten by 63 of 69 real modules. The original finding's own title called the
decomposition degenerate — it was right about that, and wrong to generalise from it.

**Survives, in weakened form: the published bands are miscalibrated on this data.** Gold
still scores **14.01**, comfortably inside Langfelder's ">10 = strongly preserved" band.
A random gene set should not be "strongly preserved" under any threshold. So the original
instruction — *never quote a Zsummary against the published bands alone* — still holds.

**The refinement is the useful part: `medianRank` discriminates where `Zsummary` does
not.** Gold's Zsummary (14.01) looks respectable in isolation; its medianRank (72 of 72)
is unambiguous. Zsummary appears inflated for everything on this data, so rank-based
comparison against gold is the statistic to report.

## Caveat

This is the bp_nolith arm only — one fixed test set, 500 permutations, no subsampling.
The 20 bp_lith subsample draws are still running. The earlier over-claim on this exact
statistic (5 draws called "STABLE" at sd 0.93, which became sd 1.86 at 15 draws) is a
standing reason not to generalise from a single arm. Treat the 63-of-69 figure as
provisional until the lith draws complete.

E=5 R=4 B=4 D=5 C=3
