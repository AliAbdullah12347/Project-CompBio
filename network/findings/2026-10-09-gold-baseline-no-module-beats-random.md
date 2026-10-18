# No module reliably beats a random gene set (degenerate decomposition)

> **SUPERSEDED 2026-10-18.** The headline claim is wrong on a real decomposition.
> Rebuilt on offset-removed expression (71 modules), gold scores Zsummary 14.01 with
> medianRank **72 of 72 - last**, and **63 of 69** real modules beat it. The 1-of-4
> figure below came from the degenerate 4-module partition this finding itself names.
>
> What survives: gold's 14.01 is still inside Langfelder's ">10 = strong" band, so the
> published thresholds are still miscalibrated here, and medianRank discriminates where
> Zsummary does not. See `network/findings/2026-10-18-modules-do-beat-gold.md`.
>
> Text below is left unaltered as the record of what was claimed.


**Date:** 2026-10-09
**Analysis:** `network/runs/base-001`
**Evidence:** 6/25 · **Interest:** 4/10  _(mostly superseded; was 20/25 and 9/10)_ · **Δ:** +1.0

## Claim
WGCNA's `gold` module is by definition a random sample of all network genes. In
control-vs-bp_nolith (n=74, 500 perms) it scores Zsummary 20.26, so Langfelder's bands
(>10 = strong) place noise in the strong band. Across 15 bp_lith draws: blue beats gold
in 9/15, turquoise 3/15, yellow and brown 0/15.

**Every Zsummary this arm reports must be quoted against gold, not the published bands.**

## Important scope limit
Computed on the DEGENERATE mergeCutHeight=0.25 decomposition (3-4 modules, one holding
95-97% of genes). Whether it holds at deepSplit=4 / mergeCutHeight=0 is what
`base-003-remerge` decides.

E=5 R=5 B=3 D=5 C=2. Interest 9: the statistic the arm rests on puts random genes in its
top interpretive band.
