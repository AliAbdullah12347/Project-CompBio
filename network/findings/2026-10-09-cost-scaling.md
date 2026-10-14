# Cost scaling, measured

**Date:** 2026-10-09
**Analysis:** Phase 0 probes + base-001 draw timings
**Evidence:** 6/25 · **Interest:** 1/10 · **Δ:** -1.4

Reference build 786 s x (n_ref/74); `modulePreservation` 5.2 s per permutation; one
50-permutation draw 262 s (sd 3 s). There is no per-draw reference rebuild - an earlier
1962 s/draw estimate was 7.5x too high and drove an over-conservative design choice
(100 draws reduced to 20).

E=4 R=0 B=1 D=1 C=0. Infrastructure, scored as such.
