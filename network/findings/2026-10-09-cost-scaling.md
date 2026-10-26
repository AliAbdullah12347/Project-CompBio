# Cost scaling, measured

**Date:** 2026-10-09 (updated 2026-10-26)
**Analysis:** Phase 0 probes + base-001 draw timings + gene-027 nolith timing
**Evidence:** 6/25 · **Interest:** 1/10 · **Δ:** -1.4

Reference build 786 s (n_ref=234, 13616 genes, bicor signed); `modulePreservation` ~12.8 s per
permutation at n_ref=234, n_test=74, 13616 genes; one 50-permutation draw ~720 s (setup ~80 s +
50×12.8 s). The earlier 262 s/draw and 5.2 s/perm figures were from Phase 0 probes on a smaller
gene set and should not be used for budgeting.

**Critical discovery (2026-10-26):** The `modulePreservation()` "calculating observed preservation
values" step scales with `nPermutations` in WGCNA 1.74. Empirically: ~80 s at nPerm=50 vs
~6460 s at nPerm=500 for a 65-module reference (12018 genes, n_ref=234, n_test=234). This is a
pre-permutation setup step — it runs once per call at full cost regardless of draw count. 500
nolith perms alone would require ~12900 s, exceeding any practical session cap. All PENDING
experiments were revised from NOLITH_PERMS=500 to NOLITH_PERMS=50 on 2026-10-26 before any
nolith results were computed (pre-specification integrity maintained; see METHODS.md 2026-10-26).

There is no per-draw reference rebuild. An earlier 1962 s/draw estimate was 7.5× too high and
drove an over-conservative design choice (100 draws reduced to 20).

E=4 R=0 B=1 D=1 C=0. Infrastructure, scored as such.
