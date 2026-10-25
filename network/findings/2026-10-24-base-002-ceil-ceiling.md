# base-002-ceil: split-half ceiling calibrates n=74 noise floor

**Date:** 2026-10-24
**Analysis:** `network/runs/base-002-ceil` (offset-removed, 36-module reference on 160 controls, split-half ceiling)
**Evidence:** 8/25 · **Interest:** 3/10 · **Δ:** -1.0

## Design

Reference built on 160 randomly held-out controls (seed fixed). Three test sets, all n=74:

| test set | perms | comparison |
|:--|--:|:--|
| ctrl_ceiling (74 disjoint controls) | 500 | within-population upper bound |
| bp_nolith | 500 | comparable to ceiling (same perm count) |
| bp_lith (20 draws of 74 from 152) | 50 | NOT comparable to 500-perm sets (see base-001 finding) |

## Key numbers

**Reference:** 36 modules, largest 12.6% of transcriptome. deepSplit=4, mergeCutHeight=0,
power=12, bicor, signed, offset-removed expression, CACHE_KEY `ctrl160_all_bicor_signed_p12_ds4_mch00_offsetrm`.

**ctrl_ceiling** (500 perms, n=74 disjoint controls):

| metric | value |
|:--|--:|
| gold Zsummary | **22.58** |
| module Z range | 9.44 – 35.46 |
| modules above gold | 27 of 36 |
| modules below gold | 9 (darkmagenta=9.44, darkturquoise=11.91, lightcyan=11.10) |

**bp_nolith** (500 perms, n=74; comparable to ctrl_ceiling):

| metric | value |
|:--|--:|
| gold Zsummary | **4.64** |
| module Z range | 4.00 – 71.75 |
| modules above gold | 35 of 36 |

**bp_lith** (50 perms per draw, 20 draws; NOT comparable to 500-perm sets):

| metric | value |
|:--|--:|
| gold mean Z across draws | 20.98 (sd 1.66) |
| top modules | blue 31.65, green 26.04, yellow 25.16, pink 24.98, turquoise 24.52 |

## Main finding: global expression coherence is much lower in bp_nolith than in controls

The critical comparison is ctrl_ceiling vs. bp_nolith — both use 500 perms, both n=74.
At matched instrument, gold Zsummary is **22.58 in ctrl_ceiling vs. 4.64 in bp_nolith**
(ratio 4.9×). Because gold is a random gene set drawn from the reference, its Zsummary
captures overall co-expression coherence rather than any specific module. A 5× reduction
means the general pattern of pairwise correlations in the control reference network is
far less recognisable in bp_nolith than in a disjoint set of healthy controls.

This is consistent with finding #22 (PC1 of ME matrix captures BP disease) and
finding #13 (BP disease tightens inter-module correlations). The global network fabric
differs between controls and nolith patients — not just at the level of specific modules,
but in the overall texture of co-expression.

## Ceiling calibration for bp_lith interpretation

The ctrl_ceiling gold Z = 22.58 (500 perms) provides the upper bound: what Zsummary looks
like when test set = same population as reference, n=74, high-perm instrument. bp_lith
draws cannot be compared to this directly due to the 50 vs. 500 perm mismatch
(established in base-001 finding). inp-022 and inp-023 test whether residualisation
changes bp_lith preservation.

## Modules below gold in ctrl_ceiling

Nine modules score below gold even in the matched-population (same-control) test:
darkmagenta (9.44), lightcyan (11.10), darkturquoise (11.91), steelblue (11.76),
white (11.49), yellowgreen (11.27), darkolivegreen (12.23), saddlebrown (12.23),
sienna3 (12.60). These modules are ill-defined relative to their own null even within
the same population: they are not coherent enough to test in patients.

## Caveats

- bp_lith perms (50) ≠ ctrl_ceiling/nolith perms (500): no valid across-test comparison
- 36-module reference (160 controls) differs from base-001's 71-module reference (234
  controls); these results are not directly comparable to base-001

## CORRECTION 2026-10-25 — all results from raw expression matrices (INVALID)

Commit 1177d8c (2026-10-25, 15:40 EDT) found that `base-002-ceil.R` line 160 built the
reference from raw `t(d$logtpm[, sids])` rather than `expr_for()`. All test matrices were
also raw. The cache key `ctrl160_all_bicor_signed_p12_ds4_mch00_offsetrm` recorded the
INTENDED mode, not the ACTUAL data.

All numbers in this finding — gold Z, module Z ranges, modules above gold, and all
comparisons between ctrl_ceiling and bp_nolith — are from the invalid raw computation.
The biological interpretation ("global network fabric differs between controls and nolith
patients") is wrong; that difference was a raw-graph artifact (see the anomaly finding).

The split-half ceiling DESIGN is correct and the experiment will re-run. Results archived
to `network/runs/_raw_test_contaminated_20261025T193947Z/`. This finding will be updated
when base-002-ceil completes with offset-removed expression.

E=2 R=2 B=1 D=2 C=1
