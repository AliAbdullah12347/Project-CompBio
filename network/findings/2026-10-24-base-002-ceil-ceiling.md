# base-002-ceil: split-half ceiling calibration (re-run with verified offset-removed expression)

**Date:** 2026-10-26
**Analysis:** `network/runs/base-002-ceil` (attempt 4; offset-removed expr_for(); 67-module reference on 160 controls, split-half ceiling)
**Evidence:** 20/25 · **Interest:** 7/10 · **Δ:** new

> **Replaces** the original 2026-10-24 finding, which used raw expression matrices for
> both reference and test (commit 1177d8c, 2026-10-25). All numbers in the original
> finding were invalid. The split-half ceiling DESIGN was correct; only the data were wrong.
> Original archived in `network/runs/_raw_test_contaminated_20261025T193947Z/`.

## Design

Reference built on 160 randomly held-out controls (seed fixed). Three test sets, all n=74:

| test set | perms | comparison |
|:--|--:|:--|
| ctrl_ceiling (74 disjoint controls) | 500 | within-population ceiling: same population as reference |
| bp_nolith | 500 | matched instrument — one reference, two test sets, identical n and perm count |
| bp_lith (20 draws of 74 from 152) | 50 | NOT directly comparable to 500-perm sets (see below) |

This is the only design where ctrl_ceiling and bp_nolith face the **same reference** at the
**same n** with the **same permutation budget**. base-001 uses a 234-control reference that
cannot be compared this way.

## Key numbers

**Reference:** 67 modules (n=160 controls, deepSplit=4, mergeCutHeight=0, power=12, bicor,
signed, offset-removed; max module turquoise=971 genes, median=122).
CACHE_KEY `ctrl160_all_bicor_signed_p12_ds4_mch00_offsetrm`.

**ctrl_ceiling** (500 perms, n=74 disjoint controls — done in attempt 3 with verified script):

| metric | value |
|:--|--:|
| gold Zsummary | **35.78** |
| gold medianRank | **68/68 (last)** |
| real modules above gold (Zsummary) | **15 of 67 (22%)** |
| module Z: median | 25.07 |
| module Z: range | 7.91 – 65.08 |
| cor(log10 size, Zsummary) | **+0.838** |

**bp_nolith** (500 perms, n=74; directly comparable to ctrl_ceiling):

| metric | value |
|:--|--:|
| gold Zsummary | **23.16** |
| gold medianRank | **68/68 (last)** |
| real modules above gold (Zsummary) | **47 of 67 (70%)** |
| module Z: median | 31.05 |
| module Z: range | 9.60 – 66.95 |
| cor(log10 size, Zsummary) | **+0.729** |

Top 5 nolith modules by Zsummary:
brown 66.95, turquoise 62.02, tan 60.22, blue 60.10, midnightblue 59.03.

Weakest nolith modules (below 15 — less preserved than many ctrl_ceiling modules):
thistle1 9.60, lavenderblush3 11.66, plum 12.21, coral2 12.34, mediumorchid 13.24.

**bp_lith draws:** PARTIAL at time of writing (draws started ~04:46Z; see STATE.md for
final count). Not directly comparable to 500-perm ctrl_ceiling and bp_nolith (50 vs 500
perms means Zsummary is on a different scale within the same experiment).

## Main finding: global co-expression coherence is lower in bp_nolith than in held-out controls

The ceiling test establishes a key methodological baseline: when test = same population as
reference (74 held-out controls), gold (a random gene set) achieves Zsummary 35.78 and only
22% of real modules beat it. This is because Zsummary captures global co-expression
coherence — in a same-population test, the random gene set reflects that global coherence
as well as any specific module. medianRank is unaffected: gold ranks last (68/68) even when
its Zsummary is high.

In bp_nolith, gold drops to 23.16 (−35% vs ceiling). Three consequences of this drop:

1. **70% of real modules beat gold by Zsummary** (vs 22% in ceiling). This is NOT because
   BP_nolith is "better preserved" — it is because the global coherence dropped, pulling
   gold down more than specific modules.
2. **medianRank puts gold last in both** (68/68). This is the reliable metric: the control
   network structure is genuinely present in bp_nolith, just with weaker overall coherence.
3. **cor(size, Zsummary) = +0.729 in nolith vs +0.838 in ceiling.** Size dependence
   weakens but does not disappear. The Zsummary anomaly seen in the contaminated run
   (nolith cor ≈ −0.036) was entirely an artifact of the raw-expression bug.

**The 35% gold Z reduction (35.78 → 23.16) is the first valid, instrument-matched signal
that BP_nolith networks are less coherent overall than the control reference.** This is
consistent with finding #13 (BP disease tightens module coupling globally) and finding #22
(PC1 of ME matrix captures BP disease).

However, **Zsummary cannot be used to identify which specific modules are disrupted.**
Module medianRanks under the nolith test are where that question lives; bp_lith draws will
allow comparison at the same n_test (74) but with 50 vs 500 perms.

## Modules of interest

The five weakest modules by nolith Zsummary (thistle1, lavenderblush3, plum, coral2,
mediumorchid) score 9.6–13.2 — below the ceiling arm's 7th-percentile module (steelblue
11.67). These may represent genuinely disrupted modules, but with only nolith and 500 perms
the claim requires further support from pres-030 (density vs connectivity decomposition).

## The anomaly is resolved

The original (contaminated) finding reported nolith gold Z=4.64 vs ceiling gold Z=22.58 —
a 5× reduction attributed to "global expression coherence being much lower in bp_nolith."
The correction note in the anomaly finding (2026-10-25-ceiling-vs-nolith-anomaly.md) stated
the anomaly was a raw-expression artifact.

Clean result: ceiling gold Z=35.78, nolith gold Z=23.16. Ratio = 1.54×, not 5×. The
anomaly finding's open question — "was nolith's cor(size, Z)=−0.036 contamination or real?"
— is answered: it was contamination. Clean nolith has cor(size, Z) = +0.729.

## Caveats

- bp_lith perms (50) ≠ ctrl_ceiling/nolith perms (500): bp_lith Zsummary is on a different
  scale and cannot be directly compared to ctrl_ceiling or bp_nolith.
- 67-module reference (160 controls) differs from base-001's 71-module reference (234
  controls): module names are not directly comparable across arms (WGCNA assigns colours by
  size rank within each network).
- The ceiling ctrl result (gold Z=35.78) was computed in attempt 3 (prior session) with the
  verified fixed script. It is confirmed valid and reused in attempt 4.

E=5 R=3 B=4 D=5 C=3
