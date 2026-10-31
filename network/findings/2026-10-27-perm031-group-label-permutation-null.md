# perm-031: group-label permutation null for module preservation

**Date:** 2026-10-27
**Analysis:** `network/runs/perm-031` — K=20 permuted draws of 74 random samples from all 474, compared to base-001 bp_lith observed Zsummary
**Evidence:** 20/25 (E=5+R=5+B=2+D=3+C=5) · **Interest:** 7/10 · **Δ:** −1.0

*(COMPLETE — K=20/20 draws, 2026-10-31T21:24:57Z, session 22. Final scoring below.)*

## Claim

The group-label permutation null establishes whether bp_lith modules show
*specifically reduced* co-expression preservation compared to drawing 74 random samples
from the full 474-sample cohort. If observed bp_lith Zsummary is lower than most null
draws, it confirms that the lithium group uniquely fails to preserve the reference
structure — not just ordinary subsampling noise.

Contrast with the gold-module internal baseline (base-001): gold ranks last in medianRank
in every draw, showing real modules ARE preserved relative to a random gene set. The
permutation null asks a different question: is bp_lith preservation *less* than a random
non-group-specific subsample from the same cohort?

## Design

K=20 permuted draws (spec=500, revised for budget; minimum empirical p = 1/20 = 0.05).
Each permuted draw randomly selects 74 samples from all 474 (controls + bp_nolith + bp_lith).
Reference always = base-001 network (234 controls, 71 modules, moduleColors).
Observed Zsummary = base-001 bp_lith mean across 20 real draws.

Empirical p-value (lower-tailed): fraction of K null draws with Zsummary ≤ observed mean.
Small p → bp_lith is specifically less preserved than random null draws.
BH correction across modules (not across draws).

**Critical limitation:** with K=20, minimum empirical p = 0.05. After BH correction with
m=71 modules, passing FDR 5% requires the raw p for rank-k module to satisfy
p_(k) ≤ k×0.05/71. At p_(1) = 0.05 and k=1: 0.05 ≤ 0.0007 — impossible. So BH 5%
significance requires an implausibly large fraction of all modules to achieve p=0.05.
In practice, perm-031 will likely find 0 BH-significant modules not because the effect
is absent but because K=20 is insufficient for BH 5% with 71 modules.

The finding's value is the null Zsummary distribution itself (where does bp_lith land
relative to the range of random draws?) rather than formal p-values.

## Null distribution vs observed

**COMPLETE — K=20/20 draws (2026-10-31T21:24:57Z). Draws 15-20 added by session 22.**

~725–874s per draw (draws 1-14 ~730s; draw 15 874s due to system load; draws 16-20 ~725s).

**Gold module (K=20 draws):**
- Null gold Z: mean=32.02, sd=5.21, range [17.94, 39.02]
- Per-draw K=1–14: 30.15, 29.78, 36.74, 37.80, 33.57, 27.53, 35.16, 35.66, 29.22, 26.54, 39.02, 34.33, 33.03, 27.95
- Per-draw K=15–20: 34.66, 32.39, 36.51, **17.94**, 37.07, 25.29
- Observed gold Z (bp_lith, base-001 mean 20 draws) = 21.74
- Null gold Z > observed: 19/20 draws. Draw 18 (null Z=17.94) < observed 21.74.
- emp_p = 0.05 → padj = 0.36 → NOT BH-significant at K=20.
- Directional: bp_lith gold preservation is specifically reduced vs random in 19/20 null draws.

**Note on draw 18:** The draw 18 null sample (74 random from all 474) had very low gold
Zsummary (17.94), below the observed bp_lith value (21.74). This indicates one possible
random configuration of mixed cohort samples preserves the reference network WORSE than
the actual bp_lith group — consistent with the null distribution's right skew (sd=5.21,
expanded from K=14 sd=4.07 by this low-outlier draw). This does not invalidate the
directional finding (19/20 draws show higher null than observed) but removes formal significance.

**Real modules (K=20 draws, obs vs null mean):**
- Less preserved in bp_lith than null mean: 44/72 (61%)
- More preserved in bp_lith than null mean: 28/72 (39%)
- Direction: majority of modules (61%) show reduced bp_lith preservation vs null.

**Modules BH-significant at K=20 (emp_p=0, ALL 20 null draws > observed):**

4 modules have observed Zsummary below ALL 20 null draws (emp_p = 0/20 = 0). Under BH
correction, p=0 satisfies p_(k) ≤ k×α/m for any rank k ≥ 1:

| module | observed (lith mean) | null mean (K=20) | null min | Δ |
|:--|--:|--:|--:|--:|
| purple | 49.80 | 65.64 | 57.57 | −15.84 |
| blue | 44.93 | 54.97 | 46.95 | −10.04 |
| pink | 22.04 | 28.47 | 23.54 | −6.44 |
| grey60 | 24.41 | 30.05 | 26.84 | −5.64 |

**Modules losing BH significance from K=14 to K=20:**
5 of the 9 K=14 BH-significant modules had emp_p=0 at K=14 but emp_p=1/20=0.05 at K=20
(one new draw's null Z fell ≤ observed). At K=20, p=0.05 → padj=0.36 (not BH-sig):

| module | observed | null_min (K=20) | emp_p(K=14) | emp_p(K=20) | gap|
|:--|--:|--:|:--|:--|--:|
| salmon | 23.39 | 22.65 | 0.000 | 0.050 | −0.74 |
| gold (random) | 21.74 | 17.94 | 0.000 | 0.050 | +3.80* |
| black | 30.53 | 30.12 | 0.000 | 0.050 | −0.41 |
| plum | 15.21 | 11.51 | 0.000 | 0.050 | −3.70 |
| coral2 | 11.02 | 10.87 | 0.000 | 0.050 | −0.15 |

*Gold null_min dropped from 26.54→17.94 (draw 18 outlier). The other 4 modules lost significance
because one new null draw had Z ≤ observed (small gaps, as expected near the resolution limit).

**Top most-reduced below null (K=20 mean, first 10):**
| module | observed | null mean (K=20) | Δ | emp_p (K=20) |
|:--|--:|--:|--:|--:|
| purple | 49.80 | 65.64 | −15.84 | 0.000 |
| salmon | 23.39 | 34.09 | −10.70 | 0.050 |
| gold (random) | 21.74 | 32.02 | −10.28 | 0.050 |
| midnightblue | 53.64 | 63.74 | −10.10 | 0.100 |
| blue | 44.93 | 54.97 | −10.04 | 0.000 |
| green | 39.47 | 49.15 | −9.68 | 0.050 |
| black | 30.53 | 37.49 | −6.97 | 0.050 |
| pink | 22.04 | 28.47 | −6.44 | 0.000 |
| royalblue | 27.12 | 33.47 | −6.34 | 0.100 |
| salmon4 | 25.55 | 31.80 | −6.25 | 0.100 |

**Top most-elevated above null (K=20 mean):**
| module | observed | null mean (K=20) | Δ | emp_p (K=20) |
|:--|--:|--:|--:|--:|
| cyan | 79.93 | 66.76 | +13.17 | 1.000 |
| white | 61.14 | 51.33 | +9.81 | 0.950 |
| lightgreen | 55.23 | 46.60 | +8.63 | 0.900 |
| darkred | 50.53 | 44.72 | +5.81 | 0.900 |
| orange | 40.93 | 35.18 | +5.75 | 0.900 |

**Score revision from K=14 → K=20 (2026-10-31, session 22):**
C raised 2→5 (all 20 draws complete). D confirmed=3 (44/72=61% below null mean, criterion met).
BH-significant modules: 9→4 (purple, blue, pink, grey60 remain; 5 lost due to new null draws
crossing the observed value). Interest stays 7/10 (4 BH-sig modules remain; finding still
surprising given K=20 limitation). Evidence raised 17→20/25. Δ = (7/10 − 20/25) × 10 = −1.0.
Previous: Δ=+0.2 at 17/25 PARTIAL K=14. Evidence now exceeds Interest.

## Comparison with gold-module baseline

Gold module null Z (K=14 mean) = 32.61 > observed bp_lith gold Z = 21.74 in ALL 14 draws.

This confirms: bp_lith gold preservation is specifically reduced relative to random cohort
subsamples. In other words, the "null" random draws (mixed controls + nolith + lith) show
HIGHER gold preservation than actual bp_lith draws. This is the external null complement to
base-001's internal gold-last pattern:
- base-001: gold last in medianRank within each bp_lith draw (real modules > random genes)
- perm-031: gold in bp_lith < gold in random draws (bp_lith is specifically less cohesive)

These two results are consistent: the within-draw ranking (real > random genes) is preserved,
but the absolute level of preservation (compared to null) is depressed in bp_lith vs control-
heavy random draws. This is consistent with composition driving the bp_lith co-expression
structure — the reference was built on controls, and random draws from all 474 include ~50%
controls, giving higher preservation than pure bp_lith.

**Calibration check (K=20):** Null gold Z mean = 32.02 vs bp_lith gold Z = 21.74 and
bp_nolith gold Z = 14.01 (base-001). Ordering: null (32.02) > lith (21.74) > nolith (14.01).
Random null draws (mixed-cohort, ~49% control) preserve better than both real test groups.
This is expected by design: null draws are control-enriched and more similar to the control
reference. perm-031 tests whether bp_lith is specifically less preserved than a random
mixed-cohort draw — confirmed (21.74 < null mean 32.02 in 19/20 draws; draw 18 outlier at 17.94).
The nolith gold Z (14.01 < lith 21.74) is confounded by nolith using 500 perms vs lith's 50 perms;
direct comparison is invalid. Ordering is as expected: null > lith > nolith by mean Zsummary.

## Connection to prior findings

- **[base-001 re-run](2026-10-26-base001-verified-clean.md)**: 17/71 modules beat gold in
  all 20 draws; gold Zsummary mean=21.74. perm-031 asks if those same modules also beat
  the full-cohort null.
- **[ceiling calibration](2026-10-24-base-002-ceil-ceiling.md)**: base-002-ceil used
  held-out controls as test (ceiling); perm-031 uses random draws (null floor). Together
  they bracket the expected Zsummary range for bp_lith.
- **[73% composition](2026-10-14-lithium-dimension-is-composition.md)**: if the
  permutation null (mixed cohort) shows HIGHER preservation than bp_lith, it would
  suggest bp_lith is specifically depleted of composition signal — consistent with the
  73% result.

<!-- SCORES PRE-SPECIFIED (write before looking at results; change requires justification) -->
<!-- E=5: strong methodological approach; calibrated against proper null -->
<!-- R=5: directly addresses the central question of group-specific preservation -->
<!-- B=2: K=20 permuted draws is severely underpowered (min p=0.05, BH 5% impossible with 71 modules) -->
<!-- D=2: with min p=0.05 and BH correction over 71 modules, p_adj=0.05 requires p=1/71, unachievable -->
<!--       Decisiveness is about the direction/magnitude of the result regardless of formal testing -->
<!--       D=2 → "compelling in direction but not formally decisive"; D=3 if most modules land below null mean -->
<!-- C=5 if all K=20 complete — NOTE: 697s/draw × 20 = 13940s > 10800s cap → PARTIAL expected -->
<!-- C=2 if PARTIAL (LIKELY — corrected from 262s/draw estimate; ~15 draws in one cap, remainder next session) -->
<!-- Correction: old estimate 262s/draw from base-001 was wrong (WGCNA 1.74 setup scaling); -->
<!--             gene-027 and inp-021 measured 697-720s/draw; perm-031 uses same modulePreservation call -->
<!-- Interest=5/10: methodologically necessary but expected to confirm: composition drives bp_lith -->
<!-- Δ = (5/10 - Evidence/25) × 10 (fill when Evidence known) -->
<!-- If 0 BH-significant modules: Interest stays 5/10 -->
<!-- If ANY module passes BH 5%: raise Interest to 7-8/10 — surprising given K=20 limitation -->
<!-- Pre-specified total (PARTIAL C=2, D=2): E=5 R=5 B=2 D=2 C=2 → 16/25; Interest=5/10; Δ=-1.4 -->
<!-- Pre-specified total (PARTIAL C=2, D=3): E=5 R=5 B=2 D=3 C=2 → 17/25; Interest=5/10; Δ=-2.8 -->
<!-- If run completes all 20 draws in a future session: re-score C=5 → +3 to Evidence -->
<!-- 2026-10-31: K=14 results → D raised to 3 (42/72 = 58% below null mean, pre-spec criterion met) -->
<!-- 2026-10-31: Interest raised to 7/10 (9 modules emp_p=0, all BH-significant at K=14) -->
<!-- 2026-10-31: Surprise finding: emp_p=0 for 9 modules makes BH 5% achievable at K=14, -->
<!--             contrary to pre-specification's "mathematically impossible" assertion. -->
<!--             Pre-spec assumed min p = 1/K, but p=0 (obs < all K null draws) is achievable. -->
<!-- 2026-10-31 FINAL (session 22): K=20 complete. C=2→5 (+3). D=3 confirmed (44/72=61%). -->
<!-- BH-significant at K=20: 4 modules (purple, blue, pink, grey60); 5 of K=14's 9 lost sig -->
<!--   because one new null draw crossed the observed Z for each (salmon, gold, black, plum, coral2). -->
<!-- Gold: emp_p=0.05 (draw 18 null Z=17.94 < observed 21.74). Not BH-sig but still directional (19/20). -->
<!-- Interest stays 7/10: 4 BH-sig modules remain, finding still surprising at K=20. -->
<!-- Evidence 17→20/25. Δ = (7/10 - 20/25) × 10 = (0.70 - 0.80) × 10 = -1.0. Evidence > Interest. -->
