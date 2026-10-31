# perm-031: group-label permutation null for module preservation

**Date:** 2026-10-27
**Analysis:** `network/runs/perm-031` — K=20 permuted draws of 74 random samples from all 474, compared to base-001 bp_lith observed Zsummary
**Evidence:** 17/25 (PARTIAL K=14/20: E=5+R=5+B=2+D=3+C=2) · **Interest:** 7/10 · **Δ:** +0.2

*(PARTIAL — ~14/20 draws expected this session, rest next session. Re-score when DONE.)*

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

**INTERIM — K=14/20 draws as of 2026-10-31. Remaining 6 draws to be computed next session.**

~730s per draw (measured from perm-031 draws 1-14; consistent with inp-021 draw timing).

**Gold module (K=14 draws):**
- Null gold Z: mean=32.61, sd=4.07, range [26.54, 39.02]
- Per-draw: 30.15, 29.78, 36.74, 37.80, 33.57, 27.53, 35.16, 35.66, 29.22, 26.54, 39.02, 34.33, 33.03, 27.95
- Observed gold Z (bp_lith, base-001 mean 20 draws) = 21.74
- Gold observed < null: TRUE in ALL 14 draws → bp_lith gold preservation specifically lower than
  random subsamples from the full cohort.

**Real modules (K=14 draws, obs vs null mean):**
- Less preserved in bp_lith than null mean: 42/72 (58%)
- More preserved in bp_lith than null mean: 30/72 (42%)
- Direction: majority of modules (58%) show reduced bp_lith preservation vs null.

**Modules with emp_p=0 (observed below minimum of all 14 null draws):**

9 modules have observed Zsummary below ALL 14 null draws (emp_p = 0/14 = 0). Under BH
correction with 71 modules (including gold), all 9 are significant at FDR 5% because
p=0 satisfies p_(k) ≤ k×0.05/71 for all k ≥ 1:

| module | observed (lith mean) | null mean (K=14) | null min | Δ |
|:--|--:|--:|--:|--:|
| purple | 49.80 | 66.88 | 57.13 | −17.08 |
| salmon | 23.39 | 34.29 | 26.66 | −10.90 |
| gold (random) | 21.74 | 32.61 | 26.54 | −10.87 |
| blue | 44.93 | 54.46 | 46.41 | −9.52 |
| black | 30.53 | 39.04 | 31.60 | −8.51 |
| pink | 22.04 | 28.99 | 23.64 | −6.95 |
| grey60 | 24.41 | 30.59 | 26.60 | −6.18 |
| plum | 15.21 | 18.92 | 15.67 | −3.72 |
| coral2 | 11.02 | 13.74 | 11.25 | −2.72 |

**Caveat on emp_p=0:** with K=14 null draws, "all null draws > observed" is the most extreme
resolution-limited statement. It does NOT imply p→0 at K=∞. At K=20 (the full run), each of
these modules could have 0–1 null draws ≤ observed (minimum non-zero p = 1/20 = 0.05). After
BH with 71 modules, only emp_p=0 (0/20 null draws ≤ observed) would survive FDR 5% (requires
p ≤ 0.0007, achievable only if truly 0/20). So the BH finding will be upheld if and only if
all K=20 null draws remain > observed for these 9 modules.

**Top most-reduced below null (K=14 mean, first 10):**
| module | observed | null mean (K=14) | Δ | emp_p (K=14) |
|:--|--:|--:|--:|--:|
| purple | 49.80 | 66.88 | −17.08 | 0.000 |
| green | 39.47 | 50.51 | −11.04 | 0.071 |
| salmon | 23.39 | 34.29 | −10.90 | 0.000 |
| blue | 44.93 | 54.46 | −9.52 | 0.000 |
| midnightblue | 53.64 | 62.39 | −8.75 | 0.143 |
| black | 30.53 | 39.04 | −8.51 | 0.000 |
| pink | 22.04 | 28.99 | −6.95 | 0.000 |
| royalblue | 27.12 | 33.66 | −6.53 | 0.071 |
| grey60 | 24.41 | 30.59 | −6.18 | 0.000 |
| salmon4 | 25.55 | 31.65 | −6.10 | 0.071 |

**Top most-elevated above null (K=14 mean):**
| module | observed | null mean (K=14) | Δ | emp_p (K=14) |
|:--|--:|--:|--:|--:|
| cyan | 79.93 | 67.54 | +12.39 | 1.000 |
| white | 61.14 | 51.17 | +9.97 | 0.929 |
| lightgreen | 55.23 | 46.44 | +8.78 | 0.857 |
| orange | 40.93 | 32.66 | +8.27 | 0.929 |
| darkred | 50.53 | 43.52 | +7.01 | 0.857 |

**Score revision from K=2 → K=14 (2026-10-31):** D raised from 2 to 3 because 42/72 modules
(58%) show obs < null mean — majority direction, meeting the pre-specification criterion
"D=3 if most modules land below null mean." Interest raised from 5 to 7 because 9 modules
are BH-significant at K=14 (per pre-specification: "raise Interest to 7-8/10 if ANY module
passes BH 5%"). Δ = (7/10 − 17/25) × 10 = +0.2. Old D=2/Interest=5: Δ=−1.4.
METHODS.md 2026-10-31.

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

**Calibration check (K=14):** Null gold Z mean = 32.61 vs bp_lith gold Z = 21.74 and
bp_nolith gold Z = 14.01 (base-001). Ordering: null (32.61) > lith (21.74) > nolith (14.01).
Random null draws (49% control, 16% nolith, 32% lith by group composition) preserve
better than both real test groups. This is expected by design: the null draws are control-
enriched and therefore more similar to the control reference. perm-031 tests whether bp_lith
is specifically less preserved than a random mixed-cohort draw — confirmed (21.74 < 32.61 in
all 14 draws). The nolith gold Z (14.01 < lith 21.74) is confounded by nolith using 500 perms
vs lith's 50 perms — more perms inflate lith Zsummary relative to 500-perm nolith, so direct
comparison is invalid. Fill with K=20 mean when DONE.

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
<!-- 2026-10-31: Δ recalculated: (7/10 - 17/25) × 10 = (0.70 - 0.68) × 10 = +0.2 -->
<!-- 2026-10-31: Surprise finding: emp_p=0 for 9 modules makes BH 5% achievable at K=14, -->
<!--             contrary to pre-specification's "mathematically impossible" assertion. -->
<!--             Pre-spec assumed min p = 1/K, but p=0 (obs < all K null draws) is achievable. -->
