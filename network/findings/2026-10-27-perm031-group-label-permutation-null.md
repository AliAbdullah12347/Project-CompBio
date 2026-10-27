# perm-031: group-label permutation null for module preservation

**Date:** 2026-10-27
**Analysis:** `network/runs/perm-031` — K=20 permuted draws of 74 random samples from all 474, compared to base-001 bp_lith observed Zsummary
**Evidence:** [FILL] · **Interest:** [FILL] · **Δ:** new

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

[FILL after perm-031 complete.

Key stats:
- null Zsummary: mean and sd per module
- observed mean (base-001 lith) per module
- direction: what fraction of modules are less preserved than null mean?
- empirical p distribution: how many modules with p=0.05 vs p=1.0?
- n_sig BH 5%: expected 0 (see limitation above)

Compare to nolith: does the null Zsummary land near bp_nolith's actual Zsummary?
(If nolith ≈ null, that confirms null is well-calibrated.)
]

## Comparison with gold-module baseline

[FILL: what is the null Zsummary for the gold module itself? If null gold Z < real module
Z, it confirms the internal baseline and external null are consistent.]

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
