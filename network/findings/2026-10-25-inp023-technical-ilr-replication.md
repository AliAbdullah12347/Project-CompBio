# inp-023: technical + ILR residualisation replicates inp-022's rescue

> **PARTIAL — 5 of 20 lith draws complete as of 2026-10-25T20:10Z.** The nolith comparison
> and draws 1–5 are done; draws 6–20 are in progress (draw 6 running; cap at 22:16 UTC;
> draws 17–20 need next session). Final scores will be computed from all 20 draws.
> Numbers below are provisional and will be revised. Do not quote them.

**Date:** 2026-10-25
**Analysis:** `network/runs/inp-023` (age+sex+rin + ILR b1–b4 residualised, 42-module reference on 234 controls, 20 lith draws)
**Evidence:** PARTIAL (see below) · **Interest:** TBD · **Δ:** TBD

## Design

Residualise all expression data for both technical covariates (age, sex, RIN) AND cell
composition (ILR b1–b4) before building the reference and test networks. This is the most
conservative form of residualisation: any preserved structure must survive both sources of
confounding. `inp-022` tested ILR-only; `inp-023` adds the three technical covariates
used in the de_v2 primary model.

Pre-specified covariates: `~ age + sex + rin + b1 + b2 + b3 + b4`. Plate and seqPCs
omitted to avoid over-fitting (stated in `config/inp-023.json` before computation).

| arm | expression | modules | instrument |
|:--|:--|--:|:--|
| base-001 (raw) | offset-removed, raw | 71 | 50 perms/draw, 20 draws |
| inp-022 (ILR) | offset-removed, ILR residualised | 44 | 50 perms/draw, 20 draws |
| inp-023 (this) | offset-removed, technical+ILR residualised | 42 | 50 perms/draw, 20 draws |

## Reference network

42 modules (deepSplit=4, mergeCutHeight=0, power=12, bicor, signed, n=234 controls after
doubly-residualising). Two fewer than inp-022 (44), consistent with residualising for
three additional technical covariates removing a small amount of additional co-expression
structure.

SFT curve (pre-computed, attempt 1):

| arm | first power R²≥0.80 | max R² | powerEstimate |
|:--|--:|--:|--:|
| base-001 (raw) | never | 0.770 | NA |
| inp-022 (ILR) | power 9 | 0.880 | 16 |
| inp-023 (this) | power 13 | 0.836 | NA |

Adding technical covariates on top of ILR degrades the SFT criterion slightly (0.880 →
0.836, no longer produces a powerEstimate) — consistent with technical covariates
capturing some genuine regulatory variation whose removal reduces the scale-free structure.

## Provisional results (4 of 20 lith draws)

**bp_nolith** (500 perms, n=74 fixed, complete):

| metric | inp-023 | inp-022 (for reference) |
|:--|--:|--:|
| gold Zsummary | 28.08 | 30.46 |
| gold medianRank | 43/43 (last) | 45/45 (last) |
| modules above gold (Zsummary) | 25 of 42 (60%) | 22 of 44 (50%) |

**bp_lith** (50 perms per draw, 5 of 20 draws):

| metric | inp-023 provisional | inp-022 final |
|:--|--:|--:|
| gold mean Z | 15.08 (sd 7.82, 5 draws) | 15.55 (sd 8.13, 20 draws) |
| gold medianRank | **43/43 in all 5 draws** | 45/45 in all 20 draws |
| modules above gold (per-draw mean) | 38.2/42 (91.0%) | 41/44 (93%) |
| draw-by-draw gold Z | 29.0, 11.7, 10.9, 12.4, 11.4 | — |

Gold ranks **last** in every draw so far, matching the inp-022 pattern. The high draw-to-draw
variance in gold Z (sd 7.82, range [10.89, 29.04]) matches inp-022 (sd 8.13). The fraction
of real modules beating gold (~91% at 5 draws) is consistent with inp-022 (93% final).

## Interpretation (provisional)

These partial results are consistent with **inp-022's rescue replicating under a stricter
residualisation.** Removing technical covariates on top of ILR composition produces a
42-module reference that behaves similarly to inp-022's 44-module reference:

- Gold always ranks last → all real modules beat the random gene set by medianRank
- ~89–93% of modules beat gold by Zsummary within each draw
- High draw-to-draw variance in gold Z → heterogeneity of co-expression coherence across
  bp_lith subsamples persists after technical+ILR residualisation

If this pattern holds across all 20 draws, the composition interpretation is robust to:
1. ILR-only residualisation (inp-022) → stable
2. Technical + ILR residualisation (inp-023) → replicates ✓ (provisional)

## Scores

**Will be assigned when all 20 draws complete.** Anticipated:
- E=5 if gold ranks last in all 20 draws and ≥39/42 modules beat gold (mean Z)
- R=5 if inp-022 pattern is replicated consistently
- B=5 (was 4 in inp-022 because inp-023 was pending; now provides the robustness check)
- D=5 if gold is always last by medianRank
- C=4 (same as inp-022; cross-arm comparison remains invalid due to different references)

At E=5 R=5 B=5 D=5 C=4 = 24/25 — the same level as findings #1 and #2 if the pattern
holds. This would make inp-022's Robustness score rise from 4 to 5 (E: 23 → 24).

**Note on inp-022 Robustness update:** per inp-022 finding's own caveat, B=4 because
"a single residualisation formula is tested; inp-023 tests technical + ILR jointly; if
inp-023 replicates, B rises to 5." If the final 20-draw summary replicates, update
inp-022 finding to E=24/25 and update RANKING.md accordingly.
