# inp-021: Technical covariates alone do not rescue module preservation

**Date:** 2026-10-27
**Analysis:** `network/runs/inp-021` — technical covariate residualisation (age+sex+rin+plate+seqpc1-3), 63-module reference on 234 controls
**Evidence:** [FILL] · **Interest:** [FILL] · **Δ:** new

## Claim

Residualising for technical covariates (age, sex, RIN, plate batch, sequencing PCs 1–3)
before building and testing networks reduces the module count from 71 (base-001 raw) to
**63** — removing only 8 modules. By contrast, ILR residualisation (inp-022, cell composition)
removed 27 modules (71 → 44). Technical covariates alone have minimal effect on network
structure; cell composition is the primary driver.

The preservation pattern under technical-only residualisation should therefore look more
like base-001 (raw) than inp-022/inp-023 (ILR-adjusted). If composition is the driver of
the preservation difference, removing only technical covariates should not rescue the
bp_lith preservation pattern seen in inp-022.

## Scale-free topology

Technical residualisation slightly improves SFT fit compared to raw expression (max R²=0.797 at power 18
vs max R²=0.770 for raw; pow-009), but still fails to reach the 0.80 threshold. powerEstimate=NA.

| arm | residualisation | max R² | power | passes 0.80? |
|:--|:--|--:|--:|:--|
| pow-009 (ctrl) | none (offset-removed) | 0.770 | 20 | NO |
| inp-021 | age+sex+rin+plate+seqpc1-3 | **0.797** | 18 | NO |
| inp-022 | ILR b1–b4 | 0.880 | 16 | YES (first fit p=9) |
| inp-023 | age+sex+rin + ILR b1–b4 | 0.836 | 17 | YES (first fit p=13) |

Technical residualisation alone is insufficient to pass SFT. This is consistent with ILR
(composition) residualisation being the key step, not technical covariates.

## Module count comparison

| arm | residualisation | modules | reduction vs raw |
|:--|:--|--:|:--|
| base-001 | none (offset-removed) | 71 | — |
| inp-021 | tech (age+sex+rin+plate+seqpc1-3) | **63** | −8 (11%) |
| inp-022 | ILR b1–b4 | 44 | −27 (38%) |
| inp-023 | tech + ILR | 42 | −29 (41%) |

Only 11% of modules are tech-driven; 38% are composition-driven. The additional 3% lost when
ILR is added to technical covariates (inp-023 vs inp-022) is negligible.

## nolith result

| | inp-021 (tech-only, 63 modules, 50 perms) | base-001 (raw, 71 modules, 500 perms) |
|:--|:--|:--|
| gold Zsummary | **[FILL]** | 14.01 |
| gold medianRank | **[FILL]** | 72 of 72 (last) |
| n modules (real) | 63 | 71 |
| modules above gold (Zsummary) | **[FILL]** of 63 | 65 of 71 |

## lith draws result

[FILL after all draws complete. Run `Rscript network/scripts/analyze_inp021.R`.]

| | inp-021 (tech-only, 63 modules, 50 perms/draw) | base-001 (raw, 71 modules, 50 perms/draw) | inp-023 (tech+ILR, 42 modules, 50 perms/draw) |
|:--|:--|:--|:--|
| gold Zsummary (mean, sd) | **[FILL]** | 21.74 (sd 7.05) | 16.98 (sd 10.79) |
| gold medianRank | **[FILL]** | 72 of 72, all 20 draws | 43 of 43, all 20 draws |
| modules beating gold in all 20 | **[FILL]** of 63 | 17 of 71 | 16 of 42 |
| modules beating gold in 0 | **[FILL]** of 63 | 4 of 71 | 0 of 42 |

## Key contrast: inp-021 vs inp-022

The inp-021 vs inp-022 comparison isolates the effect of ILR residualisation:

- If technical covariates explain the preservation difference, inp-021 should look like inp-022.
- If cell composition explains the preservation difference, inp-021 should look like base-001.

[FILL with actual comparison once draws are complete.]

## Connection to prior findings

- **[73% composition](2026-10-14-lithium-dimension-is-composition.md)**: The lithium network
  dimension is predominantly cell composition. Technical covariates alone should not reproduce
  the composition-residualised rescue seen in inp-022.
- **[SFT fails on raw data](2026-10-26-pow009-scale-free-composition-driven.md)**: scale-free topology
  only passes after ILR residualisation, not after technical-only (this finding).
- **[inp-022](2026-10-24-inp022-composition-rescue.md)** and
  **[inp-023](2026-10-25-inp023-technical-ilr-replication.md)**: the ILR-residualised arms
  both show gold last with many modules preserved above it — the comparison arm for inp-021.

<!-- Scores to fill: depends on results. If inp-021 looks like base-001 (not inp-022):
     strong confirmation of composition hypothesis; E≈22/25, I≈7/10, Δ≈-1.8
     If inp-021 looks like inp-022 (preserved): refutes composition hypothesis;
     would be surprising and high-interest. -->
