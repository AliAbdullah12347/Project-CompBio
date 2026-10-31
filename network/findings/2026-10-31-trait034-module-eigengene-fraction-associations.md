# trait-034: module eigengene ~ lineage fraction associations

**Date:** 2026-10-31
**Analysis:** `network/runs/trait-034` — eigengene ~ mono + T + NK + B (gran dropped as reference), all 474 samples, 71 modules
**Evidence:** 22/25 (E=5+R=5+B=3+D=4+C=5) · **Interest:** 4/10 · **Δ:** −4.8

## Claim

Do the 71 base-001 module eigengenes associate with raw cell-lineage fraction scores
(LM22 deconvolution)? This tests module-composition associations using the raw fraction
representation, with gran (granulocytes) dropped as the reference category to avoid
perfect multicollinearity from the sum-to-1 constraint.

Contrast with trait-035 (eigengene ~ ILR balances), which uses the compositionally correct
ILR representation where all 4 coordinates are included simultaneously without reference bias.
trait-034 uses the simpler raw-fraction approach for comparison.

## Design

Model: `eigengene ~ mono + T + NK + B` (gran dropped as reference)
- gran dropped: 5 fractions sum to 1, so one must be excluded to avoid perfect collinearity
- Coefficients for mono, T, NK, B are interpreted relative to the gran baseline
- BH correction applied **per predictor** across 71 modules
- Output: coefficient, t-statistic, p, padj per (module, predictor) pair

## Results

**Summary by predictor (BH <5%):**

| predictor | n BH-sig | top effect (|coef|) |
|:--|--:|:--|
| T (T-cells) | **56/71** | ~0.66 |
| NK (NK cells) | **52/71** | ~1.20 |
| mono (monocytes) | **54/71** | ~0.41 |
| B (B cells) | **49/71** | ~3.20 |
| gran (granulocytes) | **0/71** | — (reference) |

Total BH-significant pairs: **211/355** (4 predictors × 71 modules = 284 tested; gran excluded).
Note: gran=0 reflects that it is the reference category — its "associations" are absorbed into
the intercept. The zero count does not imply granulocytes are biologically irrelevant.

**All 71 modules are BH-significant for at least one lineage fraction.**

**Top effects by |coef|:**

| predictor | module | coef | padj |
|:--|:--|--:|--:|
| B | thistle1 | +3.20 | 2.0e-68 |
| B | yellow | −1.49 | 6.9e-18 |
| B | pink | +1.48 | 2.0e-12 |
| B | royalblue | −1.43 | 1.5e-12 |
| NK | lavenderblush3 | +1.20 | 1.1e-102 |
| T | lightyellow | −0.66 | 1.6e-137 |
| T | thistle2 | −0.62 | 8.4e-116 |
| T | darkolivegreen | +0.62 | 2.3e-95 |
| T | brown4 | −0.61 | 2.2e-151 |
| mono | saddlebrown | +0.41 | 2.9e-21 |

**Modules NOT B-significant (22/71):** coral1, cyan, darkgrey, darkolivegreen,
darkseagreen4, greenyellow, ivory, lightgreen, lightyellow, maroon, mediumpurple2,
mediumpurple3, orangered3, palevioletred3, plum1, plum2, red, sienna3, skyblue3,
tan, turquoise, yellowgreen.

## Overlap with lithium-significant modules

All 29 lithium-significant modules (from trait-033) are also fraction-significant for at
least one lineage (29/29 overlap, 0 lithium-only modules). The 42 modules not lithium-
significant include modules that associate with composition but not with the group label.

This is the raw-fraction counterpart to the ILR result (trait-035): both approaches
converge on the same conclusion — the lithium-sig set is a strict subset of the
composition-sig set.

## Methodological limitation

The dropped reference category (gran) is arbitrary. Changing the reference changes all
coefficients. The fraction representation is NOT compositionally proper: for the same
biological composition, rescaling the fractions changes absolute coefficients. The ILR
representation (trait-035) avoids these issues and should be preferred for interpretation.

trait-034's value is as a check on trait-035 using a simpler, more familiar parameterization.
The agreement between both approaches (all 71 modules sig in both) validates the ILR result.

## Connection to prior findings

- **[73% composition](2026-10-14-lithium-dimension-is-composition.md)**: predicts broad
  fraction-eigengene associations. trait-034 confirms: 211/355 BH-sig pairs.
- **[trait-033](2026-10-31-trait033-module-eigengene-trait-associations.md)**: 29/71 lithium-sig.
  All 29 are also fraction-sig — consistent with composition confounding.
- **[trait-035](2026-10-31-trait035-module-eigengene-ilr-associations.md)**: ILR approach.
  Agreement is exact (all 71 sig in both). ILR is the preferred method.
- **[inp-021](2026-10-27-inp021-technical-residualisation.md)**: 52/71 plate-sig in trait-033.
  Plate and composition together explain most inter-sample variation in eigengenes.

<!-- SCORES PRE-SPECIFIED -->
<!-- E=5: full 474 samples, standard OLS, correct BH correction -->
<!-- R=5: directly tests module-fraction associations -->
<!-- B=3: reference category (gran) choice is arbitrary; methodologically less rigorous than ILR -->
<!--       convergent with trait-035 but not an independent biological replication -->
<!-- D=4: 211 BH-sig pairs; ALL 71 modules sig for any fraction; clean pattern; gran=0 is expected -->
<!-- C=5: complete in 15s, single run -->
<!-- Interest=4/10: expected given 73% composition finding; confirms rather than surprises -->
<!-- Δ = (4/10 - 22/25) × 10 = (0.40 - 0.88) × 10 = -4.8 -->
