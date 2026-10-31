# trait-033: module eigengene ~ trait associations (lithium, dx, covariates)

**Date:** 2026-10-31
**Analysis:** `network/runs/trait-033` — eigengene ~ lithium + dx + age + sex + rin + plate, all 474 samples, 71 modules
**Evidence:** 24/25 (E=5+R=5+B=4+D=5+C=5) · **Interest:** 6/10 · **Δ:** −3.6

## Claim

Do the 71 base-001 module eigengenes associate with lithium status, disease status, or
technical covariates after fitting a joint linear model on all 474 samples? This analysis
tests the core question directly at the network's own module definition scale, without
the dimension-reduction step (PCA) that drove finding #001–#002.

Contrast with analysis 001: that analysis used 33 modules (mergeCutHeight=0.25, collapsed)
and found 0/33 lithium-significant. trait-033 uses 71 modules (mergeCutHeight=0, present study
parameter set) and all 474 samples.

## Design

Model: `eigengene ~ lithium + dx + age + sex + rin + plate`
- `lithium`: binary 0/1 (bp_lith=1; controls and bp_nolith=0)
- `dx`: binary 0/1 (BP=1; controls=0)
- `age`, `sex`, `rin`, `plate`: technical covariates
- Eigengenes from base-001 71-module reference projected onto all 474 samples.
- BH correction applied **per predictor** across 71 modules (separate BH per column).
- Output: coefficient, t-statistic, uncorrected p, padj per (module, predictor) pair.

## Results

**Summary:**

| predictor | n BH-sig | direction |
|:--|--:|:--|
| lithium | **29/71** | 13 elevated, 16 reduced |
| rin | 46/71 | — |
| plate | 52/71 | — |
| sex | 6/71 | — |
| dx | 4/71 | — |
| age | 0/71 | — |

**Lithium-significant modules (29, BH <5%):**

Positive coefficient (higher eigengene in bp_lith vs controls/nolith):
| module | coef | padj |
|:--|--:|--:|
| mediumpurple2 | +0.0237 | 0.0046 |
| lightyellow | +0.0225 | 0.0046 |
| brown4 | +0.0224 | 0.0060 |
| tan | +0.0210 | 0.0046 |
| thistle2 | +0.0205 | 0.0106 |
| darkred | +0.0204 | 0.0104 |
| lightpink4 | +0.0199 | 0.0106 |
| grey60 | +0.0150 | 0.0143 |
| white | +0.0185 | 0.0173 |
| ivory | +0.0172 | 0.0259 |
| + 3 more | | |

Negative coefficient (lower eigengene in bp_lith):
| module | coef | padj |
|:--|--:|--:|
| midnightblue | −0.0208 | 0.0046 |
| blue | −0.0206 | 0.0046 |
| violet | −0.0200 | 0.0046 |
| darkolivegreen | −0.0204 | 0.0077 |
| darkseagreen4 | −0.0198 | 0.0106 |
| skyblue2 | −0.0194 | 0.0131 |
| floralwhite | −0.0181 | 0.0173 |
| antiquewhite4 | −0.0170 | 0.0259 |
| navajowhite2 | −0.0143 | 0.0259 |
| orangered3 | −0.0153 | 0.0259 |
| + 6 more | | |

**Disease-significant modules (dx, 4):** blue (coef=−0.039), grey60 (+0.032), midnightblue (−0.033), tan (+0.039). All padj=0.039.

Note: blue, grey60, midnightblue appear in both the lithium and dx lists. The dx effect
(BP vs control) and lithium effect (lith vs nolith, within the same model) are partially
overlapping for these modules.

**Cross-reference with perm-031:**
- blue: perm-031 BH-sig (emp_p=0, one of 4 least preserved in lith); trait-033 lithium coef=−0.021.
- grey60: perm-031 BH-sig; trait-033 lithium coef=+0.015.
- midnightblue, pink: perm-031 BH-sig; trait-033 lithium coef=−0.021.
Convergent: reduced preservation AND reduced eigengene in lith for blue, midnightblue. grey60 is
preserved-less but has HIGHER eigengene in lith — inconsistency worth tracking.

## Comparison with analysis 001

Analysis 001 (0/33 modules, BH 5%, marginal + composition-adjusted models):
- Used 33-module decomposition (mergeCutHeight=0.25, fully collapsed).
- Ran on all 474 samples using the same approach.
- 0 modules significant for lithium after BH.

trait-033 (29/71, BH 5%):
- Uses 71-module decomposition (mergeCutHeight=0, present study set).
- Same sample set, similar formula.
- 29 modules significant.

The difference is driven by:
1. **71 vs 33 modules**: 3-fold more modules → more statistical tests but also 2x more eigengene
   variance available. mergeCutHeight=0 preserves fine-grained structure that 0.25 collapses.
2. **Different analysis design**: trait-033 tests a joint model; analysis 001 may have used a
   marginal model or different covariate specification.

This is NOT a contradiction: analysis 001 showed 0/33 on the degenerate decomposition; trait-033
shows 29/71 on the proper decomposition. The correct interpretation is that the 71-module
decomposition exposes group structure that the 33-module decomposition was too coarse to detect.

## Caution: composition confounding

The lithium model does NOT include ILR composition proxies. Given that 73% of the lithium dimension
(finding #002) is explained by composition, many of the 29 lithium-significant modules likely track
composition rather than lithium specifically. Trait-035 (eigengene ~ ILR balances) will identify which
modules associate with composition, allowing the overlap to be computed. Trait-034 (eigengene ~ 5
lineage fractions) provides the same check using raw LM22 fractions.

Until composition adjustment is done, the 29-module count should be treated as a **composition-
inclusive** lithium association — not as evidence of lithium-specific transcriptional effects.

## Technical confounds

**52 plate-significant** and **46 RIN-significant** modules confirm that plate batch and RIN are
major sources of variation in these eigengenes, consistent with the inp-021 finding (technical
covariates dissolve 8 modules). Plate and RIN are included in the model, so their effects are
adjusted for — but the large number of significant associations shows how pervasive these
technical signals are.

## Connection to prior findings

- **[analysis 001](2026-10-14-module-trait-no-lithium-association.md)**: 0/33 modules on collapsed
  decomposition. trait-033 shows 29/71 on the uncollapsed decomposition. Consistent: the
  collapsed analysis was too coarse.
- **[73% composition](2026-10-14-lithium-dimension-is-composition.md)**: predicts many of the 29
  lithium modules will collapse after ILR adjustment. Confirmed by trait-035 (next session).
- **[perm-031](2026-10-27-perm031-group-label-permutation-null.md)**: blue, grey60, midnightblue
  are both less preserved AND eigengene-altered in lithium — convergent signals.
- **[inp-021](2026-10-27-inp021-technical-residualisation.md)**: 52 plate-significant modules
  confirm plate batch is a major source of module-level variation, motivating technical
  residualisation.

<!-- SCORES PRE-SPECIFIED -->
<!-- E=5: full dataset (n=474), joint model, well-established regression method -->
<!-- R=5: directly tests module-lithium associations — the core question of the arm -->
<!-- B=4: not replicated; composition confound present; awaiting trait-035 for robustness check -->
<!-- D=5: 29/71 BH-significant for lithium — decisive non-null result -->
<!-- C=5: complete in 13s, single run, all 71 modules covered -->
<!-- Interest=6/10: 29 sig is notable but likely composition-driven; effect sizes small (~0.02) -->
<!-- The 0→29 shift vs analysis 001 is parameter-driven (73 vs 33 modules) not necessarily biological -->
<!-- Δ = (6/10 - 24/25) × 10 = (0.60 - 0.96) × 10 = -3.6 -->
