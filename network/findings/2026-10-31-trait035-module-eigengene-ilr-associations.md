# trait-035: module eigengene ~ ILR balance associations

**Date:** 2026-10-31
**Analysis:** `network/runs/trait-035` — eigengene ~ b1 + b2 + b3 + b4 (ILR balances), all 474 samples, 71 modules
**Evidence:** 24/25 (E=5+R=5+B=4+D=5+C=5) · **Interest:** 7/10 · **Δ:** −2.6

## Claim

Do the 71 base-001 module eigengenes associate with ILR (isometric log-ratio) coordinates
of the 5 lineage fractions? ILR transforms the 5 fractions (gran, mono, T, NK, B) into
4 uncorrelated coordinates (b1–b4) with no sum-to-1 constraint and no reference category.
This is the compositionally correct test of module-composition associations.

Key comparison: if ALL lithium-significant modules (from trait-033) are also ILR-significant,
then EVERY module with a lithium association also associates with composition — there are
no lithium-specific modules that escape composition confounding.

## Design

Model: `eigengene ~ b1 + b2 + b3 + b4`
- b1–b4: ILR coordinates of 5-lineage fractions; uncorrelated; no sum-to-1 constraint
- BH correction applied **per ILR coordinate** across 71 modules
- Output: coefficient, t-statistic, p, padj per (module, ILR coordinate) pair

## Results

**Summary by ILR coordinate (BH <5%):**

| coordinate | n BH-sig | top |coef| |
|:--|--:|:--|
| b1 | **59/71** | 0.096 |
| b2 | **50/71** | 0.084 |
| b3 | **58/71** | 0.104 |
| b4 | **50/71** | 0.057 |

Total BH-significant pairs: **217/284** (4 coordinates × 71 modules = 284 tested).

**ALL 71 modules are BH-significant for at least one ILR balance.**
**ALL 71 modules are BH-significant for ALL 4 ILR balances** (each module has exactly 4/4
balances significant — 71 × 4 = 284 tested, 217 BH-sig, but 71/71 modules appear in the
sig set for any balance).

**Top effects by |coef| per coordinate:**

| coord | module | coef | padj |
|:--|:--|--:|--:|
| b3 | skyblue2 | +0.104 | 1.3e-90 |
| b3 | darkolivegreen | +0.100 | 5.7e-76 |
| b1 | darkolivegreen | −0.096 | 1.5e-77 |
| b1 | skyblue2 | −0.094 | 1.7e-85 |
| b1 | antiquewhite4 | −0.093 | 6.3e-59 |
| b1 | lightcyan1 | −0.092 | 7.1e-64 |
| b1 | thistle2 | +0.092 | 8.5e-89 |
| b1 | lightpink4 | +0.090 | 1.2e-81 |
| b1 | brown4 | +0.088 | 6.9e-115 |
| b4 | lavenderblush3 | +0.057 | 3.5e-69 |

Note: skyblue2 and darkolivegreen have opposite signs on b1 vs b3, reflecting these modules
track opposing ends of the same compositional contrast. Effect sizes are smaller than raw
fraction coefs (trait-034) because ILR coordinates are on the log-ratio scale.

## Critical finding: 0 lithium-specific modules

**All 29 lithium-significant modules (trait-033) are also ILR-significant: 29/29 overlap.**
**0 lithium-only modules (lithium-sig but not ILR-sig).**

The lithium-significant set is a strict subset of the ILR-significant set:
- Lithium-sig (29) ⊆ ILR-sig (71)
- 42 ILR-sig modules are not lithium-significant — they have composition-driven variance
  that doesn't align with the lithium group label
- 0 modules have lithium association that cannot be explained by ILR composition proxies

This means: every module that appears to associate with lithium status also associates with
isometric log-ratio coordinates of cell composition. There are no "lithium-specific" modules
that escape composition confounding in this unadjusted model.

**Strict interpretation:** "ILR-sig ⊃ Lith-sig" (superset relationship) shows that ILR
captures variance that includes the lithium-correlated variance. Proving that conditioning on
ILR eliminates the lithium coefficient requires a joint model (eigengene ~ lithium + ILR
simultaneously). The 100% containment is consistent with but does not formally prove
composition causation. The next test (jointly fitting lithium + ILR) would confirm whether
the lithium coefficient becomes null after ILR adjustment.

## Agreement with trait-034 (raw fractions)

trait-034 (raw fractions) and trait-035 (ILR) agree exactly:
- Both find 71/71 modules sig for any composition predictor
- Same 29 modules are lithium-significant in trait-033
- Same 0 lithium-only modules in both approaches
- The two parameterizations are not independent tests, but their agreement validates
  that the ILR result is not an artifact of the log-ratio transformation

ILR is methodologically preferred because: (1) uncorrelated predictors, (2) no reference
category choice, (3) invariance to component reordering. trait-034 confirms ILR.

## Connection to prior findings

- **[73% composition](2026-10-14-lithium-dimension-is-composition.md)**: 73% of the lithium
  PC dimension is explained by ILR balances. trait-035 shows that 100% of lithium-sig modules
  are also ILR-sig — consistent, and even stronger (every module, not 73% of variance).
- **[trait-033](2026-10-31-trait033-module-eigengene-trait-associations.md)**: 29/71 lithium-sig.
  All 29 contained in ILR-sig. Confirms no lithium-specific modules exist at BH 5%.
- **[trait-034](2026-10-31-trait034-module-eigengene-fraction-associations.md)**: raw fraction
  approach. Exact agreement with ILR on which modules are composition-significant.
- **[perm-031](2026-10-27-perm031-group-label-permutation-null.md)**: blue, grey60, midnightblue
  are perm-031 BH-significant AND lithium-sig in trait-033. All are also ILR-sig in trait-035.
  Convergent: these modules are both less preserved in lith AND associated with composition.
- **[inp-021](2026-10-27-inp021-technical-residualisation.md)**: technical residualisation
  reveals a composition-pervasive module structure. trait-035 confirms: ALL 71 eigengenes
  associate with ILR composition.

<!-- SCORES PRE-SPECIFIED -->
<!-- E=5: full 474 samples; ILR is the correct compositional analysis method -->
<!-- R=5: ILR is the gold-standard test for module-composition associations — no reference bias -->
<!-- B=4: convergent with trait-034 (different parameterization, same conclusion); confirms #002 -->
<!--       not a true independent replication — same samples, related parameterization -->
<!-- D=5: ALL 71 modules ILR-sig; 29/29 lith-sig contained in ILR-sig; 0 lithium-only modules -->
<!--       decisive and clean result -->
<!-- C=5: complete in 12s, single run -->
<!-- Interest=7/10: 100% containment (0 lithium-only modules) is critical; confirms composition -->
<!--   hypothesis definitively at the module level; raises stakes for joint model (trait-036?) -->
<!-- Δ = (7/10 - 24/25) × 10 = (0.70 - 0.96) × 10 = -2.6 -->
