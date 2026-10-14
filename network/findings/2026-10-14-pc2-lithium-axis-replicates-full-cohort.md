# PC2 is a lithium axis in the full cohort too, and the lithium DE signal is module-specific

**Date:** 2026-10-14
**Analysis:** `network/analyses/011-module-annotation`
**Evidence:** 19/25 · **Interest:** 7/10 · **Δ:** -0.6

## Claim
Eigengenes over all 474 labelled samples, ANOVA of each PC on group:

| PC | var | p (BH) | control | bp_nolith | bp_lith |
|:--|--:|--:|--:|--:|--:|
| PC1 | 77.5% | 5.2e-09 | −1.43 | +0.98 | +1.73 |
| PC2 | 10.1% | 1.8e-04 | +0.22 | +0.36 | **−0.52** |

**PC1 is a diagnosis axis; PC2 is a lithium axis** — control and bp_nolith sit together,
bp_lith separates. This is the same component analyses 002/004 found lithium-associated
**within BP1**, now recovered on the full three-group cohort. Different sample set, same
axis, so PC2 is not an artifact of the within-BP1 contrast.

**Lithium DE signal is module-specific.** 6 of 33 modules are enriched for cell-type
lithium DEGs (saddlebrown OR 11.05, greenyellow 3.99), and one is significantly
*depleted* — green: 0 observed against 8.1 expected over 835 genes, BH 4.0e-3.

## Secondary
Every module's top hub has kME 0.95–0.98, a suspiciously uniform range. At mean
connectivity of 8.3% of the genome every module has a near-perfect eigengene proxy, so hub
identity carries less information on this data than it normally would and should not be
over-interpreted.

## Scoring
E=4 (Fisher+BH, ANOVA, deterministic; descriptive, no permutation null).
R=4 (independently corroborates the headline axis on a different sample set, and links the
network arm to the DE arm). B=3 (one decomposition, but a genuinely different cohort than
002/004 used). D=4. C=4.

Interest 7: the clean PC1=diagnosis / PC2=lithium split, and a module with *zero* lithium
DEGs where 8 were expected.

## Limits
Descriptive; group separation on a PC is not a calibrated test of network difference, which
is what `base-003-remerge` and `perm-031` provide. Depletion in green could reflect gene
length, expression level or detection-floor effects not modelled here.
