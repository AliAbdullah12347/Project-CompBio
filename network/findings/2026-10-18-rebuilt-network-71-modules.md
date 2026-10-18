# The corrected network has 71 modules, and its eigengenes are near-orthogonal

**Date:** 2026-10-18
**Analysis:** `network/runs/mod-000-resolution` (rebuilt on offset-removed expression)
**Evidence:** 23/25 · **Interest:** 10/10 · **Δ:** new

## Claim

Rebuilding the reference network on offset-removed expression, with **identical**
parameters (power 12, signed, bicor, deepSplit 4, minModuleSize 30), changes the answer
to the most basic question the arm asks — how many modules are there.

| at power 12, deepSplit 4 | raw `log2(TPM+1)` | offset removed |
|:--|--:|--:|
| **modules (mergeCutHeight 0)** | **33** | **71** |
| largest module | 11.68% | 9.26% |
| mean connectivity k | **1030.84** | **80.96** |
| eigengene correlation, 25th pct | +0.605 | **−0.272** |
| eigengene correlation, median | **+0.765** | **+0.011** |
| eigengene correlation, 75th pct | +0.872 | +0.328 |

## The three things that changed

**Connectivity fell 12.7-fold.** In the raw network the average gene was connected to
1,031 of 12,368 genes — 8.3% of the genome, after a power-12 soft threshold whose entire
purpose is to suppress weak edges. Offset-removed: 81 genes, 0.65%. The raw network was
dense in the way only an artifact can make it.

**The eigengenes went from collinear to orthogonal.** Raw: the 25th percentile of
eigengene correlation was **+0.605** — meaning three quarters of all module pairs
correlated above 0.6, and *no* pair was meaningfully anticorrelated. Offset-removed: the
distribution is symmetric about zero (−0.272 / +0.011 / +0.328). Modules in the corrected
network are nearly independent, which is what modules are supposed to be.

**The merge cascade disappeared.** This is the direct retraction of finding #3.

| mergeCutHeight | raw: modules (largest) | offset removed: modules (largest) |
|--:|:--|:--|
| 0 | 33 (11.7%) | 71 (9.3%) |
| 0.05 | 14 (18.7%) | 71 (9.3%) |
| 0.10 | 8 (28.3%) | 61 (9.6%) |
| 0.15 | 4 (63.9%) | 45 (12.4%) |
| 0.25 | **3 (97.5%)** | 30 (13.0%) |

On raw data the WGCNA default of 0.25 collapsed the network into a single module holding
**97.5% of genes**. That was never a property of blood co-expression; it is what happens
when every eigengene pair already correlates at ~0.77, so the merge step chains through
all of them. On corrected data the same sweep behaves ordinarily — 71 down to 30, largest
module still only 13%.

## What this invalidates

**Finding #3 ("the merge threshold is decisive: 33 modules vs 3") — retract.** It
measured the offset, not the data. mergeCutHeight sensitivity is real and known, but the
catastrophic form reported here does not occur in the corrected network.

**Finding #16 ("scale-free topology fails at every power") — retract.** Analysis 016
already showed the fit slope at power 12 going −0.603 → −1.19 and the R² curve ceasing to
be non-monotone.

**Analysis 002's structural premise is gone.** "33 modules span ~6 effective dimensions,
median inter-eigengene |r| = 0.822" described a network whose eigengenes all correlated
because of one per-sample scalar. The corrected network has 71 modules with median
correlation 0.011.

**Every module-specific finding now refers to modules that may not exist.** darkturquoise,
darkred, saddlebrown and the rest were defined by the 33-module raw decomposition. The
71-module corrected decomposition is a different partition; those labels do not carry
over, and the colour names are assigned by size rank, so same-coloured modules in the two
runs are **not** the same gene set. Any surviving claim has to be re-derived, not
re-labelled.

## What this does not touch

The headline (002's PC2-lithium-composition result) was tested against the library offset
directly in `analyses/015-referee-verification` and held (β −0.6450, p = 0.00748 with the
offset as a covariate). But it was computed on the 33-module eigengenes, so it too needs
re-deriving on the 71-module decomposition before it can be quoted as a corrected result.

E=5 R=5 B=5 D=5 C=3
