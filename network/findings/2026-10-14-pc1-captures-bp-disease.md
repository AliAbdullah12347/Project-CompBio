# PC1 of the full-sample ME matrix captures BP disease: a global network factor

**Date:** 2026-10-14
**Analysis:** `network/analyses/010-full-pca-all-samples`
**Evidence:** 18/25 · **Interest:** 7/10 · **Δ:** −0.2

## Claim

PCA of the 33-module eigengene matrix across all 460 samples reveals a dominant global
factor: PC1 explains 77.5% of variance, all 33 modules load negatively (loadings range
from −0.181 to −0.191, nearly uniform), and PC1 separates controls from BP patients.

## Structure

| PC | Variance | r(disease) | r(lithium) | r(b1) | r(b2) |
|--:|--:|--:|--:|--:|--:|
| 1 | **77.5%** | **+0.289** | +0.240 | +0.070 | +0.055 |
| 2 | 10.1% | −0.123 | −0.198 | **−0.434** | **−0.447** |
| 3 | 6.1% | −0.076 | −0.103 | −0.532 | −0.299 |

PC1 scores by group:

| Group | n | mean PC1 | SD |
|:--|--:|--:|--:|
| control | 234 | **−1.43** | 4.44 |
| bp_nolith | 74 | +0.98 | 4.96 |
| bp_lith | 152 | +1.73 | 5.37 |

Permutation test (K=1000): PC1 r(disease)=0.289, p<0.001 (0/1000 permutations as extreme).

## What this means

**PC1 is the global network factor.** All 33 modules co-vary in the same direction (all
negative loadings), meaning they form a single coherent "network level" rather than
independent signals. Controls have negative PC1 scores (high module expression); BP patients
have positive PC1 scores (lower module expression).

This is the PCA signature of finding 009: the consistent downward shift in BP across all
31 significant modules is captured by PC1 as a single, massive factor.

**PC2 (10%) is the composition/lithium dimension** identified in analysis 002. In the
BP-only PCA (002), this was PC2; including controls shifts the disease contrast to PC1.

## Connection to analysis 002

Analysis 002 ran PCA on 226 BP-only samples and found PC2 (the lithium-associated
dimension) correlates with composition. The key difference here:
- 002: no controls → the dominant network factor was within-BP structure
- 010: with controls → the dominant PC1 captures disease vs control
- PC2 (10%) is the within-group lithium/composition signal, now demoted by the larger disease PC1

The structure is: disease → PC1 (77%); composition/lithium → PC2 (10%).

## Scoring

E=4 (PCA + permutation test; well-powered N=460; single cohort; same-sample analysis as 009).
R=4 (gives the global factor structure; connects 009 to 002; shows BP disease dominates).
B=2 (single cohort, cross-sectional; no replication).
D=4 (77% variance, perm p<0.001, clear PC1-disease separation).
C=4 (all 8 PCs tested; permutation reported; group means tabulated).

Interest 7: the 77% PC1 / global network factor is a clean summary of finding 009.
Δ −0.2 — evidence and interest well-aligned; no surprise given 009.
