# samp-043: eigengene-space sample structure by group

**Date:** 2026-11-01
**Analysis:** `network/runs/samp-043` — 71 base-001 module eigengenes on all 474 samples; PCA, hclust k=3, centroid-distance outliers, per-module ANOVA by group
**Evidence:** 15/25 (E=3 R=4 B=2 D=2 C=4) · **Interest:** 4/10 · **Δ:** −2.0

## Claim
Does any cluster of samples in eigengene space track lithium?

## Results (descriptive; no composition adjustment)
- PC1 26.9%, PC2 24.5% of eigengene variance. Group centroids on PC1/PC2 lie within about 0.05 of
  each other against within-group sd ≈ 0.2 (see `group_centroids_pc12.csv`): no visible group separation.
- hclust k=3 vs group (bp_lith / bp_nolith / control; 14 samples labelled "other" excluded from the test):
  cluster 1 = 77/24/91, cluster 2 = 59/43/101, cluster 3 = 16/7/42. Chi-square on the 3×3 table:
  X²=14.14, df=4, p=0.0069 (nominal, single test, not corrected). Lith fraction: 39% / 28% / 24%
  in clusters 1/2/3. The association is real but weak, and clusters are far from pure.
- Outliers (centroid distance): 14 samples, 8 control, 1 bp_nolith, 5 bp_lith. This is proportionate to group sizes, so no group-specific outlier excess.
- 44/71 module eigengenes differ by group (ANOVA, BH 5%, unadjusted).

## Interpretation and limits
The weak lithium-enriched cluster is expected from trait-033/035: every lithium-associated eigengene
is composition-associated (0 lithium-specific modules). I did not test whether cluster membership
survives ILR/technical adjustment, so I cannot say it is anything beyond composition. The ANOVA counts overlap with
trait-033 and are not independent evidence. The eigengenes come from a reference built on controls
only, so structure in bipolar groups is measured with control-defined modules.
No p-value here comes from a label-permutation null.

<!-- SCORES: E=3 (descriptive, unadjusted, no permutation null) R=4 B=2 (restates trait-033/035) D=2 C=4 -->
