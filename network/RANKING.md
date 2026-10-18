# Arm 3 — Live Result Ranking

> **GENERATED FILE — do not hand-edit.** Rebuild with `Rscript network/build_index.R`.
> Edit the per-finding files in `network/findings/` instead. On a merge conflict here,
> regenerate rather than merge. See `CONTRIBUTING.md`.

**Findings:** 21   **Updated:** 2026-10-18

Two rankings. **Table A (Evidence, 0–25)** governs the write-up; it rewards evidential
strength so a well-powered null can score top. **Table B (Interest, 0–10)** records how
striking a result is. **Δ = (Interest/10 − Evidence/25) × 10**: positive means the story
is ahead of the proof.

## Table A — Evidence

| # | Evidence | Interest | Δ | finding | components | source |
|--:|--:|--:|--:|:--|:--|:--|
| 1 | **24/25** | 10/10 | +0.4 | [The reference network is built on an un-removed per-sample library offset](findings/2026-10-15-library-offset-artifact.md) | E=5 R=5 B=5 D=5 C=4 | network/analyses/015-referee-verification |
| 2 | **24/25** | 9/10 | -0.6 | [A lithium-associated co-expression dimension exists — and it is 73% cell composition](findings/2026-10-14-lithium-dimension-is-composition.md) | E=5 R=5 B=5 D=5 C=4 | network/analyses/002-effective-dimensionality |
| 3 | **23/25** | 10/10 | +0.8 | [The corrected network has 71 modules, and its eigengenes are near-orthogonal](findings/2026-10-18-rebuilt-network-71-modules.md) | E=5 R=5 B=5 D=5 C=3 | network/runs/mod-000-resolution (rebuilt on offset-removed expression) |
| 4 | **22/25** | 9/10 | +0.2 | [Darkturquoise is a ribosomal/translation module: the composition signal has a biological substrate](findings/2026-10-14-darkturquoise-is-ribosomal.md) | E=5 R=5 B=2 D=5 C=5 | network/analyses/006-darkturquoise-annotation |
| 5 | **22/25** | 9/10 | +0.2 | [Removing the library offset rebuilds the network into a well-behaved one](findings/2026-10-18-normalisation-rebuild-required.md) | E=5 R=5 B=4 D=5 C=3 | network/analyses/016-normalisation-decision |
| 6 | **22/25** | 7/10 | -1.8 | [The lithium dimension replicates across decompositions, and never survives composition](findings/2026-10-14-lithium-dimension-replicates.md) | E=5 R=4 B=4 D=5 C=4 | network/analyses/004-decomposition-robustness |
| 7 | **21/25** | 9/10 | +0.6 | [Real modules do beat the random gene set — once the decomposition is real](findings/2026-10-18-modules-do-beat-gold.md) | E=5 R=4 B=4 D=5 C=3 | network/runs/base-001 (offset-removed, 71-module reference) |
| 8 | **20/25** | 8/10 | +0.0 | [Darkred is an erythrocyte module: stable because RBC genes are absent from immune composition](findings/2026-10-14-darkred-is-erythrocyte-module.md) | E=4 R=0 B=2 D=4 C=5 | network/analyses/011-darkred-annotation |
| 9 | **20/25** | 5/10 | -3.0 | [The lithium-composition result is not circular: it survives removing every marker gene](findings/2026-10-14-lm22-circularity-passes.md) | E=5 R=4 B=2 D=5 C=4 | network/analyses/003-lm22-circularity |
| 10 | **19/25** | 7/10 | -0.6 | [PC2 is a lithium axis in the full cohort too, and the lithium DE signal is module-specific](findings/2026-10-14-pc2-lithium-axis-replicates-full-cohort.md) | E=4 R=4 B=3 D=4 C=4 | network/analyses/011-module-annotation |
| 11 | **18/25** | 7/10 | -0.2 | [Darkturquoise drives PC2: the lithium-associated dimension is localised to one module](findings/2026-10-14-darkturquoise-drives-pc2.md) | E=4 R=4 B=2 D=4 C=4 | network/analyses/005-pc2-module-loadings |
| 12 | **18/25** | 7/10 | -0.2 | [Darkturquoise group×composition interaction: the control-lith contrast survives proper slope correction](findings/2026-10-14-darkturquoise-interaction.md) | E=4 R=4 B=2 D=4 C=4 | network/analyses/008-darkturquoise-interaction |
| 13 | **18/25** | 6/10 | -1.2 | [No module eigengene associates with lithium — and the modules are not independent](findings/2026-10-14-module-trait-no-lithium-association.md) | E=4 R=5 B=2 D=4 C=3 | network/analyses/001-module-trait-ds4 |
| 14 | **17/25** | 8/10 | +1.2 | [Darkturquoise: control-lith contrast survives composition adjustment; nolith-lith does not](findings/2026-10-14-darkturquoise-ctrl-lith-survives-adjustment.md) | E=4 R=4 B=2 D=3 C=4 | network/analyses/007-darkturquoise-across-groups |
| 15 | **16/25** | 8/10 | +1.6 | [BP disease tightens inter-module correlations: modules move more in lockstep](findings/2026-10-14-inter-module-correlations-tighten-in-bp.md) | E=3 R=4 B=2 D=0 C=4 | network/analyses/013-inter-module-correlations |
| 16 | **6/25** | 4/10 | +1.6 | [No module reliably beats a random gene set (degenerate decomposition)](findings/2026-10-09-gold-baseline-no-module-beats-random.md) | E=5 R=5 B=3 D=5 C=2 | network/runs/base-001 |
| 17 | **6/25** | 1/10 | -1.4 | [Cost scaling, measured](findings/2026-10-09-cost-scaling.md) | E=4 R=0 B=1 D=1 C=0 | Phase 0 probes + base-001 draw timings |
| 18 | **0/25** | 0/10 | +0.0 | [The merge threshold is decisive: 33 modules vs 3](findings/2026-10-09-merge-threshold-decisive.md) | E=5 R=5 B=3 D=5 C=4 | network/runs/mod-000-resolution |
| 19 | **0/25** | 0/10 | +0.0 | [Scale-free topology fails at every power and both sample sizes](findings/2026-10-09-scale-free-fails.md) | E=4 R=2 B=3 D=4 C=2 | network/runs/base-001 (sft curve) + Phase 0 probes |
| 20 | **0/25** | 0/10 | +0.0 | [BP disease drives network-wide module differences; lithium adds nothing beyond composition](findings/2026-10-14-bp-disease-drives-network-globally.md) | E=4 R=5 B=2 D=4 C=4 | network/analyses/009-all-modules-across-groups |
| 21 | **0/25** | 0/10 | +0.0 | [PC1 of the full-sample ME matrix captures BP disease: a global network factor](findings/2026-10-14-pc1-captures-bp-disease.md) | E=4 R=4 B=2 D=4 C=4 | network/analyses/010-full-pca-all-samples |

## Table B — Interest

| # | Interest | Evidence | Δ | finding |
|--:|--:|--:|--:|:--|
| 1 | **10/10** | 24/25 | +0.4 | [The reference network is built on an un-removed per-sample library offset](findings/2026-10-15-library-offset-artifact.md) |
| 2 | **10/10** | 23/25 | +0.8 | [The corrected network has 71 modules, and its eigengenes are near-orthogonal](findings/2026-10-18-rebuilt-network-71-modules.md) |
| 3 | **9/10** | 24/25 | -0.6 | [A lithium-associated co-expression dimension exists — and it is 73% cell composition](findings/2026-10-14-lithium-dimension-is-composition.md) |
| 4 | **9/10** | 22/25 | +0.2 | [Darkturquoise is a ribosomal/translation module: the composition signal has a biological substrate](findings/2026-10-14-darkturquoise-is-ribosomal.md) |
| 5 | **9/10** | 22/25 | +0.2 | [Removing the library offset rebuilds the network into a well-behaved one](findings/2026-10-18-normalisation-rebuild-required.md) |
| 6 | **9/10** | 21/25 | +0.6 | [Real modules do beat the random gene set — once the decomposition is real](findings/2026-10-18-modules-do-beat-gold.md) |
| 7 | **8/10** | 20/25 | +0.0 | [Darkred is an erythrocyte module: stable because RBC genes are absent from immune composition](findings/2026-10-14-darkred-is-erythrocyte-module.md) |
| 8 | **8/10** | 17/25 | +1.2 | [Darkturquoise: control-lith contrast survives composition adjustment; nolith-lith does not](findings/2026-10-14-darkturquoise-ctrl-lith-survives-adjustment.md) |
| 9 | **8/10** | 16/25 | +1.6 | [BP disease tightens inter-module correlations: modules move more in lockstep](findings/2026-10-14-inter-module-correlations-tighten-in-bp.md) |
| 10 | **7/10** | 22/25 | -1.8 | [The lithium dimension replicates across decompositions, and never survives composition](findings/2026-10-14-lithium-dimension-replicates.md) |
| 11 | **7/10** | 19/25 | -0.6 | [PC2 is a lithium axis in the full cohort too, and the lithium DE signal is module-specific](findings/2026-10-14-pc2-lithium-axis-replicates-full-cohort.md) |
| 12 | **7/10** | 18/25 | -0.2 | [Darkturquoise drives PC2: the lithium-associated dimension is localised to one module](findings/2026-10-14-darkturquoise-drives-pc2.md) |
| 13 | **7/10** | 18/25 | -0.2 | [Darkturquoise group×composition interaction: the control-lith contrast survives proper slope correction](findings/2026-10-14-darkturquoise-interaction.md) |
| 14 | **6/10** | 18/25 | -1.2 | [No module eigengene associates with lithium — and the modules are not independent](findings/2026-10-14-module-trait-no-lithium-association.md) |
| 15 | **5/10** | 20/25 | -3.0 | [The lithium-composition result is not circular: it survives removing every marker gene](findings/2026-10-14-lm22-circularity-passes.md) |
| 16 | **4/10** | 6/25 | +1.6 | [No module reliably beats a random gene set (degenerate decomposition)](findings/2026-10-09-gold-baseline-no-module-beats-random.md) |
| 17 | **1/10** | 6/25 | -1.4 | [Cost scaling, measured](findings/2026-10-09-cost-scaling.md) |
| 18 | **0/10** | 0/25 | +0.0 | [The merge threshold is decisive: 33 modules vs 3](findings/2026-10-09-merge-threshold-decisive.md) |
| 19 | **0/10** | 0/25 | +0.0 | [Scale-free topology fails at every power and both sample sizes](findings/2026-10-09-scale-free-fails.md) |
| 20 | **0/10** | 0/25 | +0.0 | [BP disease drives network-wide module differences; lithium adds nothing beyond composition](findings/2026-10-14-bp-disease-drives-network-globally.md) |
| 21 | **0/10** | 0/25 | +0.0 | [PC1 of the full-sample ME matrix captures BP disease: a global network factor](findings/2026-10-14-pc1-captures-bp-disease.md) |
