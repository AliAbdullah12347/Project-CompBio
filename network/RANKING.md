# Arm 3 — Live Result Ranking

> **GENERATED FILE — do not hand-edit.** Rebuild with `Rscript network/build_index.R`.
> Edit the per-finding files in `network/findings/` instead. On a merge conflict here,
> regenerate rather than merge. See `CONTRIBUTING.md`.

**Findings:** 25   **Updated:** 2026-10-25

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
| 4 | **23/25** | 9/10 | -0.2 | [inp-022: ILR residualization rescues bp_lith co-expression preservation](findings/2026-10-24-inp022-composition-rescue.md) | E=5 R=5 B=4 D=5 C=4 | network/runs/inp-022 (ILR b1–b4 residualized, 45-module reference on 234 controls, 20 lith draws) |
| 5 | **22/25** | 9/10 | +0.2 | [Darkturquoise is a ribosomal/translation module: the composition signal has a biological substrate](findings/2026-10-14-darkturquoise-is-ribosomal.md) | E=5 R=5 B=2 D=5 C=5 | network/analyses/006-darkturquoise-annotation |
| 6 | **22/25** | 9/10 | +0.2 | [Removing the library offset rebuilds the network into a well-behaved one](findings/2026-10-18-normalisation-rebuild-required.md) | E=5 R=5 B=4 D=5 C=3 | network/analyses/016-normalisation-decision |
| 7 | **22/25** | 8/10 | -0.8 | [base-001 complete: modules beat gold under subsampling — but Zsummary is not comparable across arms](findings/2026-10-22-base001-complete-and-zsummary-not-comparable.md) | E=5 R=5 B=4 D=5 C=3 | network/runs/base-001 (offset-removed, 71-module reference, 20 draws) |
| 8 | **22/25** | 7/10 | -1.8 | [The lithium dimension replicates across decompositions, and never survives composition](findings/2026-10-14-lithium-dimension-replicates.md) | E=5 R=4 B=4 D=5 C=4 | network/analyses/004-decomposition-robustness |
| 9 | **21/25** | 9/10 | +0.6 | [Real modules do beat the random gene set — once the decomposition is real](findings/2026-10-18-modules-do-beat-gold.md) | E=5 R=4 B=4 D=5 C=3 | network/runs/base-001 (offset-removed, 71-module reference) |
| 10 | **20/25** | 8/10 | +0.0 | [Darkred is an erythrocyte module: stable because RBC genes are absent from immune composition](findings/2026-10-14-darkred-is-erythrocyte-module.md) | E=4 R=0 B=2 D=4 C=5 | network/analyses/011-darkred-annotation |
| 11 | **20/25** | 7/10 | -1.0 | [Five experiments done; the decomposition is real, but the cross-group question is still unanswerable](findings/2026-10-25-preservation-arms-not-on-one-scale.md) | E=5 R=5 B=3 D=5 C=2 | network/analyses/017-preservation-synthesis |
| 12 | **20/25** | 5/10 | -3.0 | [The lithium-composition result is not circular: it survives removing every marker gene](findings/2026-10-14-lm22-circularity-passes.md) | E=5 R=4 B=2 D=5 C=4 | network/analyses/003-lm22-circularity |
| 13 | **19/25** | 7/10 | -0.6 | [PC2 is a lithium axis in the full cohort too, and the lithium DE signal is module-specific](findings/2026-10-14-pc2-lithium-axis-replicates-full-cohort.md) | E=4 R=4 B=3 D=4 C=4 | network/analyses/011-module-annotation |
| 14 | **18/25** | 7/10 | -0.2 | [Darkturquoise drives PC2: the lithium-associated dimension is localised to one module](findings/2026-10-14-darkturquoise-drives-pc2.md) | E=4 R=4 B=2 D=4 C=4 | network/analyses/005-pc2-module-loadings |
| 15 | **18/25** | 7/10 | -0.2 | [Darkturquoise group×composition interaction: the control-lith contrast survives proper slope correction](findings/2026-10-14-darkturquoise-interaction.md) | E=4 R=4 B=2 D=4 C=4 | network/analyses/008-darkturquoise-interaction |
| 16 | **18/25** | 6/10 | -1.2 | [No module eigengene associates with lithium — and the modules are not independent](findings/2026-10-14-module-trait-no-lithium-association.md) | E=4 R=5 B=2 D=4 C=3 | network/analyses/001-module-trait-ds4 |
| 17 | **18/25** | 5/10 | -2.2 | [base-002-ceil: split-half ceiling calibrates n=74 noise floor](findings/2026-10-24-base-002-ceil-ceiling.md) | E=4 R=4 B=3 D=3 C=4 | network/runs/base-002-ceil (offset-removed, 36-module reference on 160 controls, split-half ceiling) |
| 18 | **17/25** | 8/10 | +1.2 | [Darkturquoise: control-lith contrast survives composition adjustment; nolith-lith does not](findings/2026-10-14-darkturquoise-ctrl-lith-survives-adjustment.md) | E=4 R=4 B=2 D=3 C=4 | network/analyses/007-darkturquoise-across-groups |
| 19 | **16/25** | 8/10 | +1.6 | [BP disease tightens inter-module correlations: modules move more in lockstep](findings/2026-10-14-inter-module-correlations-tighten-in-bp.md) | E=3 R=4 B=2 D=0 C=4 | network/analyses/013-inter-module-correlations |
| 20 | **6/25** | 4/10 | +1.6 | [No module reliably beats a random gene set (degenerate decomposition)](findings/2026-10-09-gold-baseline-no-module-beats-random.md) | E=5 R=5 B=3 D=5 C=2 | network/runs/base-001 |
| 21 | **6/25** | 1/10 | -1.4 | [Cost scaling, measured](findings/2026-10-09-cost-scaling.md) | E=4 R=0 B=1 D=1 C=0 | Phase 0 probes + base-001 draw timings |
| 22 | **0/25** | 0/10 | +0.0 | [The merge threshold is decisive: 33 modules vs 3](findings/2026-10-09-merge-threshold-decisive.md) | E=5 R=5 B=3 D=5 C=4 | network/runs/mod-000-resolution |
| 23 | **0/25** | 0/10 | +0.0 | [Scale-free topology fails at every power and both sample sizes](findings/2026-10-09-scale-free-fails.md) | E=4 R=2 B=3 D=4 C=2 | network/runs/base-001 (sft curve) + Phase 0 probes |
| 24 | **0/25** | 0/10 | +0.0 | [BP disease drives network-wide module differences; lithium adds nothing beyond composition](findings/2026-10-14-bp-disease-drives-network-globally.md) | E=4 R=5 B=2 D=4 C=4 | network/analyses/009-all-modules-across-groups |
| 25 | **0/25** | 0/10 | +0.0 | [PC1 of the full-sample ME matrix captures BP disease: a global network factor](findings/2026-10-14-pc1-captures-bp-disease.md) | E=4 R=4 B=2 D=4 C=4 | network/analyses/010-full-pca-all-samples |

## Table B — Interest

| # | Interest | Evidence | Δ | finding |
|--:|--:|--:|--:|:--|
| 1 | **10/10** | 24/25 | +0.4 | [The reference network is built on an un-removed per-sample library offset](findings/2026-10-15-library-offset-artifact.md) |
| 2 | **10/10** | 23/25 | +0.8 | [The corrected network has 71 modules, and its eigengenes are near-orthogonal](findings/2026-10-18-rebuilt-network-71-modules.md) |
| 3 | **9/10** | 24/25 | -0.6 | [A lithium-associated co-expression dimension exists — and it is 73% cell composition](findings/2026-10-14-lithium-dimension-is-composition.md) |
| 4 | **9/10** | 23/25 | -0.2 | [inp-022: ILR residualization rescues bp_lith co-expression preservation](findings/2026-10-24-inp022-composition-rescue.md) |
| 5 | **9/10** | 22/25 | +0.2 | [Darkturquoise is a ribosomal/translation module: the composition signal has a biological substrate](findings/2026-10-14-darkturquoise-is-ribosomal.md) |
| 6 | **9/10** | 22/25 | +0.2 | [Removing the library offset rebuilds the network into a well-behaved one](findings/2026-10-18-normalisation-rebuild-required.md) |
| 7 | **9/10** | 21/25 | +0.6 | [Real modules do beat the random gene set — once the decomposition is real](findings/2026-10-18-modules-do-beat-gold.md) |
| 8 | **8/10** | 22/25 | -0.8 | [base-001 complete: modules beat gold under subsampling — but Zsummary is not comparable across arms](findings/2026-10-22-base001-complete-and-zsummary-not-comparable.md) |
| 9 | **8/10** | 20/25 | +0.0 | [Darkred is an erythrocyte module: stable because RBC genes are absent from immune composition](findings/2026-10-14-darkred-is-erythrocyte-module.md) |
| 10 | **8/10** | 17/25 | +1.2 | [Darkturquoise: control-lith contrast survives composition adjustment; nolith-lith does not](findings/2026-10-14-darkturquoise-ctrl-lith-survives-adjustment.md) |
| 11 | **8/10** | 16/25 | +1.6 | [BP disease tightens inter-module correlations: modules move more in lockstep](findings/2026-10-14-inter-module-correlations-tighten-in-bp.md) |
| 12 | **7/10** | 22/25 | -1.8 | [The lithium dimension replicates across decompositions, and never survives composition](findings/2026-10-14-lithium-dimension-replicates.md) |
| 13 | **7/10** | 20/25 | -1.0 | [Five experiments done; the decomposition is real, but the cross-group question is still unanswerable](findings/2026-10-25-preservation-arms-not-on-one-scale.md) |
| 14 | **7/10** | 19/25 | -0.6 | [PC2 is a lithium axis in the full cohort too, and the lithium DE signal is module-specific](findings/2026-10-14-pc2-lithium-axis-replicates-full-cohort.md) |
| 15 | **7/10** | 18/25 | -0.2 | [Darkturquoise drives PC2: the lithium-associated dimension is localised to one module](findings/2026-10-14-darkturquoise-drives-pc2.md) |
| 16 | **7/10** | 18/25 | -0.2 | [Darkturquoise group×composition interaction: the control-lith contrast survives proper slope correction](findings/2026-10-14-darkturquoise-interaction.md) |
| 17 | **6/10** | 18/25 | -1.2 | [No module eigengene associates with lithium — and the modules are not independent](findings/2026-10-14-module-trait-no-lithium-association.md) |
| 18 | **5/10** | 20/25 | -3.0 | [The lithium-composition result is not circular: it survives removing every marker gene](findings/2026-10-14-lm22-circularity-passes.md) |
| 19 | **5/10** | 18/25 | -2.2 | [base-002-ceil: split-half ceiling calibrates n=74 noise floor](findings/2026-10-24-base-002-ceil-ceiling.md) |
| 20 | **4/10** | 6/25 | +1.6 | [No module reliably beats a random gene set (degenerate decomposition)](findings/2026-10-09-gold-baseline-no-module-beats-random.md) |
| 21 | **1/10** | 6/25 | -1.4 | [Cost scaling, measured](findings/2026-10-09-cost-scaling.md) |
| 22 | **0/10** | 0/25 | +0.0 | [The merge threshold is decisive: 33 modules vs 3](findings/2026-10-09-merge-threshold-decisive.md) |
| 23 | **0/10** | 0/25 | +0.0 | [Scale-free topology fails at every power and both sample sizes](findings/2026-10-09-scale-free-fails.md) |
| 24 | **0/10** | 0/25 | +0.0 | [BP disease drives network-wide module differences; lithium adds nothing beyond composition](findings/2026-10-14-bp-disease-drives-network-globally.md) |
| 25 | **0/10** | 0/25 | +0.0 | [PC1 of the full-sample ME matrix captures BP disease: a global network factor](findings/2026-10-14-pc1-captures-bp-disease.md) |
