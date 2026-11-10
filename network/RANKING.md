# Arm 3 — Live Result Ranking

> **GENERATED FILE — do not hand-edit.** Rebuild with `Rscript network/build_index.R`.
> Edit the per-finding files in `network/findings/` instead. On a merge conflict here,
> regenerate rather than merge. See `CONTRIBUTING.md`.

**Findings:** 45   **Updated:** 2026-11-10

Two rankings. **Table A (Evidence, 0–25)** governs the write-up; it rewards evidential
strength so a well-powered null can score top. **Table B (Interest, 0–10)** records how
striking a result is. **Δ = (Interest/10 − Evidence/25) × 10**: positive means the story
is ahead of the proof.

## Table A — Evidence

| # | Evidence | Interest | Δ | finding | components | source |
|--:|--:|--:|--:|:--|:--|:--|
| 1 | **24/25** | 10/10 | +0.4 | [The reference network is built on an un-removed per-sample library offset](findings/2026-10-15-library-offset-artifact.md) | E=5 R=5 B=5 D=5 C=4 | network/analyses/015-referee-verification |
| 2 | **24/25** | 9/10 | -0.6 | [A lithium-associated co-expression dimension exists — and it is 73% cell composition](findings/2026-10-14-lithium-dimension-is-composition.md) | E=5 R=5 B=5 D=5 C=4 | network/analyses/002-effective-dimensionality |
| 3 | **24/25** | 9/10 | -0.6 | [36 places scored a raw test matrix against an offset-removed reference](findings/2026-10-25-raw-test-matrices-bypassed-normalisation.md) | E=5 R=5 B=5 D=5 C=4 | network/lint_expr.sh, network/00_load.R |
| 4 | **24/25** | 9/10 | -0.6 | [The central question, asked properly: control and bipolar-off-lithium networks are not detectably different](findings/2026-10-26-matched-comparison-clean-negative.md) | E=5 R=5 B=5 D=5 C=4 | network/analyses/018-matched-ceiling-comparison |
| 5 | **24/25** | 7/10 | -2.6 | [base-002-ceil: split-half ceiling calibration (re-run with verified offset-removed expression)](findings/2026-10-24-base-002-ceil-ceiling.md) | E=5 R=4 B=5 D=5 C=5 | network/runs/base-002-ceil (attempt 4; offset-removed expr_for(); 67-module reference on 160 controls, split-half ceiling) |
| 6 | **24/25** | 7/10 | -2.6 | [trait-035: module eigengene ~ ILR balance associations](findings/2026-10-31-trait035-module-eigengene-ilr-associations.md) | E=5 R=5 B=4 D=5 C=5 | network/runs/trait-035 — eigengene ~ b1 + b2 + b3 + b4 (ILR balances), all 474 samples, 71 modules |
| 7 | **24/25** | 6/10 | -3.6 | [inp-021: technical-only residualisation — intermediate between raw and composition-adjusted](findings/2026-10-27-inp021-technical-residualisation.md) | E=5 R=5 B=4 D=5 C=5 | network/runs/inp-021 — reference rebuilt on 234 controls residualised for age+sex+rin+plate+seqpc1-3 |
| 8 | **24/25** | 6/10 | -3.6 | [trait-033: module eigengene ~ trait associations (lithium, dx, covariates)](findings/2026-10-31-trait033-module-eigengene-trait-associations.md) | E=5 R=5 B=4 D=5 C=5 | network/runs/trait-033 — eigengene ~ lithium + dx + age + sex + rin + plate, all 474 samples, 71 modules |
| 9 | **24/25** | 5/10 | -4.6 | [inp-023: technical + ILR residualisation replicates inp-022's rescue](findings/2026-10-25-inp023-technical-ilr-replication.md) | E=5 R=5 B=5 D=5 C=4 | network/runs/inp-023 (age+sex+rin + ILR b1–b4 residualised, 42-module reference on 234 controls, 20 lith draws) |
| 10 | **23/25** | 10/10 | +0.8 | [The corrected network has 71 modules, and its eigengenes are near-orthogonal](findings/2026-10-18-rebuilt-network-71-modules.md) | E=5 R=5 B=5 D=5 C=3 | network/runs/mod-000-resolution (rebuilt on offset-removed expression) |
| 11 | **23/25** | 8/10 | -1.2 | [base-001 re-run on a verified network: gold ranks last every time, but beats most modules on Zsummary](findings/2026-10-26-base001-verified-clean.md) | E=5 R=5 B=4 D=5 C=4 | network/runs/base-001 — first run with reference **and** test both verified offset-removed |
| 12 | **23/25** | 8/10 | +0.8 | [Resolved: lithium associations do not survive composition adjustment — finding 001 was right](findings/2026-11-08-composition-contradiction-resolved.md) | E=5 R=5 B=4 D=5 C=4 | network/experiments/trait-053-composition-reconcile.R → network/runs/trait-053 |
| 13 | **22/25** | 9/10 | +0.3 | [Darkturquoise is a ribosomal/translation module: the composition signal has a biological substrate](findings/2026-10-14-darkturquoise-is-ribosomal.md) | E=5 R=5 B=2 D=5 C=5 | network/analyses/006-darkturquoise-annotation |
| 14 | **22/25** | 9/10 | +0.2 | [Removing the library offset rebuilds the network into a well-behaved one](findings/2026-10-18-normalisation-rebuild-required.md) | E=5 R=5 B=4 D=5 C=3 | network/analyses/016-normalisation-decision |
| 15 | **22/25** | 8/10 | +0.9 | [perm-031 at 157 draws: the K=20 lead was a floor artifact, and there is no module-specific lithium effect](findings/2026-11-09-perm031-at-157-draws-is-a-clean-negative.md) | E=? R=? B=? D=? C=? | network/experiments/perm-031.R → network/runs/perm-031 (157 draws) |
| 16 | **22/25** | 7/10 | -1.8 | [The lithium dimension replicates across decompositions, and never survives composition](findings/2026-10-14-lithium-dimension-replicates.md) | E=5 R=4 B=4 D=5 C=4 | network/analyses/004-decomposition-robustness |
| 17 | **22/25** | 7/10 | -1.8 | [Scale-free topology fails on all groups even after offset removal; composition removal rescues it](findings/2026-10-26-pow009-scale-free-composition-driven.md) | E=2 R=4 B=4 D=4 C=5 | network/runs/pow-009 (pickSoftThreshold, all 3 groups, offset-removed expression) |
| 18 | **22/25** | 4/10 | -4.8 | [trait-034: module eigengene ~ lineage fraction associations](findings/2026-10-31-trait034-module-eigengene-fraction-associations.md) | E=5 R=5 B=3 D=4 C=5 | network/runs/trait-034 — eigengene ~ mono + T + NK + B (gran dropped as reference), all 474 samples, 71 modules |
| 19 | **21/25** | 9/10 | +0.6 | [Real modules do beat the random gene set — once the decomposition is real](findings/2026-10-18-modules-do-beat-gold.md) | E=5 R=4 B=4 D=5 C=3 | network/runs/base-001 (offset-removed, 71-module reference) |
| 20 | **21/25** | 9/10 | +1.1 | [The one module lithium moves that composition does not explain is the eosinophil module](findings/2026-11-09-eosinophil-module-survives-composition.md) | E=? R=? B=? D=? C=? |  |
| 21 | **21/25** | 6/10 | -2.4 | [LM22 exclusion does not change module preservation pattern: circularity check negative](findings/2026-10-26-gene027-lm22-circularity-preservation.md) | E=5 R=4 B=3 D=5 C=4 | network/runs/gene-027 — reference rebuilt on 12018 non-LM22 genes (350 excluded) |
| 22 | **21/25** | 6/10 | -2.4 | [pres-030: density vs connectivity preservation components](findings/2026-10-31-pres030-density-connectivity-components.md) | E=5 R=4 B=3 D=4 C=5 | network/runs/pres-030 — modulePreservation density + connectivity components separately; 71 modules; bp_nolith (50 perms) + bp_lith (20 draws × 50 perms) |
| 23 | **20/25** | 8/10 | +0.0 | [Darkred is an erythrocyte module: stable because RBC genes are absent from immune composition](findings/2026-10-14-darkred-is-erythrocyte-module.md) | E=4 R=0 B=2 D=4 C=5 | network/analyses/011-darkred-annotation |
| 24 | **20/25** | 7/10 | -1.0 | [perm-031: group-label permutation null for module preservation](findings/2026-10-27-perm031-group-label-permutation-null.md) | E=5 R=5 B=2 D=3 C=5 | network/runs/perm-031 — K=20 permuted draws of 74 random samples from all 474, compared to base-001 bp_lith observed Zsummary |
| 25 | **20/25** | 5/10 | -3.0 | [The lithium-composition result is not circular: it survives removing every marker gene](findings/2026-10-14-lm22-circularity-passes.md) | E=5 R=4 B=2 D=5 C=4 | network/analyses/003-lm22-circularity |
| 26 | **19/25** | 7/10 | -0.6 | [PC2 is a lithium axis in the full cohort too, and the lithium DE signal is module-specific](findings/2026-10-14-pc2-lithium-axis-replicates-full-cohort.md) | E=4 R=4 B=3 D=4 C=4 | network/analyses/011-module-annotation |
| 27 | **18/25** | 7/10 | +0.2 | [Darkturquoise drives PC2: the lithium-associated dimension is localised to one module](findings/2026-10-14-darkturquoise-drives-pc2.md) | E=4 R=4 B=2 D=4 C=4 | network/analyses/005-pc2-module-loadings |
| 28 | **18/25** | 7/10 | -0.2 | [Darkturquoise group×composition interaction: the control-lith contrast survives proper slope correction](findings/2026-10-14-darkturquoise-interaction.md) | E=4 R=4 B=2 D=4 C=4 | network/analyses/008-darkturquoise-interaction |
| 29 | **18/25** | 6/10 | -1.2 | [No module eigengene associates with lithium — and the modules are not independent](findings/2026-10-14-module-trait-no-lithium-association.md) | E=4 R=5 B=2 D=4 C=3 | network/analyses/001-module-trait-ds4 |
| 30 | **18/25** | 4/10 | -3.2 | [ref-003: bp_lith (n=152) as reference — preservation in bp_nolith and control](findings/2026-11-02-ref003-bplith-reference.md) | E=4 R=3 B=4 D=3 C=4 | network/runs/ref-003 — 77 modules (+gold) defined once on all 152 bp_lith; tested in bp_nolith (n=74, one run, 50 perms) and control (20 draws of 74 from 234, 50 perms each). Power 12, bicor, signed. |
| 31 | **17/25** | 8/10 | +0.8 | [Darkturquoise: control-lith contrast survives composition adjustment; nolith-lith does not](findings/2026-10-14-darkturquoise-ctrl-lith-survives-adjustment.md) | E=4 R=4 B=2 D=3 C=4 | network/analyses/007-darkturquoise-across-groups |
| 32 | **17/25** | 5/10 | -1.4 | [ann-042: module overlap with LM22 markers and cell-type lithium DEGs](findings/2026-11-01-ann042-module-lm22-celltype-overlap.md) | E=3 R=5 B=3 D=3 C=3 | network/runs/ann-042 — Fisher exact test of 71 modules against 5 gene sets (LM22 signature genes, and cell-type lithium DEG sets T/mono/gran/NK), BH across modules within each set |
| 33 | **16/25** | 8/10 | +1.6 | [BP disease tightens inter-module correlations: modules move more in lockstep](findings/2026-10-14-inter-module-correlations-tighten-in-bp.md) | E=3 R=4 B=2 D=0 C=4 | network/analyses/013-inter-module-correlations |
| 34 | **16/25** | 5/10 | -1.4 | [ann-041: module overlap with lithium / bipolar DEG lists](findings/2026-11-01-ann041-module-deg-overlap.md) | E=3 R=4 B=2 D=3 C=4 | network/runs/ann-041 — Fisher exact test of 71 base-001 modules against three DEG lists (de_final results), BH across modules within each list |
| 35 | **15/25** | 5/10 | -1.0 | [The one fully-matched comparison behaves incoherently — stop before interpreting it](findings/2026-10-25-ceiling-vs-nolith-anomaly.md) | E=4 R=2 B=3 D=3 C=3 | network/runs/base-002-ceil |
| 36 | **15/25** | 4/10 | -2.0 | [samp-043: eigengene-space sample structure by group](findings/2026-11-01-samp043-eigengene-sample-structure.md) | E=3 R=4 B=2 D=2 C=4 | network/runs/samp-043 — 71 base-001 module eigengenes on all 474 samples; PCA, hclust k=3, centroid-distance outliers, per-module ANOVA by group |
| 37 | **14/25** | 5/10 | -1.0 | [Five experiments done; the decomposition is real, but the cross-group question is still unanswerable](findings/2026-10-25-preservation-arms-not-on-one-scale.md) | E=3 R=3 B=3 D=3 C=2 | network/analyses/017-preservation-synthesis |
| 38 | **13/25** | 5/10 | -0.2 | [inp-022: ILR residualization rescues bp_lith co-expression preservation](findings/2026-10-24-inp022-composition-rescue.md) | E=3 R=5 B=2 D=2 C=1 | network/runs/inp-022 (ILR b1–b4 residualized, 44-module reference on 234 controls, 20 lith draws) |
| 39 | **10/25** | 5/10 | +1.0 | [base-001 complete: modules beat gold under subsampling — but Zsummary is not comparable across arms](findings/2026-10-22-base001-complete-and-zsummary-not-comparable.md) | E=5 R=5 B=4 D=5 C=3 | network/runs/base-001 (offset-removed, 71-module reference, 20 draws) |
| 40 | **6/25** | 4/10 | +1.0 | [No module reliably beats a random gene set (degenerate decomposition)](findings/2026-10-09-gold-baseline-no-module-beats-random.md) | E=5 R=5 B=3 D=5 C=2 | network/runs/base-001 |
| 41 | **6/25** | 1/10 | -1.4 | [Cost scaling, measured](findings/2026-10-09-cost-scaling.md) | E=4 R=0 B=1 D=1 C=0 | Phase 0 probes + base-001 draw timings + gene-027 nolith timing |
| 42 | **0/25** | 0/10 | -1.8 | [The merge threshold is decisive: 33 modules vs 3](findings/2026-10-09-merge-threshold-decisive.md) | E=5 R=5 B=3 D=5 C=4 | network/runs/mod-000-resolution |
| 43 | **0/25** | 0/10 | +1.0 | [Scale-free topology fails at every power and both sample sizes](findings/2026-10-09-scale-free-fails.md) | E=4 R=2 B=3 D=4 C=2 | network/runs/base-001 (sft curve) + Phase 0 probes |
| 44 | **0/25** | 0/10 | +1.4 | [BP disease drives network-wide module differences; lithium adds nothing beyond composition](findings/2026-10-14-bp-disease-drives-network-globally.md) | E=4 R=5 B=2 D=4 C=4 | network/analyses/009-all-modules-across-groups |
| 45 | **0/25** | 0/10 | -0.2 | [PC1 of the full-sample ME matrix captures BP disease: a global network factor](findings/2026-10-14-pc1-captures-bp-disease.md) | E=4 R=4 B=2 D=4 C=4 | network/analyses/010-full-pca-all-samples |

## Table B — Interest

| # | Interest | Evidence | Δ | finding |
|--:|--:|--:|--:|:--|
| 1 | **10/10** | 24/25 | +0.4 | [The reference network is built on an un-removed per-sample library offset](findings/2026-10-15-library-offset-artifact.md) |
| 2 | **10/10** | 23/25 | +0.8 | [The corrected network has 71 modules, and its eigengenes are near-orthogonal](findings/2026-10-18-rebuilt-network-71-modules.md) |
| 3 | **9/10** | 24/25 | -0.6 | [A lithium-associated co-expression dimension exists — and it is 73% cell composition](findings/2026-10-14-lithium-dimension-is-composition.md) |
| 4 | **9/10** | 24/25 | -0.6 | [36 places scored a raw test matrix against an offset-removed reference](findings/2026-10-25-raw-test-matrices-bypassed-normalisation.md) |
| 5 | **9/10** | 24/25 | -0.6 | [The central question, asked properly: control and bipolar-off-lithium networks are not detectably different](findings/2026-10-26-matched-comparison-clean-negative.md) |
| 6 | **9/10** | 22/25 | +0.3 | [Darkturquoise is a ribosomal/translation module: the composition signal has a biological substrate](findings/2026-10-14-darkturquoise-is-ribosomal.md) |
| 7 | **9/10** | 22/25 | +0.2 | [Removing the library offset rebuilds the network into a well-behaved one](findings/2026-10-18-normalisation-rebuild-required.md) |
| 8 | **9/10** | 21/25 | +0.6 | [Real modules do beat the random gene set — once the decomposition is real](findings/2026-10-18-modules-do-beat-gold.md) |
| 9 | **9/10** | 21/25 | +1.1 | [The one module lithium moves that composition does not explain is the eosinophil module](findings/2026-11-09-eosinophil-module-survives-composition.md) |
| 10 | **8/10** | 23/25 | -1.2 | [base-001 re-run on a verified network: gold ranks last every time, but beats most modules on Zsummary](findings/2026-10-26-base001-verified-clean.md) |
| 11 | **8/10** | 23/25 | +0.8 | [Resolved: lithium associations do not survive composition adjustment — finding 001 was right](findings/2026-11-08-composition-contradiction-resolved.md) |
| 12 | **8/10** | 22/25 | +0.9 | [perm-031 at 157 draws: the K=20 lead was a floor artifact, and there is no module-specific lithium effect](findings/2026-11-09-perm031-at-157-draws-is-a-clean-negative.md) |
| 13 | **8/10** | 20/25 | +0.0 | [Darkred is an erythrocyte module: stable because RBC genes are absent from immune composition](findings/2026-10-14-darkred-is-erythrocyte-module.md) |
| 14 | **8/10** | 17/25 | +0.8 | [Darkturquoise: control-lith contrast survives composition adjustment; nolith-lith does not](findings/2026-10-14-darkturquoise-ctrl-lith-survives-adjustment.md) |
| 15 | **8/10** | 16/25 | +1.6 | [BP disease tightens inter-module correlations: modules move more in lockstep](findings/2026-10-14-inter-module-correlations-tighten-in-bp.md) |
| 16 | **7/10** | 24/25 | -2.6 | [base-002-ceil: split-half ceiling calibration (re-run with verified offset-removed expression)](findings/2026-10-24-base-002-ceil-ceiling.md) |
| 17 | **7/10** | 24/25 | -2.6 | [trait-035: module eigengene ~ ILR balance associations](findings/2026-10-31-trait035-module-eigengene-ilr-associations.md) |
| 18 | **7/10** | 22/25 | -1.8 | [The lithium dimension replicates across decompositions, and never survives composition](findings/2026-10-14-lithium-dimension-replicates.md) |
| 19 | **7/10** | 22/25 | -1.8 | [Scale-free topology fails on all groups even after offset removal; composition removal rescues it](findings/2026-10-26-pow009-scale-free-composition-driven.md) |
| 20 | **7/10** | 20/25 | -1.0 | [perm-031: group-label permutation null for module preservation](findings/2026-10-27-perm031-group-label-permutation-null.md) |
| 21 | **7/10** | 19/25 | -0.6 | [PC2 is a lithium axis in the full cohort too, and the lithium DE signal is module-specific](findings/2026-10-14-pc2-lithium-axis-replicates-full-cohort.md) |
| 22 | **7/10** | 18/25 | +0.2 | [Darkturquoise drives PC2: the lithium-associated dimension is localised to one module](findings/2026-10-14-darkturquoise-drives-pc2.md) |
| 23 | **7/10** | 18/25 | -0.2 | [Darkturquoise group×composition interaction: the control-lith contrast survives proper slope correction](findings/2026-10-14-darkturquoise-interaction.md) |
| 24 | **6/10** | 24/25 | -3.6 | [inp-021: technical-only residualisation — intermediate between raw and composition-adjusted](findings/2026-10-27-inp021-technical-residualisation.md) |
| 25 | **6/10** | 24/25 | -3.6 | [trait-033: module eigengene ~ trait associations (lithium, dx, covariates)](findings/2026-10-31-trait033-module-eigengene-trait-associations.md) |
| 26 | **6/10** | 21/25 | -2.4 | [LM22 exclusion does not change module preservation pattern: circularity check negative](findings/2026-10-26-gene027-lm22-circularity-preservation.md) |
| 27 | **6/10** | 21/25 | -2.4 | [pres-030: density vs connectivity preservation components](findings/2026-10-31-pres030-density-connectivity-components.md) |
| 28 | **6/10** | 18/25 | -1.2 | [No module eigengene associates with lithium — and the modules are not independent](findings/2026-10-14-module-trait-no-lithium-association.md) |
| 29 | **5/10** | 24/25 | -4.6 | [inp-023: technical + ILR residualisation replicates inp-022's rescue](findings/2026-10-25-inp023-technical-ilr-replication.md) |
| 30 | **5/10** | 20/25 | -3.0 | [The lithium-composition result is not circular: it survives removing every marker gene](findings/2026-10-14-lm22-circularity-passes.md) |
| 31 | **5/10** | 17/25 | -1.4 | [ann-042: module overlap with LM22 markers and cell-type lithium DEGs](findings/2026-11-01-ann042-module-lm22-celltype-overlap.md) |
| 32 | **5/10** | 16/25 | -1.4 | [ann-041: module overlap with lithium / bipolar DEG lists](findings/2026-11-01-ann041-module-deg-overlap.md) |
| 33 | **5/10** | 15/25 | -1.0 | [The one fully-matched comparison behaves incoherently — stop before interpreting it](findings/2026-10-25-ceiling-vs-nolith-anomaly.md) |
| 34 | **5/10** | 14/25 | -1.0 | [Five experiments done; the decomposition is real, but the cross-group question is still unanswerable](findings/2026-10-25-preservation-arms-not-on-one-scale.md) |
| 35 | **5/10** | 13/25 | -0.2 | [inp-022: ILR residualization rescues bp_lith co-expression preservation](findings/2026-10-24-inp022-composition-rescue.md) |
| 36 | **5/10** | 10/25 | +1.0 | [base-001 complete: modules beat gold under subsampling — but Zsummary is not comparable across arms](findings/2026-10-22-base001-complete-and-zsummary-not-comparable.md) |
| 37 | **4/10** | 22/25 | -4.8 | [trait-034: module eigengene ~ lineage fraction associations](findings/2026-10-31-trait034-module-eigengene-fraction-associations.md) |
| 38 | **4/10** | 18/25 | -3.2 | [ref-003: bp_lith (n=152) as reference — preservation in bp_nolith and control](findings/2026-11-02-ref003-bplith-reference.md) |
| 39 | **4/10** | 15/25 | -2.0 | [samp-043: eigengene-space sample structure by group](findings/2026-11-01-samp043-eigengene-sample-structure.md) |
| 40 | **4/10** | 6/25 | +1.0 | [No module reliably beats a random gene set (degenerate decomposition)](findings/2026-10-09-gold-baseline-no-module-beats-random.md) |
| 41 | **1/10** | 6/25 | -1.4 | [Cost scaling, measured](findings/2026-10-09-cost-scaling.md) |
| 42 | **0/10** | 0/25 | -1.8 | [The merge threshold is decisive: 33 modules vs 3](findings/2026-10-09-merge-threshold-decisive.md) |
| 43 | **0/10** | 0/25 | +1.0 | [Scale-free topology fails at every power and both sample sizes](findings/2026-10-09-scale-free-fails.md) |
| 44 | **0/10** | 0/25 | +1.4 | [BP disease drives network-wide module differences; lithium adds nothing beyond composition](findings/2026-10-14-bp-disease-drives-network-globally.md) |
| 45 | **0/10** | 0/25 | -0.2 | [PC1 of the full-sample ME matrix captures BP disease: a global network factor](findings/2026-10-14-pc1-captures-bp-disease.md) |
