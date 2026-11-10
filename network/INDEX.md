# Arm 3 — Master Index

> **GENERATED FILE — do not hand-edit.** Rebuild with `Rscript network/build_catalog.R`.
> On a merge conflict here, regenerate rather than merge. See `CONTRIBUTING.md`.

**Updated:** 2026-11-10 · **Catalogued:** 186 items

## How to search

```bash
bash network/find.sh darkturquoise      # search everything for a term
grep ',finding,' network/catalog.csv    # list findings
awk -F, '$4=="RETRACTED"' network/catalog.csv
```

| file | what it is |
|:--|:--|
| [INDEX.md](INDEX.md) | this file — entry point |
| [catalog.csv](catalog.csv) | every item, one row each, grep-friendly |
| [RANKING.md](RANKING.md) | findings by Evidence and by Interest |
| [RESULTS.md](RESULTS.md) | narrative write-up |
| [METHODS.md](METHODS.md) | decision log — why each choice was made |
| [STATE.md](STATE.md) | what is running right now |
| [CONTRIBUTING.md](CONTRIBUTING.md) | file conventions, merge-conflict rules |

## Findings (46 — 42 live, 4 retracted)

| Evidence · Interest | status | finding | date |
|:--|:--|:--|:--|
| 24/25 · 9/10 | live | [A lithium-associated co-expression dimension exists — and it is 73% cell composition](findings/2026-10-14-lithium-dimension-is-composition.md) | 2026-10-14 |
| 24/25 · 10/10 | live | [The reference network is built on an un-removed per-sample library offset](findings/2026-10-15-library-offset-artifact.md) | 2026-10-15 |
| 24/25 · 7/10 | live | [base-002-ceil: split-half ceiling calibration (re-run with verified offset-removed expression)](findings/2026-10-24-base-002-ceil-ceiling.md) | 2026-10-24 |
| 24/25 · 5/10 | live | [inp-023: technical + ILR residualisation replicates inp-022's rescue](findings/2026-10-25-inp023-technical-ilr-replication.md) | 2026-10-25 |
| 24/25 · 9/10 | live | [36 places scored a raw test matrix against an offset-removed reference](findings/2026-10-25-raw-test-matrices-bypassed-normalisation.md) | 2026-10-25 |
| 24/25 · 9/10 | live | [The central question, asked properly: control and bipolar-off-lithium networks are not detectably different](findings/2026-10-26-matched-comparison-clean-negative.md) | 2026-10-26 |
| 24/25 · 6/10 | live | [inp-021: technical-only residualisation — intermediate between raw and composition-adjusted](findings/2026-10-27-inp021-technical-residualisation.md) | 2026-10-27 |
| 24/25 · 6/10 | live | [trait-033: module eigengene ~ trait associations (lithium, dx, covariates)](findings/2026-10-31-trait033-module-eigengene-trait-associations.md) | 2026-10-31 |
| 24/25 · 7/10 | live | [trait-035: module eigengene ~ ILR balance associations](findings/2026-10-31-trait035-module-eigengene-ilr-associations.md) | 2026-10-31 |
| 23/25 · 10/10 | live | [The corrected network has 71 modules, and its eigengenes are near-orthogonal](findings/2026-10-18-rebuilt-network-71-modules.md) | 2026-10-18 |
| 23/25 · 8/10 | live | [base-001 re-run on a verified network: gold ranks last every time, but beats most modules on Zsummary](findings/2026-10-26-base001-verified-clean.md) | 2026-10-26 |
| 23/25 · 8/10 | live | [Resolved: lithium associations do not survive composition adjustment — finding 001 was right](findings/2026-11-08-composition-contradiction-resolved.md) | 2026-11-08 |
| 22/25 · 9/10 | CORRECTED | [Darkturquoise is a ribosomal/translation module: the composition signal has a biological substrate](findings/2026-10-14-darkturquoise-is-ribosomal.md) | 2026-10-14 |
| 22/25 · 7/10 | live | [The lithium dimension replicates across decompositions, and never survives composition](findings/2026-10-14-lithium-dimension-replicates.md) | 2026-10-14 |
| 22/25 · 9/10 | live | [Removing the library offset rebuilds the network into a well-behaved one](findings/2026-10-18-normalisation-rebuild-required.md) | 2026-10-18 |
| 22/25 · 7/10 | live | [Scale-free topology fails on all groups even after offset removal; composition removal rescues it](findings/2026-10-26-pow009-scale-free-composition-driven.md) | 2026-10-26 |
| 22/25 · 4/10 | live | [trait-034: module eigengene ~ lineage fraction associations](findings/2026-10-31-trait034-module-eigengene-fraction-associations.md) | 2026-10-31 |
| 22/25 · 8/10 | live | [perm-031 at 157 draws: the K=20 lead was a floor artifact, and there is no module-specific lithium effect](findings/2026-11-09-perm031-at-157-draws-is-a-clean-negative.md) | 2026-11-09 |
| 21/25 · 9/10 | CORRECTED | [Real modules do beat the random gene set — once the decomposition is real](findings/2026-10-18-modules-do-beat-gold.md) | 2026-10-18 |
| 21/25 · 6/10 | live | [LM22 exclusion does not change module preservation pattern: circularity check negative](findings/2026-10-26-gene027-lm22-circularity-preservation.md) | 2026-10-26 |
| 21/25 · 6/10 | live | [pres-030: density vs connectivity preservation components](findings/2026-10-31-pres030-density-connectivity-components.md) | 2026-10-31 |
| 21/25 · 9/10 | live | [The one module lithium moves that composition does not explain is the eosinophil module](findings/2026-11-09-eosinophil-module-survives-composition.md) | 2026-11-09 |
| 20/25 · 8/10 | live | [Darkred is an erythrocyte module: stable because RBC genes are absent from immune composition](findings/2026-10-14-darkred-is-erythrocyte-module.md) | 2026-10-14 |
| 20/25 · 5/10 | live | [The lithium-composition result is not circular: it survives removing every marker gene](findings/2026-10-14-lm22-circularity-passes.md) | 2026-10-14 |
| 20/25 · 7/10 | live | [perm-031: group-label permutation null for module preservation](findings/2026-10-27-perm031-group-label-permutation-null.md) | 2026-10-27 |
| 19/25 · 7/10 | live | [PC2 is a lithium axis in the full cohort too, and the lithium DE signal is module-specific](findings/2026-10-14-pc2-lithium-axis-replicates-full-cohort.md) | 2026-10-14 |
| 18/25 · 7/10 | CORRECTED | [Darkturquoise drives PC2: the lithium-associated dimension is localised to one module](findings/2026-10-14-darkturquoise-drives-pc2.md) | 2026-10-14 |
| 18/25 · 7/10 | CORRECTED | [Darkturquoise group×composition interaction: the control-lith contrast survives proper slope correction](findings/2026-10-14-darkturquoise-interaction.md) | 2026-10-14 |
| 18/25 · 6/10 | live | [No module eigengene associates with lithium — and the modules are not independent](findings/2026-10-14-module-trait-no-lithium-association.md) | 2026-10-14 |
| 18/25 · 4/10 | live | [ref-003: bp_lith (n=152) as reference — preservation in bp_nolith and control](findings/2026-11-02-ref003-bplith-reference.md) | 2026-11-02 |
| 17/25 · 8/10 | live | [Darkturquoise: control-lith contrast survives composition adjustment; nolith-lith does not](findings/2026-10-14-darkturquoise-ctrl-lith-survives-adjustment.md) | 2026-10-14 |
| 17/25 · 5/10 | live | [ann-042: module overlap with LM22 markers and cell-type lithium DEGs](findings/2026-11-01-ann042-module-lm22-celltype-overlap.md) | 2026-11-01 |
| 16/25 · 8/10 | live | [BP disease tightens inter-module correlations: modules move more in lockstep](findings/2026-10-14-inter-module-correlations-tighten-in-bp.md) | 2026-10-14 |
| 16/25 · 5/10 | live | [ann-041: module overlap with lithium / bipolar DEG lists](findings/2026-11-01-ann041-module-deg-overlap.md) | 2026-11-01 |
| 15/25 · 5/10 | CORRECTED | [The one fully-matched comparison behaves incoherently — stop before interpreting it](findings/2026-10-25-ceiling-vs-nolith-anomaly.md) | 2026-10-25 |
| 15/25 · 4/10 | live | [samp-043: eigengene-space sample structure by group](findings/2026-11-01-samp043-eigengene-sample-structure.md) | 2026-11-01 |
| 14/25 · 5/10 | live | [Five experiments done; the decomposition is real, but the cross-group question is still unanswerable](findings/2026-10-25-preservation-arms-not-on-one-scale.md) | 2026-10-25 |
| 14/25 · 6/10 | live | [cent-040's "52 of 71 modules differ" is one global density shift counted 52 times](findings/2026-11-10-centrality-difference-is-one-global-shift.md) | 2026-11-10 |
| 13/25 · 5/10 | CORRECTED | [inp-022: ILR residualization rescues bp_lith co-expression preservation](findings/2026-10-24-inp022-composition-rescue.md) | 2026-10-24 |
| 10/25 · 5/10 | CORRECTED | [base-001 complete: modules beat gold under subsampling — but Zsummary is not comparable across arms](findings/2026-10-22-base001-complete-and-zsummary-not-comparable.md) | 2026-10-22 |
| 6/25 · 1/10 | live | [Cost scaling, measured](findings/2026-10-09-cost-scaling.md) | 2026-10-09 |
| 6/25 · 4/10 | live | [No module reliably beats a random gene set (degenerate decomposition)](findings/2026-10-09-gold-baseline-no-module-beats-random.md) | 2026-10-09 |
| 0/25 · 0/10 | RETRACTED | [The merge threshold is decisive: 33 modules vs 3](findings/2026-10-09-merge-threshold-decisive.md) | 2026-10-09 |
| 0/25 · 0/10 | RETRACTED | [Scale-free topology fails at every power and both sample sizes](findings/2026-10-09-scale-free-fails.md) | 2026-10-09 |
| 0/25 · 0/10 | RETRACTED | [BP disease drives network-wide module differences; lithium adds nothing beyond composition](findings/2026-10-14-bp-disease-drives-network-globally.md) | 2026-10-14 |
| 0/25 · 0/10 | RETRACTED | [PC1 of the full-sample ME matrix captures BP disease: a global network factor](findings/2026-10-14-pc1-captures-bp-disease.md) | 2026-10-14 |

## Analyses (21)

| id | status | output | title |
|:--|:--|:--|:--|
| 001-module-trait-ds4 | has results | 3 result files | 001 — Module-trait relationships at the usable decomposition |
| 002-effective-dimensionality | has results | 2 result files | 002 — Effective dimensionality, and the arm's question answered |
| 003-lm22-circularity | has results | 3 result files | 003 — LM22 circularity check |
| 004-decomposition-robustness | has results | 2 result files | 004 — Does the lithium dimension replicate across decompositions? |
| 005-pc2-module-loadings | has results | 3 result files | 005-pc2-module-loadings |
| 006-darkturquoise-annotation | has results | 6 result files | 006-darkturquoise-annotation |
| 007-darkturquoise-across-groups | has results | 8 result files | 007-darkturquoise-across-groups |
| 008-darkturquoise-interaction | has results | 7 result files | 008-darkturquoise-interaction |
| 009-all-modules-across-groups | has results | 3 result files | 009-all-modules-across-groups |
| 010-full-pca-all-samples | has results | 5 result files | 010-full-pca-all-samples |
| 011-darkred-annotation | has results | 8 result files | 011-darkred-annotation |
| 011-module-annotation | has results | 4 result files | 005 — Module annotation, sample structure, hub genes |
| 012-full-module-annotation | ERROR | 4 result files | 012-full-module-annotation |
| 013-inter-module-correlations | has results | 8 result files | 013-inter-module-correlations |
| 014-effect-sizes | has results | 3 result files | 014-effect-sizes |
| 015-referee-verification | has results | 2 result files | 015-referee-verification |
| 016-normalisation-decision | has results | 2 result files | 016-normalisation-decision |
| 017-preservation-synthesis | has results | 2 result files | 017-preservation-synthesis |
| 018-matched-ceiling-comparison | has results | 1 result files | 018-matched-ceiling-comparison |
| 019-bpd-module-connection | has results | 2 result files | 019-bpd-module-connection |
| inp023_final_analysis.R | no results | 0 result files | inp023_final_analysis.R |

## Notes — infrastructure and defects (6)

| date | file | what happened |
|:--|:--|:--|
| 2026-10-15 | 2026-10-15-session-cap-starvation.md | Defect: rows whose cap equals SESSION_CAP can be BLOCKED while making progress |
| 2026-10-18 | 2026-10-18-logout-kills-eduroam.md | The real cause of the 3-day outage: logging out disconnects eduroam |
| 2026-10-18 | 2026-10-18-three-days-lost-to-wifi.md | Three days lost: every wake started a session into a dead network |
| 2026-10-25 | 2026-10-25-reference-build-verification.md | Reference builds now verify the matrix, not the label |
| 2026-10-26 | 2026-10-26-how-to-read-the-preservation-rows.md | How to read the 28 preservation rows: the perm split is by design, the comparison was not |
| 2026-11-08 | 2026-11-08-stale-running-path-ignored-progress.md | Three experiments were retired while holding 57 checkpoints between them |

## Run directories (61)

| run | status | size |
|:--|:--|:--|
| _raw_matrix_archive_20261018T211754Z | ARCHIVED (raw matrix) | 5 files |
| _raw_matrix_archive_20261018T211754Z | ARCHIVED (raw matrix) | 5 files |
| _raw_test_contaminated_20261025T193947Z | live | 59 files |
| ann-041 | live | 2 files |
| ann-042 | live | 2 files |
| ann-044 | live | 6 files |
| base-001 | live | 37 files |
| base-002-ceil | live | 6 files |
| base-003-remerge | live | 6 files |
| bmind-046 | live | 3 files |
| bmind-047 | live | 108 files |
| cent-040 | live | 5 files |
| comm-038 | live | 7 files |
| comm-039 | live | 0 files |
| corr-005 | live | 7 files |
| corr-006 | live | 49 files |
| driver_detached.log | live | 0 files |
| driver_session.out | live | 0 files |
| eos-052 | live | 6 files |
| filt-044 | live | 4 files |
| gene-025 | live | 49 files |
| gene-026 | live | 49 files |
| gene-027 | live | 7 files |
| gene-028 | live | 49 files |
| gene-029 | live | 49 files |
| grp-050 | live | 16 files |
| grp-051 | live | 27 files |
| hub-036 | live | 27 files |
| hub-037 | live | 23 files |
| inp-021 | live | 6 files |
| inp-022 | live | 15 files |
| inp-023 | live | 6 files |
| inp-024 | live | 48 files |
| mod-000-resolution | live | 2 files |
| mod-013 | live | 49 files |
| mod-014 | live | 49 files |
| mod-015 | live | 49 files |
| mod-016 | live | 49 files |
| mod-017 | live | 49 files |
| mod-018 | live | 49 files |
| mod-019 | live | 49 files |
| mod-020 | live | 49 files |
| net-007 | live | 49 files |
| net-008 | live | 49 files |
| perm-031 | live | 143 files |
| perm-032 | live | 44 files |
| pow-009 | live | 2 files |
| pow-010 | live | 49 files |
| pow-011 | live | 49 files |
| pow-012 | live | 49 files |
| pres-030 | live | 4 files |
| ref-002 | live | 87 files |
| ref-003 | live | 6 files |
| ref-004 | live | 50 files |
| samp-043 | live | 8 files |
| synth-048 | live | 4 files |
| tob-045 | live | 13 files |
| trait-033 | live | 3 files |
| trait-034 | live | 2 files |
| trait-035 | live | 2 files |
| trait-053 | live | 2 files |

## Queue experiments (52)

Status: DONE=48  RUNNING=4

Full detail in `RUNLOG/queue.tsv`; one row per experiment in `catalog.csv`.

| id | family | status | spec |
|:--|:--|:--|:--|
| ann-041 | annot | DONE | module overlap with lithium + bipolar DEG lists (de_final/results/); Fisher exact, BH across modules |
| ann-042 | annot | DONE | module overlap with LM22 markers and cell-type DE results; Fisher exact, BH |
| ann-044 | geneset | DONE | RE-RUN THE MODULE ANNOTATION ON THE CURRENT NETWORK. analyses/012-full-module-annotation ran GO over-representation on t |
| base-001 | baseline | DONE | ref=control;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions);stat=Zsumma |
| base-002-ceil | baseline | DONE | SPLIT-HALF CEILING. ref=160 randomly held-out controls; test=the OTHER 74 controls (disjoint), plus bp_nolith(74) and bp |
| base-003-remerge | baseline | DONE | THE KEY RE-RUN. Repeat base-001 exactly, but at the parameters mod-000-resolution found usable: deepSplit=4, mergeCutHei |
| bmind-046 | celltype | DONE | REBUILD bMIND profiles: Rscript de_analysis/scripts/06_bmind_profiles.R (~3.1h). CHECKPOINT. Own window. |
| bmind-047 | celltype | RUNNING | per-lineage networks from bMIND profiles (requires bmind-046 DONE) |
| cent-040 | central | DONE | degree, betweenness, eigenvector centrality; between-group differences |
| comm-038 | community | RUNNING | Louvain on adjacency; compare to WGCNA modules (ARI, NMI) |
| comm-039 | community | RUNNING | Leiden on adjacency; compare to Louvain and to WGCNA (ARI, NMI) |
| corr-005 | corr | DONE | ref=control;corr=pearson;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| corr-006 | corr | DONE | ref=control;corr=spearman;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| filt-044 | filtsens | DONE | reuse 31 filters from de_v2/scripts/13_filter_sweep_wb.R; which modules survive all? |
| gene-025 | geneset | DONE | ref=control;genes=top2000_most_variable |
| gene-026 | geneset | DONE | ref=control;genes=top5000_most_variable |
| gene-027 | geneset | DONE | ref=control;genes=all_minus_547_LM22 ** CIRCULARITY CHECK ** |
| gene-028 | geneset | DONE | ref=control;genes=lithium_DEGs_only (de_final/results/) |
| gene-029 | geneset | DONE | ref=control;genes=explicit_non_DEGs |
| hub-036 | hub | DONE | kME hub identity; stability across 100 draws; do hubs change between groups? |
| hub-037 | hub | RUNNING | differential connectivity per gene between groups, same permutation null |
| inp-021 | input | DONE | ref=control;input=residualised_technical(age+sex+rin+plate+seqpc1-3);draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (sp |
| inp-022 | input | DONE | ref=control;input=residualised_composition(ILR b1-b4) ** KEY QUESTION **;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW |
| inp-023 | input | DONE | ref=control;input=residualised_both(technical+ILR) ** KEY QUESTION **;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (s |
| inp-024 | input | DONE | ref=control;input=voom_weighted_from_counts;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| mod-000-resolution | moddetect | DONE | BLOCKING DIAGNOSTIC - revised. base-001's log shows cutreeDynamic found 18 MODULES, which mergeCloseModules at mergeCutH |
| mod-013 | moddetect | DONE | ref=control;deepSplit=0;minModuleSize=30;mergeCutHeight=0.25 |
| mod-014 | moddetect | DONE | ref=control;deepSplit=2;minModuleSize=30;mergeCutHeight=0.25 |
| mod-015 | moddetect | DONE | ref=control;deepSplit=4;minModuleSize=30;mergeCutHeight=0.25 |
| mod-016 | moddetect | DONE | ref=control;deepSplit=2;minModuleSize=20;mergeCutHeight=0.25 |
| mod-017 | moddetect | DONE | ref=control;deepSplit=2;minModuleSize=50;mergeCutHeight=0.25 |
| mod-018 | moddetect | DONE | ref=control;deepSplit=2;minModuleSize=100;mergeCutHeight=0.25 |
| mod-019 | moddetect | DONE | ref=control;deepSplit=2;minModuleSize=30;mergeCutHeight=0.15 |
| mod-020 | moddetect | DONE | ref=control;deepSplit=2;minModuleSize=30;mergeCutHeight=0.35 |
| net-007 | nettype | DONE | ref=control;corr=bicor;net=unsigned;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| net-008 | nettype | DONE | ref=control;corr=bicor;net=signed_hybrid;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| perm-031 | perm | DONE | group-label permutation null K=500; draws matched observed vs permuted; BH across MODULES not draws |
| perm-032 | perm | DONE | group-label permutation null K=1000 if K=500 wall-clock allows |
| pow-009 | power | DONE | pickSoftThreshold powerVector=1:20 per group; record scale-free R2 curves (curves ARE the result) |
| pow-010 | power | DONE | ref=control;power=6 fixed;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans session |
| pow-011 | power | DONE | ref=control;power=8 fixed;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans session |
| pow-012 | power | DONE | ref=control;power=12 fixed;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessio |
| pres-030 | preserve | DONE | ref=control;report density and connectivity components SEPARATELY; disagreement is informative |
| ref-002 | reference | DONE | ref=bp_nolith;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| ref-003 | reference | DONE | ref=bp_lith;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| ref-004 | reference | DONE | ref=consensus(all3);corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| samp-043 | samplestruct | DONE | eigengene-space separation by group; outlier detection; hclust; does any cluster track lithium? |
| synth-048 | synthesis | DONE | cross-run synthesis: which conclusions are stable across all variants? write to RESULTS.md |
| tob-045 | tobacco | DONE | repeat inp-022 across several of the 20 tobacco imputations; does any conclusion move? |
| trait-033 | modtrait | DONE | eigengene ~ lithium, dx, age, sex, RIN, plate; BH across modules |
| trait-034 | modtrait | DONE | eigengene ~ 5 lineage fractions (gran mono T NK B) |
| trait-035 | modtrait | DONE | eigengene ~ ILR balances b1-b4 |

