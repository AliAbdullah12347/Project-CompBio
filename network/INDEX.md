# Arm 3 — Master Index

> **GENERATED FILE — do not hand-edit.** Rebuild with `Rscript network/build_catalog.R`.
> On a merge conflict here, regenerate rather than merge. See `CONTRIBUTING.md`.

**Updated:** 2026-10-18 · **Catalogued:** 94 items

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

## Findings (20 — 16 live, 4 retracted)

| Evidence · Interest | status | finding | date |
|:--|:--|:--|:--|
| 24/25 · 9/10 | live | [A lithium-associated co-expression dimension exists — and it is 73% cell composition](findings/2026-10-14-lithium-dimension-is-composition.md) | 2026-10-14 |
| 24/25 · 10/10 | live | [The reference network is built on an un-removed per-sample library offset](findings/2026-10-15-library-offset-artifact.md) | 2026-10-15 |
| 23/25 · 10/10 | live | [The corrected network has 71 modules, and its eigengenes are near-orthogonal](findings/2026-10-18-rebuilt-network-71-modules.md) | 2026-10-18 |
| 22/25 · 9/10 | CORRECTED | [Darkturquoise is a ribosomal/translation module: the composition signal has a biological substrate](findings/2026-10-14-darkturquoise-is-ribosomal.md) | 2026-10-14 |
| 22/25 · 7/10 | live | [The lithium dimension replicates across decompositions, and never survives composition](findings/2026-10-14-lithium-dimension-replicates.md) | 2026-10-14 |
| 22/25 · 9/10 | live | [Removing the library offset rebuilds the network into a well-behaved one](findings/2026-10-18-normalisation-rebuild-required.md) | 2026-10-18 |
| 20/25 · 9/10 | live | [No module reliably beats a random gene set (degenerate decomposition)](findings/2026-10-09-gold-baseline-no-module-beats-random.md) | 2026-10-09 |
| 20/25 · 8/10 | live | [Darkred is an erythrocyte module: stable because RBC genes are absent from immune composition](findings/2026-10-14-darkred-is-erythrocyte-module.md) | 2026-10-14 |
| 20/25 · 5/10 | live | [The lithium-composition result is not circular: it survives removing every marker gene](findings/2026-10-14-lm22-circularity-passes.md) | 2026-10-14 |
| 19/25 · 7/10 | live | [PC2 is a lithium axis in the full cohort too, and the lithium DE signal is module-specific](findings/2026-10-14-pc2-lithium-axis-replicates-full-cohort.md) | 2026-10-14 |
| 18/25 · 7/10 | CORRECTED | [Darkturquoise drives PC2: the lithium-associated dimension is localised to one module](findings/2026-10-14-darkturquoise-drives-pc2.md) | 2026-10-14 |
| 18/25 · 7/10 | CORRECTED | [Darkturquoise group×composition interaction: the control-lith contrast survives proper slope correction](findings/2026-10-14-darkturquoise-interaction.md) | 2026-10-14 |
| 18/25 · 6/10 | live | [No module eigengene associates with lithium — and the modules are not independent](findings/2026-10-14-module-trait-no-lithium-association.md) | 2026-10-14 |
| 17/25 · 8/10 | live | [Darkturquoise: control-lith contrast survives composition adjustment; nolith-lith does not](findings/2026-10-14-darkturquoise-ctrl-lith-survives-adjustment.md) | 2026-10-14 |
| 16/25 · 8/10 | live | [BP disease tightens inter-module correlations: modules move more in lockstep](findings/2026-10-14-inter-module-correlations-tighten-in-bp.md) | 2026-10-14 |
| 6/25 · 1/10 | live | [Cost scaling, measured](findings/2026-10-09-cost-scaling.md) | 2026-10-09 |
| 0/25 · 0/10 | RETRACTED | [The merge threshold is decisive: 33 modules vs 3](findings/2026-10-09-merge-threshold-decisive.md) | 2026-10-09 |
| 0/25 · 0/10 | RETRACTED | [Scale-free topology fails at every power and both sample sizes](findings/2026-10-09-scale-free-fails.md) | 2026-10-09 |
| 0/25 · 0/10 | RETRACTED | [BP disease drives network-wide module differences; lithium adds nothing beyond composition](findings/2026-10-14-bp-disease-drives-network-globally.md) | 2026-10-14 |
| 0/25 · 0/10 | RETRACTED | [PC1 of the full-sample ME matrix captures BP disease: a global network factor](findings/2026-10-14-pc1-captures-bp-disease.md) | 2026-10-14 |

## Analyses (17)

| id | status | output | title |
|:--|:--|:--|:--|
| 001-module-trait-ds4 | has results | 5 result files | 001 — Module-trait relationships at the usable decomposition |
| 002-effective-dimensionality | has results | 3 result files | 002 — Effective dimensionality, and the arm's question answered |
| 003-lm22-circularity | has results | 4 result files | 003 — LM22 circularity check |
| 004-decomposition-robustness | has results | 3 result files | 004 — Does the lithium dimension replicate across decompositions? |
| 005-pc2-module-loadings | has results | 3 result files | 005-pc2-module-loadings |
| 006-darkturquoise-annotation | has results | 6 result files | 006-darkturquoise-annotation |
| 007-darkturquoise-across-groups | has results | 8 result files | 007-darkturquoise-across-groups |
| 008-darkturquoise-interaction | has results | 7 result files | 008-darkturquoise-interaction |
| 009-all-modules-across-groups | has results | 3 result files | 009-all-modules-across-groups |
| 010-full-pca-all-samples | has results | 5 result files | 010-full-pca-all-samples |
| 011-darkred-annotation | has results | 8 result files | 011-darkred-annotation |
| 011-module-annotation | has results | 5 result files | 005 — Module annotation, sample structure, hub genes |
| 012-full-module-annotation | ERROR | 4 result files | 012-full-module-annotation |
| 013-inter-module-correlations | has results | 8 result files | 013-inter-module-correlations |
| 014-effect-sizes | has results | 3 result files | 014-effect-sizes |
| 015-referee-verification | has results | 5 result files | 015-referee-verification |
| 016-normalisation-decision | has results | 2 result files | 016-normalisation-decision |

## Notes — infrastructure and defects (2)

| date | file | what happened |
|:--|:--|:--|
| 2026-10-15 | 2026-10-15-session-cap-starvation.md | Defect: rows whose cap equals SESSION_CAP can be BLOCKED while making progress |
| 2026-10-18 | 2026-10-18-three-days-lost-to-wifi.md | Three days lost: every wake started a session into a dead network |

## Run directories (4)

| run | status | size |
|:--|:--|:--|
| _raw_matrix_archive_20261018T211754Z | ARCHIVED (raw matrix) | 18 files |
| _raw_matrix_archive_20261018T211754Z | ARCHIVED (raw matrix) | 18 files |
| base-001 | live | 57 files |
| mod-000-resolution | live | 15 files |

## Queue experiments (51)

Status: DONE=1  PENDING=50

Full detail in `RUNLOG/queue.tsv`; one row per experiment in `catalog.csv`.

| id | family | status | spec |
|:--|:--|:--|:--|
| ann-041 | annot | PENDING | module overlap with lithium + bipolar DEG lists (de_final/results/); Fisher exact, BH across modules |
| ann-042 | annot | PENDING | module overlap with LM22 markers and cell-type DE results; Fisher exact, BH |
| base-001 | baseline | PENDING | ref=control;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions);stat=Zsumma |
| base-002-ceil | baseline | PENDING | SPLIT-HALF CEILING. ref=160 randomly held-out controls; test=the OTHER 74 controls (disjoint), plus bp_nolith(74) and bp |
| base-003-remerge | baseline | PENDING | THE KEY RE-RUN. Repeat base-001 exactly, but at the parameters mod-000-resolution found usable: deepSplit=4, mergeCutHei |
| bmind-046 | celltype | PENDING | REBUILD bMIND profiles: Rscript de_analysis/scripts/06_bmind_profiles.R (~3.1h). CHECKPOINT. Own window. |
| bmind-047 | celltype | PENDING | per-lineage networks from bMIND profiles (requires bmind-046 DONE) |
| cent-040 | central | PENDING | degree, betweenness, eigenvector centrality; between-group differences |
| comm-038 | community | PENDING | Louvain on adjacency; compare to WGCNA modules (ARI, NMI) |
| comm-039 | community | PENDING | Leiden on adjacency; compare to Louvain and to WGCNA (ARI, NMI) |
| corr-005 | corr | PENDING | ref=control;corr=pearson;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| corr-006 | corr | PENDING | ref=control;corr=spearman;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| filt-044 | filtsens | PENDING | reuse 31 filters from de_v2/scripts/13_filter_sweep_wb.R; which modules survive all? |
| gene-025 | geneset | PENDING | ref=control;genes=top2000_most_variable |
| gene-026 | geneset | PENDING | ref=control;genes=top5000_most_variable |
| gene-027 | geneset | PENDING | ref=control;genes=all_minus_547_LM22 ** CIRCULARITY CHECK ** |
| gene-028 | geneset | PENDING | ref=control;genes=lithium_DEGs_only (de_final/results/) |
| gene-029 | geneset | PENDING | ref=control;genes=explicit_non_DEGs |
| hub-036 | hub | PENDING | kME hub identity; stability across 100 draws; do hubs change between groups? |
| hub-037 | hub | PENDING | differential connectivity per gene between groups, same permutation null |
| inp-021 | input | PENDING | ref=control;input=residualised_technical(age+sex+rin+plate+seqpc1-3);draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (sp |
| inp-022 | input | PENDING | ref=control;input=residualised_composition(ILR b1-b4) ** KEY QUESTION **;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW |
| inp-023 | input | PENDING | ref=control;input=residualised_both(technical+ILR) ** KEY QUESTION **;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (s |
| inp-024 | input | PENDING | ref=control;input=voom_weighted_from_counts;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| mod-000-resolution | moddetect | DONE | BLOCKING DIAGNOSTIC - revised. base-001's log shows cutreeDynamic found 18 MODULES, which mergeCloseModules at mergeCutH |
| mod-013 | moddetect | PENDING | ref=control;deepSplit=0;minModuleSize=30;mergeCutHeight=0.25 |
| mod-014 | moddetect | PENDING | ref=control;deepSplit=2;minModuleSize=30;mergeCutHeight=0.25 |
| mod-015 | moddetect | PENDING | ref=control;deepSplit=4;minModuleSize=30;mergeCutHeight=0.25 |
| mod-016 | moddetect | PENDING | ref=control;deepSplit=2;minModuleSize=20;mergeCutHeight=0.25 |
| mod-017 | moddetect | PENDING | ref=control;deepSplit=2;minModuleSize=50;mergeCutHeight=0.25 |
| mod-018 | moddetect | PENDING | ref=control;deepSplit=2;minModuleSize=100;mergeCutHeight=0.25 |
| mod-019 | moddetect | PENDING | ref=control;deepSplit=2;minModuleSize=30;mergeCutHeight=0.15 |
| mod-020 | moddetect | PENDING | ref=control;deepSplit=2;minModuleSize=30;mergeCutHeight=0.35 |
| net-007 | nettype | PENDING | ref=control;corr=bicor;net=unsigned;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| net-008 | nettype | PENDING | ref=control;corr=bicor;net=signed_hybrid;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| perm-031 | perm | PENDING | group-label permutation null K=500; draws matched observed vs permuted; BH across MODULES not draws |
| perm-032 | perm | PENDING | group-label permutation null K=1000 if K=500 wall-clock allows |
| pow-009 | power | PENDING | pickSoftThreshold powerVector=1:20 per group; record scale-free R2 curves (curves ARE the result) |
| pow-010 | power | PENDING | ref=control;power=6 fixed;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans session |
| pow-011 | power | PENDING | ref=control;power=8 fixed;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans session |
| pow-012 | power | PENDING | ref=control;power=12 fixed;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessio |
| pres-030 | preserve | PENDING | ref=control;report density and connectivity components SEPARATELY; disagreement is informative |
| ref-002 | reference | PENDING | ref=bp_nolith;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| ref-003 | reference | PENDING | ref=bp_lith;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| ref-004 | reference | PENDING | ref=consensus(all3);corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions) |
| samp-043 | samplestruct | PENDING | eigengene-space separation by group; outlier detection; hclust; does any cluster track lithium? |
| synth-048 | synthesis | PENDING | cross-run synthesis: which conclusions are stable across all variants? write to RESULTS.md |
| tob-045 | tobacco | PENDING | repeat inp-022 across several of the 20 tobacco imputations; does any conclusion move? |
| trait-033 | modtrait | PENDING | eigengene ~ lithium, dx, age, sex, RIN, plate; BH across modules |
| trait-034 | modtrait | PENDING | eigengene ~ 5 lineage fractions (gran mono T NK B) |
| trait-035 | modtrait | PENDING | eigengene ~ ILR balances b1-b4 |

