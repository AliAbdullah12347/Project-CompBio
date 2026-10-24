# Arm 3 — Run State

**Updated:** 2026-10-24

## Where the science stands

**The arm's central question is answered provisionally.** Three reinforcing lines:

1. **The lithium network dimension is 73% cell composition** (finding #2, E=24/25). This
   replicates across all 10 tested decompositions (finding #7) and cannot be explained by
   circularity in the composition estimator (finding #10). It never survives composition adjustment.

2. **inp-022 (DONE 2026-10-24): ILR residualisation rescues bp_lith module preservation.**
   At matched instrument (50 perms/draw, 20 draws), gold Z rises from 4.06 (base-001 raw)
   to 15.55 after removing ILR b1–b4 (+284%). 41/44 (93%) modules beat gold mean Z; gold
   ranked last by medianRank in all 20 draws. The reference module count drops from 71 to
   44 — ~27 raw modules were composition artifacts. Finding:
   `findings/2026-10-24-inp022-composition-rescue.md` (E=23/25, I=9/10).

3. **base-002-ceil (DONE 2026-10-23): global co-expression differs between controls and
   nolith even at random-gene level.** At matched 500-perm instrument, gold Z=22.58 in
   disjoint controls vs 4.64 in bp_nolith — 5× difference. This is consistent with BP
   disease driving a global network shift (PC1=77.5%, finding #22).

**Read `RANKING.md` first.** Rebuilt from `findings/`; 24 findings as of 2026-10-24.

## Infrastructure

- Scheduler: **LaunchDaemon** `com.aaylab.networkarm.daemon` (system domain, survives
  logout; auth verified with a stripped environment). The earlier LaunchAgent died at the
  login window and cost 4 days 18 hours — see `METHODS.md` 2026-10-14.
- `com.aaylab.networkarm` (agent) retained as redundancy; `wake.sh`'s lock prevents overlap.
- Hourly watchdog `com.aaylab.networkarm.watchdog` kickstarts the agent if no wake is
  recorded for 6 h. Note it is a LaunchAgent and is therefore **also** dead when logged out;
  the daemon is the real protection.

## Queue

DONE (5): mod-000-resolution, base-001, base-003-remerge, base-002-ceil, inp-022.
RUNNING (1): **inp-023** (technical + ILR residualisation, claimed 2026-10-24T00:40Z).
PENDING: 45 rows.

Order after inp-023: pow-009 (SFT curves, fast), gene-027 (LM22 circularity check), inp-021
(technical-only residualisation), perm-031, trait-033/034/035, pres-030, samp-043, ann-041,
ann-042, ref-002/003/004, corr-005/006, net-007/008, pow-010/011/012, mod-013 through
mod-020, inp-024, gene-025/026/028/029, perm-032, hub-036/037, comm-038/039, cent-040,
filt-044, tob-045, samp-043, bmind-046/047 (require bMIND profiles, ~3.1 h to compute),
synth-048. Full order in `RUNLOG/queue.tsv`.

All 45 PENDING experiments have scripts. All scripts use `deepSplit=4, mergeCutHeight=0`
(builder scripts) or point to the v2 cache key `control_all_bicor_signed_p12_ds4_mch00`
(reader scripts). Fixed 2026-10-14.

## Parameters every downstream row should use

`power = 12`, `deepSplit = 4`, `mergeCutHeight = 0` → 33 modules, largest 11.68% of genes.
**Never `mergeCutHeight = 0.25`**: it collapses 33 modules to 3 with one holding 97.5% of
the transcriptome, because the median eigengene pair correlates 0.739 against a 0.75
threshold.

Quote every Zsummary against WGCNA's `gold` random-gene module, never against Langfelder's
published bands — gold scores 20.26 here, so the bands place noise in "strongly preserved".

## Analyses done (session 5–6, 2026-10-14)

- **001-module-trait-ds4**: No module eigengene associates with lithium (0/33, BH). Modules
  not independent (median |r|=0.822).
- **002-effective-dimensionality**: PC2 associates with lithium (β=-0.908, p=0.0015);
  73.3% cell composition; composition attenuates by 68.1%.
- **003-lm22-circularity**: Dropping all 350 LM22 markers: association and composition
  share move <1%. Not circular.
- **004-decomposition-robustness**: Lithium dimension (PC2) replicates in 6/10 decomposition
  cells, including 4/4 usable (≥10 modules, largest <25%). 0/6 survive composition
  adjustment. r(n_modules, composition share)=+0.879.
- **005-pc2-module-loadings**: Darkturquoise (140 genes, loading=0.455, r_pc2=0.843) is
  the only module with |loading|>0.3. PC2 is dominated by one module.
- **013-inter-module-correlations**: Inter-module correlation structure is largely preserved
  across groups (Mantel r=0.934, ctrl vs nolith), but BP disease TIGHTENS module coupling:
  median |r| ctrl=0.763 → nolith=0.790 → lith=0.845. Permutation p=0.024 shows this is
  more than random subsampling. Mechanistically explains PC1=77.5% variance. E=16/25, I=8/10.
- **011-darkred-annotation**: Darkred (164 genes) is an erythrocyte module (hub gene ANK1,
  ankyrin-1). LM22 has no erythrocyte cell type, so darkred has near-zero composition
  correlations (r_b2=−0.097 vs darkturquoise r_b2=−0.316). This mechanistically explains
  why darkred alone shows no group effect (FDR=0.981). E=20/25, I=8/10, Δ=+0.4.
- **006-darkturquoise-annotation**: Darkturquoise is a ribosomal/translation module (E=22/25,
- **007-darkturquoise-across-groups**: Extended to 460 samples (all groups). Control-lith
  contrast survives composition adjustment (p_bonf=0.001); nolith-lith does NOT (p_bonf=0.339).
  Composition-darkturquoise correlation weaker in controls (r_b2=−0.13) than BP (r_b2=−0.44).
- **010-full-pca-all-samples**: PC1=77.5% across 460 samples; all 33 modules load uniformly
  negative; r(PC1,disease)=0.289, perm p<0.001. PC2 (10%)=composition/lithium. Global
  network factor captures BP disease effect (ctrl mean PC1=-1.43, BP mean ≈+1.3).
- **009-all-modules-across-groups**: 31/33 ctrl-nolith sig; 32/33 ctrl-lith sig; 0/33
  nolith-lith sig (after additive composition adjustment, 460 samples). BP disease per se
  drives network change. Darkturquoise uniquely ctrl-lith (not ctrl-nolith) — composition.
  E=19/25, I=9/10, Δ=+1.4.
- **008-darkturquoise-interaction**: b2 slope differs between groups (F(2,451)=3.53, p=0.030);
  bp_lith vs control interaction = −0.041 (p=0.008). Control-lith contrast survives proper
  group-specific slope model (p_bonf=0.002). nolith-lith NS under all specs (p_bonf=0.509).
  I=9/10). 35 BP / 47 CC / 12 MF GO terms FDR<0.05. Top CC: cytosolic ribosome (29/92 genes,
  FDR=4.9e-32). Hub genes: RPL30, RPL35, RPL6, RPL27, RPS3A. Also contains Complex III/IV
  subunits. Story complete: high B-cell fraction → elevated ribosomal gene expression →
  darkturquoise eigengene high → PC2 positive → appears to track lithium (via composition).

## Pending task for next session

1. Check inp-023 DONE status — it mirrors inp-022's design (technical+ILR vs ILR-only)
2. Write finding for inp-023 comparing to inp-022 numbers
3. Rebuild `Rscript network/build_index.R`, commit, push
4. Run driver for next row (`bash network/run_next.sh`)

## Layout

New work goes in new files: `analyses/<NNN>-<slug>/`, `runs/<id>/`, `findings/<date>-<slug>.md`.
`RANKING.md`, `ranking.csv`, `results_index.csv` are **generated** — rebuild with
`Rscript network/build_index.R`, never hand-edit. See `CONTRIBUTING.md` for conflict
resolution from the laptop.
