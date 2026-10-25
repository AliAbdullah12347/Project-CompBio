# Arm 3 — Run State

**Updated:** 2026-10-25T20:10Z

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
   `findings/2026-10-24-inp022-composition-rescue.md` (**E=12/25 after correction** — the
   cross-arm Zsummary comparison is invalid due to different reference networks; within-arm
   gold comparison stands; see correction note in the finding file).

3. **inp-023 (PARTIAL, 2026-10-25 session 7): technical + ILR residualisation, 5/20 draws**
   done. Gold ranks 43/43 (last) in all 5 draws; mean 91% of modules beat gold (sd=7.8% Z).
   Consistent with inp-022 within-arm pattern. Draw 6 running; cap expires ~22:16 UTC.

4. **base-002-ceil anomaly (2026-10-25 — EXPLAINED, NOT BIOLOGY):** the anomaly was
   caused by a code bug: `base-002-ceil.R` used raw expression matrices for both reference
   and test, bypassing `expr_for()`. The "dense raw-graph signature" at n=160 vs n=74
   explains cor(size,Z) flipping from +0.911 to -0.036. Script fixed in commit 1177d8c;
   experiment reset to PENDING. See finding file and METHODS.md for full correction.

**Read `RANKING.md` first.** Rebuilt from `findings/`; 27 findings as of 2026-10-25 (including
partial inp-023 finding and corrected ceiling anomaly finding).

## Infrastructure

- Scheduler: **LaunchDaemon** `com.aaylab.networkarm.daemon` (system domain, survives
  logout; auth verified with a stripped environment). The earlier LaunchAgent died at the
  login window and cost 4 days 18 hours — see `METHODS.md` 2026-10-14.
- `com.aaylab.networkarm` (agent) retained as redundancy; `wake.sh`'s lock prevents overlap.
- Hourly watchdog `com.aaylab.networkarm.watchdog` kickstarts the agent if no wake is
  recorded for 6 h. Note it is a LaunchAgent and is therefore **also** dead when logged out;
  the daemon is the real protection.

## Queue

DONE (2): mod-000-resolution, inp-022.
PARTIAL (1): **inp-023** (technical + ILR residualisation, 2 prior attempts, 5/20 lith
  draws done; claimed session 7 at 2026-10-25T19:16Z; **42 modules** in doubly-residualised
  reference; draw 6 running; cap expires 22:16Z; draws 17-20 need next session).
PENDING: 48 rows (includes base-001, base-002-ceil, base-003-remerge — reset by commit
  1177d8c after raw test matrix bug found).

Order after inp-023: **base-001** (reset to PENDING — raw matrix bug fix, needs re-run
  with correct expr_for() test matrices; ~2.8h for 20 draws), then base-002-ceil (reset,
  ~6h), base-003-remerge (reads from base-001, fast ~1s), pow-009 (SFT curves, ~15 min),
  gene-027 (LM22 circularity, ~2.8h), inp-021, perm-031, trait-033/034/035, pres-030
  (**KEY** for base-002-ceil re-run diagnosis), samp-043, ann-041, ann-042, ref-002/003/004,
  corr-005/006, net-007/008, pow-010/011/012, mod-013 through mod-020, inp-024,
  gene-025/026/028/029, perm-032, hub-036/037, comm-038/039, cent-040, filt-044, tob-045,
  bmind-046/047 (require bMIND profiles, ~3.1 h to compute), synth-048. Full order in
  `RUNLOG/queue.tsv`.

Key caches available:
- `ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm` (base-001) — EXISTS ✓
- `ref_control_all_bicor_signed_p12_composition_ilr_b1b4_ds4_mch00_offsetrm` (inp-022) — EXISTS ✓
- `ref_control_all_bicor_signed_p12_technical_ilr_age_sex_rin_b1b4_ds4_mch00_offsetrm` (inp-023) — EXISTS ✓
- ctrl160 reference (base-002-ceil) — gitignored (will rebuild in ~27 min when needed)
- bp_nolith reference (ref-002) — NOT YET BUILT
- bp_lith reference (ref-003) — NOT YET BUILT
- technical-only (inp-021): needs age+sex+rin+plate+seqpc1-3 residualisation — NOT YET BUILT

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

After inp-023 PARTIAL (driver auto-commits after cap):
1. **Update finding** at `network/findings/2026-10-25-inp023-technical-ilr-replication.md`
   with stats from all draws done at session start (expected 5–16):
   - Run `Rscript network/cache/compute_inp023_stats.R` for current counts
   - Compare to inp-022: gold rank always last? ≥89% modules above gold?
2. **If 20/20 draws complete** (next session should finish 17-20, ~52 min):
   - Score fully: expected E=5 R=5 B=5 D=5 C=4 = 24/25
   - Update inp-022 finding: B score 4→5, total 12→13/25 (robustness confirmed)
     IMPORTANT: inp-022 total is 12/25 (cross-arm comparison invalid). B rising 4→5
     brings it to 13/25. Update carefully; do NOT conflate with pre-correction 23/25.
3. **Rebuild RANKING.md** via `Rscript network/build_index.R`
4. **Update RESULTS.md** with inp-023 completed run summary
5. **Commit + push**
6. **Run driver** `bash network/run_next.sh` → base-001 (reset, ~2.8h) → base-002-ceil
   (~6h) → base-003-remerge (~1s) → pow-009 (15 min) → gene-027 → inp-021

NOTE: base-002-ceil anomaly is EXPLAINED (code bug, not biology). The pres-030 follow-up
(density vs connectivity) is still valid for the re-run results when base-002-ceil
completes; do not run pres-030 against the archived contaminated results.

## inp-023 session 7 progress notes (2026-10-25)

- Driver claimed at 19:16Z UTC (15:16 EDT), cap expires 22:16Z (18:16 EDT)
- TOM cache loaded in ~18 min; nolith and draws 1-3 skipped (checkpointed)
- Draw 4 done at 15:29 EDT, draw 5 done ~15:42 EDT (~13 min each); draw 6 running
- 5-draw stats: gold Z mean=15.08 (sd=7.82), always last (5/5), 91% modules above gold
- Finding file: `findings/2026-10-25-inp023-technical-ilr-replication.md` (PARTIAL)
- **base-002-ceil anomaly corrected:** initial investigation found extraction code correct
  (both true). Concurrent commit 1177d8c found the real cause: raw expression matrices.
  METHODS.md entry and ceiling anomaly finding both corrected this session.

## inp-023 reference network (pre-computed, for finding reference)

**42 modules** in doubly-residualised reference (age+sex+rin + ILR b1–b4; deepSplit=4,
mergeCutHeight=0, power=12, bicor, signed). TOM cache written by attempt 1 and reused
by attempt 2 (TOM took 18 min in attempt 2 = loading cache + clustering).

Module count progression:
| arm | residualisation | modules | reduction vs raw |
|:--|:--|--:|:--|
| base-001 | none (raw) | 71 | — |
| inp-022 | ILR b1–b4 | 44 | −27 (38%) |
| inp-023 | age+sex+rin + ILR b1–b4 | **42** | −29 (41%) |

SFT curve comparison:
| arm | residualisation | first power R²≥0.80 | max R² | powerEstimate |
|:--|:--|--:|--:|--:|
| base-001 | none (raw) | never | 0.770 | NA |
| inp-022 | ILR b1–b4 | power 9 (R²=0.802) | 0.880 | 16 |
| inp-023 | age+sex+rin + ILR b1–b4 | power 13 (R²=0.806) | 0.836 | NA |

## Layout

New work goes in new files: `analyses/<NNN>-<slug>/`, `runs/<id>/`, `findings/<date>-<slug>.md`.
`RANKING.md`, `ranking.csv`, `results_index.csv` are **generated** — rebuild with
`Rscript network/build_index.R`, never hand-edit. See `CONTRIBUTING.md` for conflict
resolution from the laptop.
