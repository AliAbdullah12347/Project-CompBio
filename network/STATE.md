# Arm 3 — Run State

**Updated:** 2026-10-26T09:15Z

## Where the science stands

**The arm's central question is answered.** Current picture:

1. **The lithium network dimension is 73% cell composition** (finding #2, E=24/25). This
   replicates across all 10 tested decompositions (finding #7) and cannot be explained by
   circularity in the composition estimator (finding #10). It never survives composition adjustment.

2. **base-001 (DONE 2026-10-26): 71-module reference, corrected test matrices.**
   Gold Z nolith=14.01 (72/72 last, 65/71 above). Gold Z lith mean=21.74 (sd 7.05;
   72/72 last in all 20 draws; 17/71 beat gold in all 20 draws). See
   `findings/2026-10-26-base001-verified-clean.md` (E=23/25).

3. **inp-022 (DONE 2026-10-24): ILR residualisation within-arm gold always last.**
   Gold Z lith=15.55 (sd 8.13); 41/44 (93%) modules above gold; gold ranked last in
   all 20 draws. 44-module reference (27 raw modules were composition artifacts).
   `findings/2026-10-24-inp022-composition-rescue.md` (E=12/25 — cross-arm comparison
   invalid; within-arm gold comparison stands).

4. **inp-023 (PARTIAL, PENDING): technical + ILR residualisation, 7/20 draws.**
   Gold ranks 43/43 (last) in all 7 draws; mean 92.5% of modules beat gold (sd Z=6.84).
   Consistent with inp-022 within-arm pattern. Queued next.

5. **base-002-ceil (DONE 2026-10-26T07:46Z): split-half ceiling with corrected matrices.
   COMPLETE — 20/20 lith draws.** 67-module reference (n=160 controls).
   - ctrl_ceiling: gold Z=35.78, 68/68 last, 15/67 real modules above gold.
   - bp_nolith: gold Z=23.16, 68/68 last, 47/67 real modules above gold.
   - bp_lith (20 draws): gold Z mean=22.22 (sd=7.34, range 11.89–34.49), always 68/68 last,
     23/67 modules above gold in ALL 20 draws, mean 47.0/67 above gold per draw.
   - **KEY NEGATIVE:** lith and nolith are nearly indistinguishable (module Z correlation
     r=0.93; nolith median Z=31.05 vs lith mean Z=32.28). Both are below ceiling.
   `findings/2026-10-24-base-002-ceil-ceiling.md` (E=24/25).

**Read `RANKING.md` first.** Rebuilt from findings/; 30 findings as of 2026-10-26.

## Infrastructure

- Scheduler: **LaunchDaemon** `com.aaylab.networkarm.daemon` (system domain, survives
  logout; auth verified with a stripped environment). The earlier LaunchAgent died at the
  login window and cost 4 days 18 hours — see `METHODS.md` 2026-10-14.
- `com.aaylab.networkarm` (agent) retained as redundancy; `wake.sh`'s lock prevents overlap.
- Hourly watchdog `com.aaylab.networkarm.watchdog` kickstarts the agent if no wake is
  recorded for 6 h. Note it is a LaunchAgent and is therefore **also** dead when logged out;
  the daemon is the real protection.

## Queue

DONE (5): mod-000-resolution, inp-022, **base-001** (re-run 2026-10-26 with corrected
  expr_for() test matrices; 20/20 lith draws, 3603s), **base-003-remerge** (1s reader),
  **base-002-ceil** (DONE 2026-10-26T07:46Z; 20/20 lith draws, all groups complete, 5357s).
RUNNING (1): **inp-023** (started 2026-10-26T09:13Z; rebuilding TOM cache ~18 min, then
  draws 8–20 ×~12 min = ~156 min; expected done ~12:30Z; cap 12:13Z — may be PARTIAL).
  Bug fixed: assert_expr_mode incompatible with residualised matrices (commit 7bf5747).
PENDING (45): pow-009, gene-027,
  inp-021, perm-031, trait-033/034/035, pres-030 (**KEY** for base-002-ceil diagnosis),
  samp-043, ann-041, ann-042, ref-002/003/004, corr-005/006, net-007/008,
  pow-010/011/012, mod-013–mod-020, inp-024, gene-025/026/028/029, perm-032,
  hub-036/037, comm-038/039, cent-040, filt-044, tob-045, bmind-046/047, synth-048.

Key caches available:
- `ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm` (base-001) — EXISTS ✓
- `ref_control_all_bicor_signed_p12_composition_ilr_b1b4_ds4_mch00_offsetrm` (inp-022) — EXISTS ✓
- `ref_control_all_bicor_signed_p12_technical_ilr_age_sex_rin_b1b4_ds4_mch00_offsetrm` (inp-023) — EXISTS ✓
- ctrl160 reference (base-002-ceil) — EXISTS in cache ✓ (loaded at 03:17Z start; 67 modules)
- bp_nolith reference (ref-002) — NOT YET BUILT
- bp_lith reference (ref-003) — NOT YET BUILT
- technical-only (inp-021): needs age+sex+rin+plate+seqpc1-3 residualisation — NOT YET BUILT

All 45 PENDING experiments have scripts. All scripts use `deepSplit=4, mergeCutHeight=0`
(builder scripts) or point to the v2 cache key `control_all_bicor_signed_p12_ds4_mch00`
(reader scripts). Fixed 2026-10-14.

## Parameters every downstream row should use

`power = 12`, `deepSplit = 4`, `mergeCutHeight = 0`. With offset-removed expression from
`expr_for()` this yields **71 modules** in the n=234 control reference (base-001 verified),
NOT 33. The earlier "33 modules" came from mod-000-resolution which used the pre-offset-removal
expression. Never use 33 as the expected module count downstream — it is from a superseded analysis.
**Never `mergeCutHeight = 0.25`**: it collapses the structure with one module holding 97.5%
because median eigengene-pair correlation 0.739 is just below the 0.75 merge threshold.

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

1. **When inp-023 completes (20/20 draws):**
   - Score fully: expected E=5 R=5 B=5 D=5 C=4 = 24/25
   - Update inp-022 finding: R score 4→5, total 12→13/25 (robustness confirmed by replication)
   - IMPORTANT: do NOT conflate with pre-correction 23/25
2. **Driver continues:** pow-009 (SFT curves, ~15 min) → gene-027 (LM22, ~2.8h)
3. **pres-030** (density vs connectivity decomposition) is queued; run against the
   validated base-002-ceil results (attempt 4, corrected matrices). Do NOT use archived
   contaminated results.

## Session 7 notes (2026-10-25)

- inp-023 claimed at 19:16Z; process died after draw 7 (at 16:08 EDT); 7 draws kept
- base-001 claimed at 20:17Z UTC; re-ran draws 16-20 (each ~12 min); **DONE** 00:17Z
- base-003-remerge ran immediately after; **DONE** 1s (reader script only)
- base-002-ceil attempt 4 started 00:17Z; completed ceiling ctrl checkpoint by 21:57 local;
  cap at 03:17Z; requeued PARTIAL with 4 new checkpoints
- **base-002-ceil anomaly:** corrected in commit 1177d8c — raw expression was cause, not
  biology. METHODS.md and ceiling anomaly finding both updated.

## Session 8 notes (2026-10-26)

- base-002-ceil attempt 4 started 03:17Z; loaded 67-module cache (n=160 ref), ceiling
  checkpoint skipped (done in prior attempt); nolith 500 perms in progress
- base-002-ceil final attempt started 06:17Z; **DONE 07:46Z** (5357s). All 20 lith draws complete.
- KEY NEGATIVE: lith and nolith are nearly indistinguishable in module preservation
  (r=0.93, mean module Z: nolith=31.05 vs lith=32.28). Both well below ceiling (gold Z 23 vs 35.78).
- finding file updated; base-002-ceil scored E=24/25 (raised from 20 to 24 with complete data)
- inp-023 started at 09:13Z after assert_expr_mode bug fix (commit 7bf5747)
  TOM rebuilding (~18 min), then draws 8–20 (~156 min); cap 12:13Z

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
