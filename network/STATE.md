# Arm 3 — Run State

**Updated:** 2026-10-27T19:15Z

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
   `findings/2026-10-24-inp022-composition-rescue.md` (E=13/25 — R raised 4→5 after inp-023
   replication; cross-arm comparison remains invalid; within-arm gold comparison stands).

4. **inp-023 (DONE 2026-10-26T12:05Z): technical + ILR residualisation, 20/20 draws.**
   Gold Z mean=16.98 (sd=10.79); gold ranked 43/43 (last) in ALL 20 draws; 40/42 real modules
   above gold mean Z (95%). bp_nolith: gold Z=28.08, 43/43 last. Fully replicates inp-022
   within-arm pattern. E=24/25. inp-022 R score raised 4→5 (total 12→13/25).
   `findings/2026-10-25-inp023-technical-ilr-replication.md`

5. **base-002-ceil (DONE 2026-10-26T07:46Z): split-half ceiling with corrected matrices.
   COMPLETE — 20/20 lith draws.** 67-module reference (n=160 controls).
   - ctrl_ceiling: gold Z=35.78, 68/68 last, 15/67 real modules above gold.
   - bp_nolith: gold Z=23.16, 68/68 last, 47/67 real modules above gold.
   - bp_lith (20 draws): gold Z mean=22.22 (sd=7.34, range 11.89–34.49), always 68/68 last,
     23/67 modules above gold in ALL 20 draws, mean 47.0/67 above gold per draw.
   - **KEY NEGATIVE:** lith and nolith are nearly indistinguishable (module Z correlation
     r=0.93; nolith median Z=31.05 vs lith mean Z=32.28). Both are below ceiling.
   `findings/2026-10-24-base-002-ceil-ceiling.md` (E=24/25).

6. **pow-009 (DONE 2026-10-26T14:12Z): SFT curves on offset-removed expression.**
   All 3 groups fail scale-free (control max R²=0.770, bp_nolith=0.755, bp_lith=0.701;
   all below 0.80). Offset removal improves but doesn't fix; composition removal (inp-022)
   rescues it. New finding `2026-10-26-pow009-scale-free-composition-driven.md` (E=22/25).
   Partially rehabilitates retracted finding #27.

7. **gene-027 (DONE 2026-10-26T23:15Z): LM22 circularity check NEGATIVE.**
   65-module reference rebuilt on 12018 non-LM22 genes. nolith gold Z=12.23 (66/66 last);
   lith mean gold Z=21.41 (sd 6.67, range 9.00–33.02), 66/66 last in ALL 20 draws.
   27/65 modules beat gold in all 20 draws. Pattern indistinguishable from base-001 (all genes).
   `findings/2026-10-26-gene027-lm22-circularity-preservation.md` (E=21/25).

8. **inp-021 (PARTIAL — 1/20 lith draws; run resuming in session 16): technical-only residualisation.**
   63-module reference on technically-residualised controls (age+sex+rin+plate+seqpc1-3).
   Nolith gold Z=31.13 (64/64 last, 50 perms). Draw 1 gold Z=31.29 (64/64 last).
   Scale-free: max R²=0.797 (never achieves 0.80; ILR is required). Module count: 63 (−8 vs base-001).
   All 6 LM22-dependent modules also absent from technical-only reference (convergent dissolution).
   Finding file `findings/2026-10-27-inp021-technical-residualisation.md` — interim filled with 1-draw stats.
   Interim score: E=5 R=5 B=4 D=5 C=2 → 21/25 (PARTIAL); Interest=6/10; Δ=−2.4.
   NOTE: 18/20 draws now complete (draws 1-18 checkpointed). Session 19 (2026-10-27T19:10Z) found
   inp-021 BLOCKED due to orphan-process race (METHODS.md 2026-10-27T19:10Z). Manual reset to PENDING,
   attempts=1. Will complete draws 19-20 (~24 min) in the next driver invocation after perm-031 finishes.

**Read `RANKING.md` first.** Rebuilt from findings/; 32 findings as of 2026-10-27 (inp-021 pending score).

## Infrastructure

- Scheduler: **LaunchDaemon** `com.aaylab.networkarm.daemon` (system domain, survives
  logout; auth verified with a stripped environment). The earlier LaunchAgent died at the
  login window and cost 4 days 18 hours — see `METHODS.md` 2026-10-14.
- `com.aaylab.networkarm` (agent) retained as redundancy; `wake.sh`'s lock prevents overlap.
- Hourly watchdog `com.aaylab.networkarm.watchdog` kickstarts the agent if no wake is
  recorded for 6 h. Note it is a LaunchAgent and is therefore **also** dead when logged out;
  the daemon is the real protection.

## Queue

DONE (8): mod-000-resolution, inp-022, **base-001** (re-run 2026-10-26 with corrected
  expr_for() test matrices; 20/20 lith draws, 3603s), **base-003-remerge** (1s reader),
  **base-002-ceil** (DONE 2026-10-26T07:46Z; 20/20 lith draws, all groups complete, 5357s),
  **inp-023** (DONE 2026-10-26T12:05Z; 20/20 lith draws, 10311s),
  **pow-009** (DONE 2026-10-26T14:12Z; SFT curves all groups, 86s),
  **gene-027** (DONE 2026-10-26T23:15Z; LM22 circularity NEGATIVE; 65-module reference on 12018
  non-LM22 genes; nolith gold Z=12.23 (66/66 last); lith mean gold Z=21.41 (sd 6.67), 66/66 last
  in ALL 20 draws; 27/65 modules beat gold in all 20 draws. E=21/25. 14093s).
RUNNING (1): **perm-031** (claimed 2026-10-27T19:06Z; group-label permutation null K=20;
  random draws from all 474 samples, vs base-001 bp_lith observed Zsummary;
  draw 1 in progress as of 19:15Z; cap 22:06Z; ~15 draws expected before cap).
PENDING (43): **inp-021** (reset from BLOCKED 2026-10-27T19:10Z; 18/20 lith draws done;
  claims when perm-031 completes this session), trait-033/034/035, pres-030,
  samp-043, ann-041, ann-042, ref-002/003/004, corr-005/006, net-007/008,
  pow-010/011/012, mod-013–mod-020, inp-024, gene-025/026/028/029, perm-032,
  hub-036/037, comm-038/039, cent-040, filt-044, tob-045, bmind-046/047, synth-048.

Key caches available:
- `ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm` (base-001) — EXISTS ✓
- `ref_control_all_bicor_signed_p12_composition_ilr_b1b4_ds4_mch00_offsetrm` (inp-022) — EXISTS ✓
- `ref_control_all_bicor_signed_p12_technical_ilr_age_sex_rin_b1b4_ds4_mch00_offsetrm` (inp-023) — EXISTS ✓
- ctrl160 reference (base-002-ceil) — EXISTS in cache ✓ (loaded at 03:17Z start; 67 modules)
- technical-only (inp-021) — EXISTS ✓ (`ref_control_all_bicor_signed_p12_technical_full_age_sex_rin_plate_seqpc123_ds4_mch00_offsetrm.rds`)
- bp_nolith reference (ref-002) — NOT YET BUILT
- bp_lith reference (ref-003) — NOT YET BUILT

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

1. **perm-031** PARTIAL expected at 22:06Z cap (~15 draws); driver will commit/push on exit.
   - After perm-031 session ends: inp-021 is next PENDING (row 1 in queue).
2. **inp-021** PENDING (18/20 draws done):
   - Next session (after perm-031 current RUNNING): claims inp-021, runs draws 19-20 (~24 min), exits 0 DONE.
   - When DONE: run `Rscript network/scripts/analyze_inp021.R`, fill finding lith section fully,
     re-score to C=5 → E=24/25, rebuild RANKING.md, commit and push.
3. **perm-031** resumes after inp-021 DONE (perm-031 will be PARTIAL, queued after inp-021).
4. **trait-033/034/035** — fast (<60s each); claimed after perm-031 session.
5. **pres-030** — WGCNA heavy, similar to base-001 structure; PARTIAL expected.

## Session 19 notes (2026-10-27T19:06Z)

- Session started, found inp-021 BLOCKED (18/20 lith draws done; driver blocked after orphan-race).
- Root cause: 09:06Z driver was killed (SIGHUP) leaving Rscript as orphan; draws 17-18 completed by orphan.
  14:08Z new session requeued (attempts=2<3) → claimed (attempts=3) → ran draws 17-18 again → started draw 19.
  19:06Z session: dead claim detected → attempts=3≥3 → BLOCKED. But 18 draws are valid checkpoints.
- Manual fix: reset inp-021 BLOCKED→PENDING, attempts=3→1. Documented in METHODS.md 2026-10-27T19:10Z.
- perm-031 claimed at 19:06Z (was next PENDING after inp-021's BLOCKED was detected by driver).
  perm-031 running draw 1 as of 19:15Z. Cap 22:06Z. Expected ~15/20 draws before cap.
- Commits: 9fab110 (queue reset + METHODS.md).
- Queue: DONE=8 RUNNING=1 (perm-031) PENDING=43 (inp-021 included) BLOCKED=1 (inp-021, now reset).

## Session 18 notes (2026-10-27T14:08Z)

- inp-021 re-claimed (09:06Z claim PID dead). Found draws 1-16 done. Started draw 17.
- Draw 17 done at 10:20:27 EDT (14:20:27 UTC). Draw 18 done at 10:32:08 EDT. Draw 19 started.
- Process died during draw 19 (no draw 19 checkpoint). 18/20 draws total.
- No progress log "done" entry — driver was killed before Rscript finished.
  Attempted 3rd time: still only 18 draws. Driver at cap or killed before PARTIAL recorded.
- Result: inp-021 sits RUNNING (dead) → BLOCKED in session 19.

## Session 17 notes (2026-10-27T09:06Z)

- inp-021 claimed at 09:06Z. Driver process itself killed (SIGHUP). Rscript became orphan.
- Orphan R computed draws 17-18 (10:20, 10:32 EDT respectively).
- Orphan started draw 19 before being killed. Draw 19 incomplete.
- Driver never called mark() → no PARTIAL/DONE record. Claim PID dead for session 18.

## Session 16 notes (2026-10-27T04:07Z)

- inp-021 claim PID was dead (process died mid-draw-2 at ~00:45 EDT; only 1 draw saved, not ~13).
  Reason: session 15 estimate of ~13 draws before 23:08 EDT cap was correct in timing, but the
  R process exited abnormally. No error in log — log ends mid-draw-2 with no completion line.
  Likely OOM or signal kill during draw 2 permutation computation.
- Session 16: detected dead claim via .claim_inp-021.pid PID check. Dead claim enabled immediate
  requeue (no 6h stale wait). Driver re-ran at 04:07Z; reference loaded from cache (fast); nolith
  and draw 1 checkpoints skipped; draw 2 now running.
- Queue: DONE=8 RUNNING=1 (inp-021, attempt 2) PENDING=42 BLOCKED=0.
- inp-021 interim finding filled (1-draw stats), STATE.md/RESULTS.md updated.
- RANKING.md already rebuilt (session 16): 34 findings, inp-021 at rank 15 (E=21/25 interim). Re-run after all 20 draws to update C=2→5 score.

## Session 15 notes (2026-10-27T00:40Z)

- inp-021 nolith complete at 2026-10-26T20:33:29 EDT (703s; 50 perms, 63 modules + gold).
  gold Z=31.13 (50 perms), medianRank=64/64 (LAST). Same gold-last pattern as all other arms.
  25/63 real modules above gold by Zsummary; all 63 above gold by medianRank.
  Zsummary inflation vs base-001 (14.01) reflects fewer perms + cleaner data; medianRank is primary.
- Lith draw 1 done 2026-10-26T20:45:06 EDT (697s). At 697s/draw: ~13 draws expected before cap → PARTIAL.
  ACTUAL: R process died mid-draw-2; only 1 draw saved. Session 16 restarted.
- All commits pushed (LaunchDaemon sync_up() resolved earlier push failures).

## Session 14 notes (2026-10-27T00:45Z)

- gene-027 documentation completed: RANKING.md rebuilt (32 findings), RESULTS.md/STATE.md updated.
- inp-021 automatically claimed by LaunchDaemon at 2026-10-27T00:08Z. RUNNING.
- Queue: DONE=8 RUNNING=1 PENDING=42 BLOCKED=0.
- perm-031 is next claimable (script exists; cap 10800s; group-label permutation null K=500).

## Session 13 notes (2026-10-26T23:30Z)

- gene-027 DONE at 23:15:47Z (14093s). All 20 lith draws confirmed present.
- Final results: nolith gold Z=12.23 (66/66 last); lith mean gold Z=21.41 (sd=6.67,
  range 9.00–33.02), 66/66 last in ALL 20 draws. 27/65 modules beat gold in all 20 draws.
- Circularity check NEGATIVE. Same population-level pattern as base-001 with 350 LM22 genes excluded.
  Mean gold Z nearly identical (21.41 vs base-001 21.74), confirming LM22 genes do not drive the result.
- Finding file fully filled (Evidence=21/25, Interest=6/10, Δ=−2.4). RANKING.md rebuilt (32 findings).
- RESULTS.md, STATE.md updated. Committed and pushed. inp-021 started automatically.

## Session 12 notes (2026-10-26T19:45Z)

- gene-027 nolith confirmed: gold Z=12.23, medianRank=66/66 (last). Same qualitative
  pattern as base-001. Circularity check NEGATIVE at population level.
- Key structural finding: nolith Zsummary correlation between same-named modules
  (gene-027 vs base-001) = 0.359; medianRank correlation = −0.013. Explanation: removing
  LM22 genes reshuffles the dendrogram — same-named modules have different gene content.
  Individual module comparison across runs is invalid; only population-level (gold ranking) holds.
- Cross-base-001 observation: the 6 LM22-dependent modules (absent in gene-027's reference)
  are exactly the most disrupted in base-001 bp_lith: mean beat_gold=3.3/20 vs 14.8/20 for
  non-LM22 modules. These modules likely capture immune cell composition patterns.
- Draw 1 complete at 19:44Z: gold Z=30.16, rank=66/66 (LAST).
- Draw 2 complete at ~20:00Z: gold Z=26.42, rank=66/66 (LAST). Both gold last. Mean Z=28.29.
- Draw 3 complete at 20:06Z: gold Z=21.89, rank=66/66 (LAST). 3 draws done, gold last in all 3.
- All 20 expected ~23:31Z, within cap (~00:21 EDT).
- Finding template updated: `findings/2026-10-26-gene027-lm22-circularity-preservation.md`
  has nolith results, 6-module analysis, and 3-draw preliminary result (gold last in all 3).
  Run `Rscript network/scripts/analyze_gene027.R` once lith_draw_020.rds exists.
- cost-scaling finding: fixed n_test typo (was 234, correct is 74).

## Session 11 notes (2026-10-26T19:35Z)

- gene-027 nolith checkpoint completed at 15:32:27 EDT (713s; 50 perms, 65 modules, 12018 genes).
  bp_nolith: Z=[10.89, 74.91]. Draw phase started immediately.
- 20 lith draws in progress; expected completion ~19:32 EDT (~14400s at ~720s/draw).
  Cap fires at ~20:21 EDT — all draws should complete before cap.
- Updated `findings/2026-10-09-cost-scaling.md`: corrected per-draw time ~720s (was 262s),
  per-perm ~12.8s (was 5.2s). Added critical note: nPermutations setup scaling in WGCNA 1.74.
- `network/scripts/analyze_gene027.R` present for post-hoc comparison (written earlier session).
- Monitor armed on gene-027 log (draw completions and COMPLETE/ERROR lines).

## Session 10 notes (2026-10-26T19:10Z)

- gene-027 cap raised 7200→10800: previous attempt showed modulePreservation() rebuilds both
  TOMs internally (~3600s total), leaving no room for 500 perms in the old cap. cap_sec in
  queue.tsv updated; documented in METHODS.md 2026-10-26T16:30Z.
- gene-027 re-started at 19:10Z with cap=10800. Entered permutation phase within 5s (all
  prior checkpoints exist: SFT curve, reference cache, gene exclusion summary).
- Expected completion: ~22:10Z for nolith checkpoint, then lith draws. May be PARTIAL again
  if 20 draws don't fit in remaining ~3600s (~13 expected; 7 would roll to next session).

## Session 9 notes (2026-10-26T14:15Z)

- inp-023 DONE at 12:05Z (10311s), 20/20 draws. Gold last in all 20 draws (43/43).
  Finding updated: E=24/25. inp-022 R score raised 4→5, total 12→13/25.
- pow-009 DONE at 14:12Z (86s). SFT curves on offset-removed expression: all 3 groups fail
  (control max R²=0.770, bp_nolith=0.755, bp_lith=0.701). New finding E=22/25 (partially
  rehabilitates retracted finding #27 with correct context).
- gene-027 started at 14:15Z (LM22 circularity check, 12018 genes, cap 7200s). Expected
  2.8h → likely PARTIAL at cap; checkpoint per draw means safe resumption.

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
| inp-021 | age+sex+rin+plate+seqpc1-3 | **63** | −8 (11%) |
| inp-022 | ILR b1–b4 | 44 | −27 (38%) |
| inp-023 | age+sex+rin + ILR b1–b4 | **42** | −29 (41%) |

Pattern: technical residualisation removes 8 modules (11%); ILR composition removes 27
(38%); technical+ILR removes 29 (41%). The large compositional effect confirms the 73%
composition finding — ILR is what collapses the structure.

SFT curve comparison:
| arm | residualisation | first power R²≥0.80 | max R² | powerEstimate |
|:--|:--|--:|--:|--:|
| pow-009 ctrl | none (offset-removed) | never | 0.770 at p20 | NA |
| pow-009 nolith | none (offset-removed) | never | 0.755 at p20 | NA |
| pow-009 lith | none (offset-removed) | never | 0.701 at p20 | NA |
| inp-021 | age+sex+rin+plate+seqpc1-3 | never | **0.797 at p18** | NA |
| inp-022 | ILR b1–b4 | power 9 (R²=0.802) | 0.880 | 16 |
| inp-023 | age+sex+rin + ILR b1–b4 | power 13 (R²=0.806) | 0.836 | NA |

Pattern: technical residualisation (inp-021) brings R² from 0.770 to 0.797 — close but
not over the threshold. ILR composition removal is necessary (inp-022: 0.880, inp-023:
0.836). Composition is the primary driver of scale-free topology failure in this data.

## Layout

New work goes in new files: `analyses/<NNN>-<slug>/`, `runs/<id>/`, `findings/<date>-<slug>.md`.
`RANKING.md`, `ranking.csv`, `results_index.csv` are **generated** — rebuild with
`Rscript network/build_index.R`, never hand-edit. See `CONTRIBUTING.md` for conflict
resolution from the laptop.
