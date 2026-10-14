# Arm 3 — Run State

**Updated:** 2026-10-14

## Where the science stands

**The arm's founding question has a robust provisional answer.** A lithium-associated
co-expression dimension exists (PC2 of module eigengenes, permutation p=0.0015), is 73%
cell composition (73.3% at 33 modules, BH survives across all comparisons), replicates
across all 10 usable decompositions in the deepSplit × mergeCutHeight grid, and **does not
survive composition adjustment in any decomposition**. The monotone relationship between
decomposition quality and composition share (r=+0.879) suggests a real effect being
resolved more sharply as the network improves, not a parameter artifact.

Key: `findings/2026-10-14-lithium-dimension-is-composition.md` (E=24/25, I=9/10),
`findings/2026-10-14-lithium-dimension-replicates.md` (E=22/25).

The next question is whether module structure is preserved across groups at all in these
data (`base-001` v2 / `base-003-remerge`), and whether residualising for composition changes
that (`inp-022`, `inp-023`).

**Read `RANKING.md` first.** Rebuilt from `findings/`; 8 findings as of 2026-10-14.

## Infrastructure

- Scheduler: **LaunchDaemon** `com.aaylab.networkarm.daemon` (system domain, survives
  logout; auth verified with a stripped environment). The earlier LaunchAgent died at the
  login window and cost 4 days 18 hours — see `METHODS.md` 2026-10-14.
- `com.aaylab.networkarm` (agent) retained as redundancy; `wake.sh`'s lock prevents overlap.
- Hourly watchdog `com.aaylab.networkarm.watchdog` kickstarts the agent if no wake is
  recorded for 6 h. Note it is a LaunchAgent and is therefore **also** dead when logged out;
  the daemon is the real protection.

## Queue

`mod-000-resolution` DONE. `base-001` **RUNNING v2** (deepSplit=4, mergeCutHeight=0,
33-module reference from mod-000-resolution labels fast-path; nolith 500-perm preservation
started ~19:15Z 2026-10-14; 20 lith draws follow). ETA completion ~21:30Z.

Order: **`base-001`** (RUNNING, 33-module v2) → **`base-003-remerge`** (reads base-001
results, ~2 min) → `base-002-ceil` → `inp-022` → `inp-023`.

All 50 PENDING experiments have scripts. All scripts use `deepSplit=4, mergeCutHeight=0`
(builder scripts) or point to the v2 cache key `control_all_bicor_signed_p12_ds4_mch00`
(reader scripts). Fixed 2026-10-14.

## Parameters every downstream row should use

`power = 12`, `deepSplit = 4`, `mergeCutHeight = 0` → 33 modules, largest 11.68% of genes.
**Never `mergeCutHeight = 0.25`**: it collapses 33 modules to 3 with one holding 97.5% of
the transcriptome, because the median eigengene pair correlates 0.739 against a 0.75
threshold.

Quote every Zsummary against WGCNA's `gold` random-gene module, never against Langfelder's
published bands — gold scores 20.26 here, so the bands place noise in "strongly preserved".

## Analyses done (session 5, 2026-10-14)

- **001-module-trait-ds4**: No module eigengene associates with lithium (0/33, BH). Modules
  not independent (median |r|=0.822).
- **002-effective-dimensionality**: PC2 associates with lithium (β=-0.908, p=0.0015);
  73.3% cell composition; composition attenuates by 68.1%.
- **003-lm22-circularity**: Dropping all 350 LM22 markers: association and composition
  share move <1%. Not circular.
- **004-decomposition-robustness**: Lithium dimension (PC2) replicates in 6/10 decomposition
  cells, including 4/4 usable (≥10 modules, largest <25%). 0/6 survive composition
  adjustment. r(n_modules, composition share)=+0.879.

## Layout

New work goes in new files: `analyses/<NNN>-<slug>/`, `runs/<id>/`, `findings/<date>-<slug>.md`.
`RANKING.md`, `ranking.csv`, `results_index.csv` are **generated** — rebuild with
`Rscript network/build_index.R`, never hand-edit. See `CONTRIBUTING.md` for conflict
resolution from the laptop.
