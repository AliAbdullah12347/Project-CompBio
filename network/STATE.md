# Arm 3 — Run State

**Updated:** 2026-10-14

## Where the science stands

**The arm's question has a provisional answer.** At the module-eigengene level in BP1
(n=226): a lithium-associated dimension exists in the co-expression structure (PC2,
permutation p = 0.0015) and it is **73% cell composition**; adjusting for composition
attenuates the association by 68% to p = 0.081. See
`findings/2026-10-14-lithium-dimension-is-composition.md`.

That is the project's founding hypothesis measured on the network side. It is one
decomposition and one sample — `base-003-remerge` asks the preservation version of the
same question and is queued next.

**Read `RANKING.md` first.** It is generated from `findings/`; six findings so far.

## Infrastructure

- Scheduler: **LaunchDaemon** `com.aaylab.networkarm.daemon` (system domain, survives
  logout; auth verified with a stripped environment). The earlier LaunchAgent died at the
  login window and cost 4 days 18 hours — see `METHODS.md` 2026-10-14.
- `com.aaylab.networkarm` (agent) retained as redundancy; `wake.sh`'s lock prevents overlap.
- Hourly watchdog `com.aaylab.networkarm.watchdog` kickstarts the agent if no wake is
  recorded for 6 h. Note it is a LaunchAgent and is therefore **also** dead when logged out;
  the daemon is the real protection.

## Queue

`mod-000-resolution` DONE. `base-001` PARTIAL at 18/20 draws, on the degenerate
decomposition — its numbers are valid only for that parameter set.

Order: `base-001` (finish) → **`base-003-remerge`** (the key re-run at deepSplit=4,
mergeCutHeight=0) → `base-002-ceil` → `inp-022` → `inp-023`.

## Parameters every downstream row should use

`power = 12`, `deepSplit = 4`, `mergeCutHeight = 0` → 33 modules, largest 11.68% of genes.
**Never `mergeCutHeight = 0.25`**: it collapses 33 modules to 3 with one holding 97.5% of
the transcriptome, because the median eigengene pair correlates 0.739 against a 0.75
threshold.

Quote every Zsummary against WGCNA's `gold` random-gene module, never against Langfelder's
published bands — gold scores 20.26 here, so the bands place noise in "strongly preserved".

## Layout

New work goes in new files: `analyses/<NNN>-<slug>/`, `runs/<id>/`, `findings/<date>-<slug>.md`.
`RANKING.md`, `ranking.csv`, `results_index.csv` are **generated** — rebuild with
`Rscript network/build_index.R`, never hand-edit. See `CONTRIBUTING.md` for conflict
resolution from the laptop.
