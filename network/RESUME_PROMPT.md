You are resuming the unattended network-analysis arm (Arm 3) on a dedicated iMac.
No human will answer you. Never ask a question — record blockers and work around them.

## Read first, in this order
1. `IMAC_README.md` — orientation, verified facts (§5), traps (§6), hard rules (§8).
2. `network/STATE.md` — what is done, running, queued, blocked.
3. `tail -30 network/RUNLOG/progress.jsonl` — where the last session stopped.
4. `network/METHODS.md` — decisions already made. Do not re-litigate them.

## Scope — absolute
- **Write ONLY inside `network/`.** Everything else in this repository is read-only to
  you. Read anything; change nothing outside. (`IMAC_README.md` §8 calls the repo root
  `Implementation/`; no such directory exists — the root is this checkout.)
- **Stage with `git add network/` — NEVER `git add -A`.** `run_next.sh` enforces this
  and unstages outsiders, but do not rely on the guard.
- **This repository is PUBLIC.** Everything you commit is world-readable within minutes.
  Never write a credential, token or private path into any file under `network/`.
- `python index/build_index.py --query "..."` is **allowed** — it is read-only and is the
  discovery path `IMAC_README.md` §7 prescribes first. Do **not** run it without
  `--query`: a bare run rewrites `INDEX.db`/`INDEX.csv` at the repository root.
- Add no new R packages. The approved set is installed (WGCNA, igraph, impute,
  preprocessCore, GO.db, AnnotationDbi, dynamicTreeCut, fastcluster, limma, edgeR).
  If you need another, record it in `network/BLOCKED.md` and work without it.

## The work loop
The driver will NOT run an experiment whose script does not exist — it releases the row
and exits 2 so you can write it. Nothing is ever retired because you had not got to it.

```bash
bash network/run_next.sh --peek     # next row: id, family, cap_sec, spec, script, exists
# if exists=NO, write network/experiments/<id>.R for that spec, then:
bash network/run_next.sh            # claims, runs under its cap, records, commits, pushes
```

Exit codes: `0` ran, `2` needs a script (write it, then re-run), `3` queue empty.

Repeat until roughly 80% of the usage window is spent. Write
`network/config/<id>.json` — every parameter, the seed, the git SHA, package versions —
**before** computing anything.

`bmind-047` needs `bmind-046` finished first. Its script must check for the bMIND
profiles and stop with a clear message if they are absent; the row is then requeued
automatically rather than producing a wrong result.

## Non-negotiable design points (reasoning in IMAC_README.md §6)
- Define modules ONCE on a fixed reference group; measure those same modules everywhere
  else. Re-detecting per draw makes pooling meaningless. (§6.2)
- Subsample larger groups to n=74. `bp_nolith` is never subsampled — it already is 74.
- Draws are NOT independent (two draws of 74 from 152 share ~36 people). No Fisher, no
  Stouffer, no averaging p-values. (§6.4)
- `Zsummary`/`medianRank` are descriptive — no calibrated null, no multiple-testing
  correction, never pushed through a normal CDF. (§6.3)
- P-values come only from a group-label permutation null. Correct **across modules
  within a comparison** using BH. When you aggregate across several comparisons or
  families, do **not** pool them into one BH — use **Benjamini–Bogomolov**, and reuse the
  validated `de_v2/scripts/mtc.R` rather than writing your own or adding `qvalue`/`IHW`,
  which are deliberately absent. (§6.6, §6.8)
- **Ensembl ids in `data/cohort_474/counts_474.tsv.gz` carry version suffixes**
  (`ENSG00000000003.10`); everything else in the repo is unversioned. Strip with
  `strsplit(x, ".", fixed = TRUE)` — a regex needs `"\\."`. Getting this wrong silently
  produced empty intersections for 45 minutes once. `00_load.R` already returns
  unversioned ids, so prefer it. (§6.5)
- `enableWGCNAThreads()` or WGCNA runs single-threaded. (§6.9)
- **Always pass `permutedStatisticsFile` to `modulePreservation()`.** It defaults to
  `savePermutedStatistics = TRUE` writing `permutedStats-actualModules.RData` into the
  working directory — the repository root, outside your scope, under one name shared by
  every experiment. Two runs would overwrite each other's null, and
  `loadPermutedStatistics = TRUE` could read a different experiment's permutations and
  report a wrong p-value. Always:
  `permutedStatisticsFile = file.path("network/runs", id, "permutedStats.RData")`
- **NEVER trust `pickSoftThreshold$powerEstimate` on this data.** Measured at the full
  12,368 genes, n=74: the scale-free fit never reaches R²≥0.80 at any power in 1:20, and
  R² is *negative* at low powers (control −0.97, −0.97, −0.93 at powers 1–3).
  `powerEstimate` returns **1** for control and bp_nolith, **2** for bp_lith — the
  function failing, not a threshold. A network at power 1 is an unthresholded correlation
  matrix and collapses to ~2 modules. **Use power 14** (WGCNA's default for signed,
  n>40) unless a spec fixes another on purpose. Guard in code: reject any chosen power
  < 4, fall back to 14, and log that it happened. Record the full R² curve every time —
  per catalogue item 4 the curves are themselves a result, and the failure of scale-free
  topology in this data is a finding worth reporting.
- Pre-specify every threshold in the config before computing the result.

## Checkpointing
Anything over ~15 minutes writes partial results as it goes — one `.rds` or `.csv` per
draw, per permutation block, per group — and skips work whose checkpoint already exists.
A resumed experiment must not restart from zero.

## Spend the last slice writing
Update `network/RESULTS.md`, `network/results_index.csv`, `network/STATE.md`, and append
to `network/METHODS.md`. Commit and push. Then top the queue back up so it never empties;
if the catalogue is exhausted, invent defensible experiments and record why.

## Reporting
Report what you find, including nothing. A well-documented negative result is a complete
deliverable. Never fabricate a number, never silently skip a check, never report a
partial run as complete. If a result contradicts a verified fact in `IMAC_README.md` §5,
stop that analysis and write the discrepancy down in full rather than adjusting the
pipeline to match.
