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

Exit codes: `0` ran, `2` needs a script (write it, then re-run), `3` nothing claimable.

`--peek` is strictly read-only — it never changes the queue, so it is safe to call as
often as you like. It therefore does NOT requeue: if it reports `QUEUE_EMPTY_PENDING`
while showing `PARTIAL` or `RUNNING` counts, run `bash network/run_next.sh` once anyway —
the driver requeues stale `RUNNING` and retryable `PARTIAL` rows before it claims.

Repeat until roughly 80% of the usage window is spent. Write
`network/config/<id>.json` — every parameter, the seed, the git SHA, package versions —
**before** computing anything.

`bmind-047` needs `bmind-046` finished first. Its script must check for the bMIND
profiles and stop with a clear message if they are absent; the row is then requeued
automatically rather than producing a wrong result.

## Budget reality — measured, not estimated

All figures below are **observed on this machine**, not extrapolated:

| quantity | measured | source |
|:--|:--|:--|
| reference network build | 786 s x (n_ref / 74) | base-001: 47 min at n=234 |
| `modulePreservation` | **5.2 s per permutation** | base-001: 262 s per 50-perm draw |
| one subsample draw (50 perms) | **262 s (4.4 min)**, sd ~3 s | base-001, 9 draws, range 261–269 s |

**There is NO per-draw reference rebuild.** `modulePreservation` takes the reference data
once; only the n=74 test network is built per draw. An earlier estimate of 1962 s/draw
assumed otherwise and was **7.5x too high** — it is what drove `base-001` down from 100
draws to 20. Do not repeat that reasoning.

Cost of one `draws=D;perms_per_draw=50` row with an n=234 reference:
`2360 s (build, paid once per cache key) + 2620 s (nolith, 500 perms) + D x 262 s`

- D=20 -> **~2.8 h**  (all 11 such rows: ~31 h, comfortably inside the ~61 h available)
- D=100 -> **~8.7 h** (all 11: ~89 h, roughly 1.5x oversubscribed)

So the real tradeoff is **breadth vs. precision of the across-draw spread**, and it is a
scientific choice, not a budget constraint:

- **20 draws** fits every row, and gives a spread good enough to see whether a shift
  exceeds draw-to-draw noise (base-001 at 5 draws already separated turquoise from the rest).
- **100 draws** is the catalogue spec and estimates the spread far better, but roughly
  five rows then consume the run.

Reuse of the cached reference across rows sharing a cache key saves the 2360 s build each
time and is what makes either option viable. Record in `config.json` which cache entry a
run used.

If you change a row's draw count, say why in `METHODS.md` **before** computing, and state
which of breadth or precision you traded away.

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
  matrix and collapses to ~2 modules. **Use power 12**, the value WGCNA's FAQ recommends for a *signed*
  network with n>40 when no power satisfies the scale-free criterion (the FAQ table is
  signed: <20->18, 20-30->16, 30-40->14, >40->12; unsigned/hybrid: 9/8/7/6). An earlier
  draft of this file said 14, which is the 30-40 row - wrong for n=74. Unless a spec fixes
  another power on purpose, use 12. Guard in code: reject any chosen power < 4, fall back
  to 12, and log that it happened. Record the full R² curve every time —
  per catalogue item 4 the curves are themselves a result, and the failure of scale-free
  topology in this data is a finding worth reporting.
- Pre-specify every threshold in the config before computing the result.

## Checkpointing
Anything over ~15 minutes writes partial results as it goes — one `.rds` or `.csv` per
draw, per permutation block, per group — and skips work whose checkpoint already exists.
A resumed experiment must not restart from zero.

## Rank every result, every time — this is not optional

After **each** experiment, update `network/RANKING.md` and `network/ranking.csv`. There
are **two** rankings and you maintain both:

**Table A — Evidence (0–25).** Five components 0–5: Evidence, Relevance, Robustness,
Decisiveness, Completeness. The rubric is in RANKING.md — read it before scoring. It is
built so a well-powered **null scores high** and a striking result resting on one draw
scores low, because this arm is exploratory and `IMAC_README.md` forbids making something
look significant.

**Table B — Interest (0–10).** How striking, surprising or quotable the finding is. This
is the honest answer to "did we find anything cool" and it is worth recording.

Also record **Δ = (Interest/10 − Evidence/25) × 10**. The gap is the point: Interest well
above Evidence means *be suspicious and go get more evidence before anyone quotes it*;
Evidence well above Interest means *solid, unexciting, still worth reporting*.

**Interest never decides what runs next.** Queue order is scientific priority, set above.
Do not reorder the queue, extend an analysis, or choose a variant because it scored high
on Interest. Rank it, note it, move on. If a high-Interest result tempts you to chase it,
that is precisely the moment to write the Δ down and continue with the queue.

Rules for both tables: score what the run demonstrated, not what it suggests; re-rank
against everything already there so the tables stay live; never quietly re-score an old
result to make a new one look consistent — if a score changes, say so in Notes and in
`METHODS.md`; cap `PARTIAL` rows at C=2 and rank them on what they actually produced.

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
