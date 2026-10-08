# MASTER PROMPT — unattended network-analysis run

Paste the block below into Claude Code on the iMac, in the `Implementation/`
checkout, as the **first** message of a fresh session. Everything after the
horizontal rule is the prompt.

Before you paste it, make sure the Mac will stay awake and signed in to GitHub:

```bash
caffeinate -dimsu &
git -C ~/Project-CompBio push --dry-run origin main   # must succeed without prompting
```

If that push asks for a password, stop and set up a credential helper or SSH key
first. The run cannot recover from an auth prompt with nobody there.

---

You are running unattended on a dedicated iMac for several days. **No human will
see your output or answer a question after the first twenty minutes.** Plan for
that from the first line.

**Read `IMAC_README.md` in the repository root before anything else.** It tells you
where the data is, what is deliberately missing and how to regenerate it, the
verified facts your work must reproduce, the traps that have already cost this
project time, and what to do when you cannot find something. Everything below
assumes you have read it.

Your job is the **network-analysis arm (Arm 3)**: do gene co-expression networks
differ across healthy / bipolar-off-lithium / bipolar-on-lithium, and does any
difference survive removing cell composition? The arm is **exploratory**. A
well-documented negative or inconclusive result is a complete deliverable. Making
something look significant is not.

## Phase 0 — the only time you may ask for anything (first 20 minutes)

Do these in order, then ask **one** consolidated approval message containing every
permission and install you will need for the whole run. After that message, never
ask again — if something later needs approval you did not request, record it in
`network/BLOCKED.md` and work around it.

1. Confirm the repo loads: `Rscript network/00_load.R` must print
   `12368 genes x 474 samples | control 234, bp_nolith 74, bp_lith 152`.
   If it does not, stop and report — nothing else is worth doing.
2. Inventory the environment: R version, installed packages, CPU cores, free disk,
   free RAM. Check at minimum `WGCNA impute preprocessCore GO.db AnnotationDbi
   dynamicTreeCut fastcluster limma edgeR igraph`.
3. Time one small experiment end to end (one group, 2,000 genes, one preservation
   run) and record the wall clock. You need this number to size everything else.
4. Verify `git push` works without prompting.
5. **Then** send one message listing: packages to install, anything needing disk or
   network access, and your proposed per-experiment time cap. Wait for the reply.

## Phase 1 — build the machinery before you build any science

A single long session will die — to a crash, a usage limit, a reboot, or a power
cut. Do not rely on one. Build this instead:

**A queue.** `network/RUNLOG/queue.tsv`, tab-separated, columns:
`id  family  status  attempts  started  finished  wall_sec  note`
Status is one of `PENDING RUNNING DONE PARTIAL BLOCKED`. Seed it with at least 40
experiments from the catalogue below before you run any of them, so a crashed
session always has something queued to pick up.

**A driver.** `network/run_next.sh` — claims the oldest `PENDING` row, marks it
`RUNNING` with a timestamp, runs it, marks the outcome, commits, pushes. At start
it resets any `RUNNING` row older than 6 hours back to `PENDING` and increments
`attempts`; at `attempts >= 3` it marks `BLOCKED` with the reason and moves on.
A row must never be able to wedge the queue.

**A schedule.** A cron entry firing every five hours so each invocation starts in a
fresh usage window:

```
0 */5 * * * cd $HOME/Project-CompBio && /usr/local/bin/claude -p "$(cat network/RESUME_PROMPT.md)" --dangerously-skip-permissions >> network/RUNLOG/cron.log 2>&1
```

Write `network/RESUME_PROMPT.md` yourself: a short standing instruction telling a
fresh session to read `IMAC_README.md`, read `network/STATE.md`, and work the queue
until the window is nearly spent. `--dangerously-skip-permissions` is appropriate
here precisely because this is a dedicated machine doing one job; confirm it in your
Phase 0 message.

**A heartbeat.** Append one JSON line to `network/RUNLOG/progress.jsonl` at every
step boundary — never less often than every 10 minutes of compute — with
`{ts, run_id, phase, pct, msg}`. Long R loops write a line per iteration. If the
machine dies, the last line must tell the next session exactly where it was.

**Checkpointing inside experiments.** Anything over ~15 minutes writes partial
results as it goes (one `.rds` or `.csv` per draw, per permutation block, per
group) and skips work whose checkpoint already exists. A resumed experiment must
not restart from zero. This project has already lost 45 minutes of fits to a
missing checkpoint — do not repeat it.

**A time cap.** Every experiment gets a wall-clock cap. On exceeding it: write what
you have, mark `PARTIAL` with how far it got, move to the next. Never let one
experiment consume a whole window.

**Push constantly.** After every finished experiment, and at least hourly
regardless: `git pull --rebase` then `git add -A && git commit && git push`. If the
push fails, retry three times with backoff, then log it and keep working — never
let a network failure stop the science. Each commit message says what finished and
what it found, in one line each.

## Phase 2 — the science

Run as many distinct analyses as you can. Breadth is the goal; prefer ten honest
variations over one perfect one. The catalogue below is a **starting point, not a
limit** — extend it with anything defensible you think of, and record your reasoning
when you add something.

**Core design points that are not negotiable** (details and reasoning in
`IMAC_README.md` §6):

- Define modules **once** on a fixed reference group, then measure those same
  modules everywhere else. Re-detecting modules per draw makes pooling meaningless.
- Subsample the two larger groups to n = 74. BP1-off-lithium is **never**
  subsampled — it already is 74.
- Draws are **not** independent (two draws of 74 from the 152 share ~36 people).
  Never combine them with Fisher, Stouffer, or by averaging p-values.
- `Zsummary` and `medianRank` are **descriptive**. Neither has a calibrated null,
  so neither can be multiple-testing corrected. If you want p-values, build them
  from a **group-label permutation** null, matching the number of draws on the
  observed and permuted sides, and correct across **modules** (BH), not draws.
- Pre-specify every threshold in the run's config **before** computing the result.
  Reporting only the variant that worked would invalidate the arm.

**Catalogue.** Cross these freely; each distinct combination is an experiment.

1. **Reference group**: controls / bp_nolith / bp_lith as the module-defining set;
   also WGCNA consensus modules across all three.
2. **Correlation**: Pearson, Spearman, biweight midcorrelation.
3. **Network type**: signed, unsigned, signed hybrid.
4. **Soft threshold**: scan `powerVector = 1:20`; also fixed 6, 8, 12. Record the
   scale-free R² curve per group — the curves themselves are a result.
5. **Module detection**: `deepSplit` 0–4, `minModuleSize` 20/30/50/100,
   `mergeCutHeight` 0.15/0.25/0.35.
6. **Input matrix**: raw log2TPM; residualised for technical covariates; residualised
   for composition (ILR b1–b4); residualised for both; voom-weighted from counts.
7. **Gene set**: all 12,368; top 2,000/5,000 most variable; the 547 LM22 markers
   **excluded** (circularity check); lithium DEGs only; explicitly non-DEG genes.
8. **Preservation**: Zsummary and medianRank, plus the density and connectivity
   components separately — they can disagree, and that disagreement is informative.
9. **Permutation test**: group-label null as specified above, K ≥ 500.
10. **Module–trait relationships**: eigengene correlation with lithium, diagnosis,
    age, sex, RIN, plate, each of the 5 lineage fractions, each ILR balance.
11. **Hub genes**: intramodular connectivity (kME); hub identity stability across
    the 100 draws; do hubs change between groups?
12. **Differential connectivity**: per-gene connectivity differences between groups,
    with the same permutation null.
13. **Community detection outside WGCNA**: Louvain and Leiden on the adjacency;
    compare partitions to WGCNA modules and to each other (ARI, NMI).
14. **Centrality**: degree, betweenness, eigenvector; between-group differences.
15. **Module annotation**: overlap of each module with the lithium and bipolar DEG
    lists in `de_final/results/`, with LM22 markers, and with the cell-type DE
    results. Fisher exact with BH across modules.
16. **Sample-level structure**: do samples separate by group in eigengene space?
    Outlier detection; hierarchical clustering; does any clustering track lithium?
17. **Gene-filter sensitivity**: reuse the 31 filters from the DE sweep
    (`de_v2/scripts/13_filter_sweep_wb.R`) and ask which modules survive all of them.
18. **Tobacco imputation sensitivity**: repeat a residualised analysis across several
    of the 20 imputations; does any conclusion move?
19. **Cell-type-specific networks** from the bMIND profiles. These need a ~3.1 h
    rebuild (`Rscript de_analysis/scripts/06_bmind_profiles.R`) — give it a window
    of its own, checkpoint it, and only then build networks per lineage.

Prioritise: get one complete, correct baseline through to a written result **first**
(controls as reference, signed, bicor, all genes, 100 draws, descriptive Zsummary).
Then vary. A broad sweep with no finished baseline is worth much less.

## Phase 3 — output that someone can actually read

Both machine-parsable and human-readable, every time.

**Per experiment**, `network/runs/<id>/`:
- `config.json` — every parameter, the seed, the git SHA, the package versions.
- `result.json` — the headline numbers, flat keys, no nesting beyond one level.
- `*.csv` — per-module and per-gene tables. Consistent column names across all
  runs: `run_id, module, n_genes, statistic, value, group, comparison`.
- `README.md` — five to fifteen lines: what was asked, what was done, what came
  out, what it means, what would change the answer. Written for a reader who has
  not seen the code.
- `log.txt` — the full console log.

**Repository level**, rebuilt after every experiment:
- `network/RESULTS.md` — one table of every run: id, family, key parameters,
  headline number, verdict, link to the folder. Newest first. This is the file a
  human opens first, so keep it tight and scannable.
- `network/results_index.csv` — the same thing, parsable, one row per run.
- `network/STATE.md` — what is done, what is running, what is queued, what is
  blocked, and what the next session should pick up. Overwritten each run.
- `network/BLOCKED.md` — append-only. Every gap, failure and workaround, with
  enough detail to act on without this conversation.
- `network/METHODS.md` — append-only decision log. Every methodological choice,
  why you made it, and what you rejected. This is the file the paper is written
  from; treat it as the real deliverable.

**Reporting rules.** Report what you find, including nothing. State uncertainty as
uncertainty. If a result contradicts a verified fact in `IMAC_README.md` §5, stop
that analysis and write the discrepancy down in full rather than adjusting the
pipeline to match. If an experiment fails, the failure and its cause go in
`BLOCKED.md` and the run continues with the next item. Never fabricate a number,
never quietly skip a check, and never report a partial run as complete.

## Phase 4 — standing behaviour for every later session

Each time you wake:

1. Read `IMAC_README.md`, then `network/STATE.md`, then the last 30 lines of
   `network/RUNLOG/progress.jsonl`.
2. Reset stale `RUNNING` rows, as above.
3. Work the queue until roughly 80% of the usage window is spent.
4. Spend the last slice writing: update `RESULTS.md`, `STATE.md`, `METHODS.md`,
   rebuild `python index/build_index.py`, commit, push.
5. Before stopping, **top the queue back up** so it never empties. If you are out
   of catalogue items, invent defensible new ones and write down why they are worth
   running.

Keep going until the queue is exhausted and no new experiment can be justified. If
you genuinely reach that point, write a final summary in `network/RESULTS.md`
saying so, and what you would do with more time.
