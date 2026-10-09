# Arm 3 — Run State

**Updated:** 2026-10-09T02:50:31Z
**Status:** SETUP COMPLETE, VERIFIED — scheduler NOT yet loaded, no experiment has run

## Environment (all verified, not assumed)
- R 4.2.2, 4 cores, 16 GB RAM, ~189 GB free.
- All ten packages installed and loading: WGCNA 1.74, igraph 2.3.4, impute 1.72.3,
  preprocessCore 1.60.2, GO.db 3.16.0, AnnotationDbi 1.60.2, dynamicTreeCut 1.63.1,
  fastcluster 1.3.0, limma 3.54.2, edgeR 3.40.2. `enableWGCNAThreads()` → 3 workers.
- `Rscript network/00_load.R` prints the required line exactly:
  `12368 genes x 474 samples | control 234, bp_nolith 74, bp_lith 152`
- Headless CLI authenticates (`claude -p` → CLI_OK) using a one-year token in
  `~/.arm3_token` (0600, outside this public repository).
- `bash network/wake.sh --check` passes with no side effects.

## Measured costs (not estimated — see METHODS.md)
| stage | exponent in n_genes | at 12,368 genes |
|:--|:--|:--|
| TOMsimilarity | 2.88 | ~12 min |
| modulePreservation, per permutation | 0.80 | ~8 s |

K=500 permutations ≈ 1.1 h at the full gene set. Peak RSS 2.05 GB at 6,000 genes;
~5–6 GB projected at full size. **Never run two experiments concurrently.**

## Known scientific constraint, decided before any run
Scale-free topology **fails in this data**: R² never reaches 0.80 at any power in 1:20,
in any group, and is negative at low powers. `powerEstimate` returns 1 (control,
bp_nolith) or 2 (bp_lith) — the function failing, not a threshold.
**Use power 14.** Record the R² curve every time; the failure is itself a result.

## Machinery
- `network/run_next.sh` — driver. `--peek` inspects the next row without claiming.
  Exit 0 ran, 2 needs a script, 3 queue empty. No status is absorbing except BLOCKED,
  and BLOCKED needs MAX_ATTEMPTS genuine failures.
- `network/wake.sh` — session launcher. 4 h cap, token validation, auth-failure alert,
  zero-throughput watchdog. `--check` for side-effect-free testing.
- `network/RESUME_PROMPT.md` — standing instruction.
- Scheduler: launchd `com.aaylab.networkarm`, 18000 s (5 h), RunAtLoad true.
  **NOT LOADED.** Loading it starts the run.

## Live ranking
`network/RANKING.md` (human) and `network/ranking.csv` (parsable) hold a running
ranking of every result, rebuilt after each experiment. The rubric rewards evidential
strength, not how striking a finding looks — a well-powered null can score top.

## Done
(no experiment has run)

## Next
1. `base-001` — the baseline: ref=control, bicor, signed, all genes, 100 draws,
   descriptive Zsummary + medianRank. **Nothing else is trusted until this completes.**
2. `inp-022` / `inp-023` — composition-residualised. These are the arm's real question.
3. Breadth across the remaining 45 rows.

The session must write `network/experiments/<id>.R` before the driver will run a row.
Use `--peek` to see the spec.

## Blocked
See `network/BLOCKED.md`. Nothing outstanding that prevents the run from starting.

## Scope
Write only inside `network/`. Stage with `git add network/`, never `-A`.
**This repository is PUBLIC** — never commit a credential or private path.
`build_index.py --query` is allowed (read-only); a bare run is not.
