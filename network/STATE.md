# Arm 3 — Run State

**Updated:** 2026-10-09T02:12:50Z
**Git SHA at setup:** a931461
**Status:** SETUP COMPLETE — not yet started

## What exists
- Queue seeded: 48 experiments, 20 families, all `PENDING`.
- Driver: `network/run_next.sh` (claims, runs, caps, marks, commits, pushes).
- Scheduler: launchd `com.aaylab.networkarm`, 5 h interval, **not yet loaded**.
- Standing prompt: `network/RESUME_PROMPT.md`.
- Data verified: `12368 genes x 474 samples | control 234, bp_nolith 74, bp_lith 152`.

## Done
(nothing — no experiment has run)

## Running
(none)

## Next
1. `base-001` — the baseline. ref=control, bicor, signed, all genes, 100 draws,
   descriptive Zsummary + medianRank. **Must finish before any variant is trusted.**
2. Then `inp-022` / `inp-023` (composition-residualised) — these are the arm's
   actual question.
3. Then breadth across the remaining families.

## Blocked
See `network/BLOCKED.md`. One item outstanding at setup: the headless CLI
needs `claude login` before any scheduled session can do anything.

## Scope reminder
Write only inside `network/`. Stage with `git add network/`, never `-A`.
Do not run `index/build_index.py` — it writes at the repository root.
