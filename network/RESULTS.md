# Arm 3 — Results

**Updated:** 2026-10-14
**Runs completed:** 1 of 50 queued — the run was halted by a logout, not by a failure.

## Status

The run executed two sessions on 2026-10-09/10 and then stopped for 4 days 18 hours
because the Mac was logged out. GUI LaunchAgents do not run at the login window, so both
the scheduler and its watchdog were dead. See `METHODS.md` 2026-10-14. Nothing was lost:
all checkpoints, logs and results are intact and committed.

## Completed

| id | wall | verdict | folder |
|:--|--:|:--|:--|
| `mod-000-resolution` | 2938 s | **Arm unblocked.** Merge threshold is decisive; `deepSplit=4, mergeCutHeight=0` gives 33 modules, largest 11.68% | [network/runs/mod-000-resolution](runs/mod-000-resolution) |

## Partial

| id | progress | note |
|:--|:--|:--|
| `base-001` | 18/20 draws, nolith complete (500 perms) | Computed on the **degenerate** `mergeCutHeight=0.25` decomposition (3–4 modules, largest 95–97%). Valid as a record of that parameter set; **must be re-run** at the recommended parameters before anything is concluded. |

## Headline findings so far

All five findings to date are about **method, not biology**, and they form one chain:

1. **Scale-free topology fails** at every power 1–20, in all three groups, at both n=74 and
   n=234. `pickSoftThreshold` returns `powerEstimate = 1` — the function failing, not a
   threshold. Max signed R² is 0.547 at power 20.
2. **The network is extraordinarily dense**: mean connectivity 1,031 at power 12 — each
   gene tied to 8.3% of the genome.
3. **Module eigengenes are therefore mostly correlated** (median pair 0.739–0.765), so
   WGCNA's default `mergeCutHeight = 0.25` (merge above 0.75) cascades and collapses
   33 modules into 3, the largest holding 97.5% of genes.
4. On that degenerate decomposition, **no module reliably beats a random gene set**:
   WGCNA's own `gold` control scores Zsummary 20.26, and the best real module beats it in
   only 9 of 15 subsamples.
5. **A usable decomposition exists**: `deepSplit=4, mergeCutHeight=0` gives 33 modules,
   largest 11.68%, median size 248.

**Nothing here yet speaks to lithium.** Finding 4 is a property of the parameter set in
finding 3, and finding 5 is the fix. The next experiment is the one that matters.

## Next

1. **Re-run `base-001` at `deepSplit=4, mergeCutHeight=0`.** This tests whether the gold
   result survives a usable decomposition — the difference between "WGCNA preservation is
   uninformative on this data" and "that parameter set was".
2. `base-002-ceil` — the split-half control ceiling.
3. `inp-022` / `inp-023` — the arm's actual question, composition-residualised.

## The question

Do gene co-expression networks differ across healthy / bipolar-off-lithium /
bipolar-on-lithium, and does any difference survive removing cell composition?

Exploratory. A well-documented negative or inconclusive result is a complete deliverable.
