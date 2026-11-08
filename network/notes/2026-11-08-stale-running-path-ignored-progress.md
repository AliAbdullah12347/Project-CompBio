# Three experiments were retired while holding 57 checkpoints between them

`requeue()` retires a stale RUNNING row at `MAX_ATTEMPTS` with
`"dead RUNNING after N attempts, no progress"` — **without looking at its output.**

The progress refund (`run_next.sh:355-361`: non-empty, non-stray files newer than a
marker) runs only on the clean-exit PARTIAL path. A row whose driver is hard-killed never
reaches it: it stays RUNNING, goes stale, and is blocked on the attempt counter alone.

Found on 2026-11-08. All three BLOCKED rows had substantial work banked:

| row | .rds | .csv | how far it got |
|:--|--:|--:|:--|
| `ref-002` | 30 | 2 | ctrl draw **10 of 20** |
| `ref-004` | 18 | 1 | — |
| `corr-006` | 9 | 2 | — |

This is the same family as the 2026-10-25 fix, which made a dead claimer free its row
*faster*. That fix addressed latency; it did not add a progress check. The stale path
still discarded the evidence.

## Fix

`rows_with_work()` lists run directories holding at least one non-empty, non-stray
`.rds`/`.csv` — the same test the clean path uses. `requeue()` now refunds an attempt and
returns such a row to PENDING instead of retiring it. A row with nothing but zero-byte
files still blocks, so the absorbing state remains reachable.

Sandbox-verified both branches: a row with one real checkpoint → PENDING with the
attempt refunded; a row with only a zero-byte file → BLOCKED.

Applied with write-temp-then-rename, because a driver was executing the script.

## Also

`corr-006`'s cached reference was built under `corType = "spearman"`, which WGCNA does
not support — that is why the row kept dying. The local fix rank-transforms and uses
`cor` instead. The stale cache is quarantined so it rebuilds under the corrected code.
