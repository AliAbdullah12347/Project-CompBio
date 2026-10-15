# Defect: rows whose cap equals SESSION_CAP can be BLOCKED while making progress

**Found:** 2026-10-15, before the unattended multi-day stretch.
**Affected rows:** `perm-031`, `perm-032`, `bmind-046` (all had `cap_sec = 14400`).
**Severity:** would have retired three experiments as BLOCKED despite every attempt
producing valid checkpoints. `perm-031` is queue position 8, so it would have hit
within ~2 days.

## Mechanism

`run_next.sh:74-77` justifies the BLOCKED threshold like this:

> The driver refunds attempts when a run produced real checkpoints, so by the time
> `$4` reaches maxa here the row has genuinely failed that many times without
> producing anything.

That invariant holds **only when `run_next.sh` exits normally.** The refund is computed
in the normal exit path. But `wake.sh:118-123` enforces the session cap by killing the
whole process tree:

```sh
pkill -TERM -P $CPID; kill -TERM $CPID
sleep 15; pkill -KILL -P $CPID; kill -KILL $CPID
```

`run_next.sh` is a descendant of `$CPID`, so it dies there. Its only trap is
`trap 'rm -rf "$LOCKDIR"' EXIT` — the lock is released, but the row is never moved to
PARTIAL and the refund never runs. The row stays RUNNING, `requeue()` later sees it as
stale, and increments toward `MAX_ATTEMPTS` **with no refund**. Three cycles -> BLOCKED.

For most rows this is rare: the row's own cap fires first, the driver exits cleanly, and
the refund applies. But when `cap_sec == SESSION_CAP`, the session kill is *guaranteed*
to win — the session clock starts a minute or two before the row is even claimed. Those
rows can never exit cleanly, so they can never earn a refund.

## Fix applied

Lowered `cap_sec` 14400 -> 10800 for the three affected rows, leaving ~1 h of headroom
inside the 4 h session. The row's own cap now fires first, the driver records PARTIAL,
the progress refund applies, and results are committed before the session ends. All three
scripts are checkpointed, so a PARTIAL simply resumes next session.

Data-only change (`network/RUNLOG/queue.tsv`). `run_next.sh` was deliberately **not**
edited: pid 67870 was executing it at the time, and bash reads scripts incrementally, so
editing a running script risks corrupting the live execution.

## Better fix, for when the driver is idle

Give `run_next.sh` a `trap ... TERM INT` that marks the claimed row PARTIAL and applies
the refund before exiting. That removes the invariant's dependence on a clean exit and
protects every row, not just the three with the pathological cap. Until then, keep the
rule: **`cap_sec` must be at least ~1 h below `SESSION_CAP`.**
