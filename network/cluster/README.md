# Arm 3 on the research cluster

Everything left in `network/NEXT_STEPS.md` that is bulk compute, moved off the
iMac. Roughly 113 hours of serial work finishing in about four hours wall clock.

## The three jobs, in order

| job | what | pods × CPU | wall clock |
|:--|:--|:--|:--|
| `arm3-prep` | rebuild the reference cache + compile the R library, publish both to S3 | 1 × 16 | ~30 min |
| `arm3-queue` | drain the 31 PENDING queue rows that have scripts | 8 × 3 | ~3 h |
| `arm3-perm031` | perm-031 from K=20 to K=200 (180 new draws) | 12 × 2 | ~3 h |

`arm3-queue` and `arm3-perm031` both fail fast without the cache, so prep must
finish first. They then run concurrently.

```bash
export KUBECONFIG=./kubeconfig-user-aabdullah.yaml
kubectl create -f network/cluster/arm3-prep-job.yaml
# wait for the cache to appear, then:
kubectl create -f network/cluster/arm3-queue-job.yaml
kubectl create -f network/cluster/arm3-perm031-job.yaml
```

## Why prep exists at all

`network/cache/` is gitignored, so a clone has no
`ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds`. That file is the
base-001 reference decomposition, and **seven scripts `quit(status=1)` the
moment it is absent** — perm-031, hub-036, hub-037, comm-038, comm-039,
cent-040, tob-045. It exists on the iMac only, because base-001 built it there.

`prep_cache.R` rebuilds it from base-001.R's parameters verbatim and **refuses
to write it unless every module size reproduces base-001 exactly**
(`reference_module_sizes.csv` is tracked, which makes that check possible). A
reference that silently disagreed would score every downstream row against a
decomposition the arm's existing results never used — worse than no rebuild.

Prep also tars the compiled R library to S3, so the 20 fan-out pods untar it in
a minute instead of each compiling WGCNA from source for ten.

## Why perm-031 had to be re-run

At K=20 the smallest attainable empirical p is 1/20 = 0.05. After BH across 71
modules, nothing can clear FDR **regardless of the biology** — so the existing
result (bp_lith below the null in 19 of 20 draws, p = 0.05, BH = 0.36) is
uninterpretable rather than negative. At K=200 the floor is 0.005, which can
clear BH, and the test can finally answer either way.

Watch whether `purple`, `blue`, `pink` and `grey60` — the four below all 20
draws now — stay below the null at K=200. If they do, that is the arm's first
multiple-testing-survivable lithium-specific structural result.

## How work is divided without coordination

A job's pods cannot reach cluster services, so there is no lock to take.
Instead both fan-out jobs use a **static partition on pod index**, read from
the pod's own name via the downward API:

- `arm3-queue` sorts PENDING rows by declared cap descending and deals them
  round-robin, so the long rows spread across pods rather than piling on one.
- `arm3-perm031` gives each pod a contiguous slice of draws 21–200.

A static deal can neither double-run a row nor drop one.

Every pod uploads each checkpoint as soon as it lands (perm-031 syncs every
120 s in the background), so a pod lost at hour three still contributes
everything it finished.

## Getting results back

```bash
python network/cluster/s3ls.py              # everything in the bucket prefix
python network/cluster/s3ls.py runs/        # just finished rows
python network/cluster/collect.py           # status: what landed vs what was asked for
python network/cluster/collect.py --fetch   # download and unpack into network/runs/
```

`collect.py` reconciles S3 against `network/RUNLOG/queue.tsv`, reports which
PENDING rows still have no result, checks the perm-031 draw range for gaps, and
prints the last line of each row's log.

## After the jobs finish

1. `python network/cluster/collect.py --fetch`
2. Aggregate perm-031 across all 200 draws — the pods ran in slice mode and
   exited without aggregating, by design:
   ```bash
   PERM_K=200 Rscript network/experiments/perm-031.R perm-031
   ```
   With every draw present it skips straight to the empirical p-values.
3. Update the queue statuses, then `Rscript network/build_index.R` and
   `Rscript network/build_catalog.R` to regenerate `RANKING.md` and the
   indices — never merge those, regenerate them.
4. Write findings into `network/findings/` as usual.

## Operational notes

- **Quota.** The guaranteed quota (24 CPU, 64Gi) is shared with the prediction
  arm, which can hold all of it. All three jobs run at `priorityClassName:
  burst`, which has its own quota and a 24 h ceiling on the cpu queue. If a job
  sits `Pending` with a `FailedCreate` quota event, that is why.
- **Per-job CPU ceiling is 24**, summed over every container of every replica.
  Not per pod.
- **Pods and their logs are deleted a minute after the job ends.** Anything you
  want must be in S3 before the job exits; each job's `trap ... EXIT` uploads
  its log even on failure.
- **Dry-run before submitting:** `kubectl create --dry-run=server -f <file>`
  checks every cluster rule without running anything.
