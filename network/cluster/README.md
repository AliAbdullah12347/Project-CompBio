# Arm 3 on the research cluster

Everything in `network/NEXT_STEPS.md` that is bulk compute, moved off the iMac.

## State on 9 Nov 2026, 12:00 UTC

| job | what | pods | status |
|:--|:--|--:|:--|
| `arm3-prep` | rebuild the base-001 reference cache, publish the R library | 1 | Completed 8 Nov |
| `arm3-queue` | first drain attempt, 31 rows over 6 pods | 6 | **Failed** — 2 rows of 31 |
| `arm3-perm031` | perm-031 from K=20 to K=200 | 12 | Running, 80+ draws in |
| `arm3-run-01..07` | the re-drain: 29 rows, one per pod | 29 | Running since 11:45 |

## What the first drain actually did, and why

`arm3-queue` reported 30 rows uploaded, which looked like success. It was not:
the per-row return codes in `status/queue-*.tsv` say **2 rows finished, 19 hit
their time cap, 9 errored, 2 never ran**. Three independent causes.

### 1. WGCNA ran 31 threads inside a 4-CPU cgroup

Every timed-out row's log ends on the same line — `..Working with set 1 as
reference set`, inside `modulePreservation` — and every one of them says
`Allowing parallel execution with up to 31 working processes`.

`parallel::detectCores()` reports the **host's** 32 cores, not the 4 the pod is
allowed, and `enableWGCNAThreads()` takes `nCores - 1`. So each row ran 31
threads against a 4-CPU quota and spent its time context-switching: **17 to 25
minutes per preservation call**, against the ~80 s per 50 permutations recorded
in `gene-027.R`'s header. At 21 calls per row that is 6 to 9 hours against a
3 h cap.

`write_rprofile` in `lib.sh` corrects `detectCores()` to the pod's allocation.
R reads `.Rprofile` from the working directory at startup, so the correction is
in place before `library(WGCNA)`. This is not a workaround — inside a cgroup
`detectCores()` is simply wrong, and R has no way to see the limit.

### 2. Rows were five deep on a shared pod

A validating policy caps a job at **24 CPU in all** — 6 pods × 4 CPU — and at
**86400 s** on the burst priority class. Thirty-one rows through six pods is
five rows per pod behind a 3 h per-row `timeout`.

But the burst quota allows **200 pods across 50 jobs with no CPU ceiling of its
own**. Parallelism comes from more jobs, not bigger ones. `mkjobs.py` emits one
job per six rows, one row per pod, each row holding the full 24 h window.

### 3. Nine rows failed in seconds, for nine ordinary reasons

| row | cause | fix |
|:--|:--|:--|
| `comm-038`, `comm-039` | `ref_net$moduleColors` is unnamed — base-001 built it with `labels2colors(cutreeDynamic(...))` | both rows name their own copy from `rownames(d$logtpm)` |
| `tob-045` | logical mask passed to `expr_samples()`, which takes ids | `colnames(d$logtpm)[keep]` |
| `net-008` | `TOMType = "signed hybrid"`, which WGCNA rejects | `TOMType = "signed"` |
| `bmind-046` | `de_analysis/scripts` missing from the sparse checkout | added |
| `inp-024` | `limma` not in the prebuilt library | `setup_r` tops up missing packages |
| `filt-044` | `counts_474.tsv.gz` is untracked | ships via S3, like the cache |
| `hub-037`, `bmind-047` | read a sibling row's output, dealt to a different pod | chained onto one pod with `+` |
| `cent-040` | its pod was killed before reaching it | re-running |
| `ann-044` | **had no script** — the drain only deals rows with one | written 9 Nov |
| `net-007` | left `RUNNING` by the iMac with no output; the drain filtered on `PENDING` | re-running |

**The colours are not named on purpose.** base-001 cached them positionally and
the preservation rows consume them that way, so naming them *in the cache*
would change what `modulePreservation` matches on. The two rows that need gene
ids attach them to their own copy.

### A note on the leftover pod

`queue-6.log` says `pod 6/8` while the others say `pod N/6`: a pod of the
first, OOMKilled 8-pod submission outlived the job delete and duplicated three
rows. A static deal cannot double-run a row *within* a job, but it can across
two generations of one. Confirm `kubectl get pods` is clear before resubmitting
under a name you have used before.

## The re-drain

```bash
export KUBECONFIG=./kubeconfig-user-aabdullah.yaml
python network/cluster/mkjobs.py            # write the yaml, submit nothing
python network/cluster/mkjobs.py --submit   # write and create
python network/cluster/mkjobs.py --only "net-007" --start 7 --submit   # one row
```

`restore_run` is the cheapest of the fixes. Every long row checkpoints per draw
and skips draws whose checkpoint exists, and the timed-out rows uploaded their
checkpoints before dying — so the re-drain **resumes**: `gene-025` with 18 of
20 draws already done, `gene-029` with 9, `mod-013` with 8. Without it the
resubmission would repeat that work and time out in the same place.

## Why perm-031 had to be re-run

At K=20 the smallest attainable empirical p is 1/20 = 0.05, so after BH across
71 modules nothing can clear FDR regardless of the biology — the existing
result (bp_lith below the null in 19 of 20 draws, p = 0.05, BH = 0.36) is
uninterpretable rather than negative. At K=200 the floor is 0.005.

Watch whether `purple`, `blue`, `pink` and `grey60` — the four below all 20
draws now — stay below the null. If they do, that is the arm's first
multiple-testing-survivable lithium-specific structural result.

The 12 pods each own a contiguous slice of draws 21–200 and upload every
checkpoint as it lands, so a pod lost at hour 20 still contributes what it
finished. They were at ~80 draws after 12 h, which will not reach 180 inside
the 24 h ceiling; aggregate over whatever draws exist. **Draws are exchangeable
and the subset is determined by which slices finished, not by any outcome, so a
partial set is a valid null** — it only raises the p-floor to 1/(K+1).

## Getting results back

```bash
python network/cluster/s3ls.py              # everything under the prefix
python network/cluster/s3ls.py runs/        # just finished rows
python network/cluster/collect.py           # what landed vs what was asked for
python network/cluster/collect.py --fetch   # download and unpack into network/runs/
```

`collect.py` counts uploads; it does **not** know whether a row succeeded. For
that read `status/*.tsv`, where each line is `<row> rc=<code> secs=<n>`:
`rc=124` is the `timeout`, `rc=1` is an R error, `rc=0` is a real completion.
The 8 Nov run looked like 30 successes by upload count and was 2 by return code.

## After the jobs finish

1. `python network/cluster/collect.py --fetch`
2. Check `status/*.tsv` return codes before believing any row finished.
3. Aggregate perm-031 over the draws that exist — the pods ran in slice mode
   and exited without aggregating, by design.
4. **Re-run `synth-048`.** It completed on 8 Nov, but from `6/25` robustness
   rows, so its "71/72 modules stable across ≥80% of variants" describes six
   variants. It must be re-run once `arm3-run-01..05` land.
5. Update `network/RUNLOG/queue.tsv`, then `Rscript network/build_index.R` and
   `Rscript network/build_catalog.R` to regenerate `RANKING.md` and the
   indices — never merge those, regenerate them.
6. Write findings into `network/findings/`.

## Operational notes

- **Per-job ceiling is 24 CPU**, summed over every container of every replica,
  and **86400 s** `activeDeadlineSeconds` on burst. Both are enforced by a
  validating admission policy, so `kubectl create --dry-run=server` catches
  them without running anything. Always dry-run first.
- **Burst has no CPU or memory quota** — only 200 pods and 50 jobs for the
  namespace. More jobs is the lever.
- **22Gi per pod.** The rows that build a full 12,368² TOM need adjacency +
  TOM + dissTOM + the `as.dist` copy, about 6 GB before R's overhead; 5 of 8
  pods were OOMKilled at 7Gi.
- **Pods and their logs are deleted a minute after the job ends.** Anything you
  want must be in S3 first; each job's `trap ... EXIT` uploads its log even on
  failure. A pod killed by the OOM killer never runs its trap, which is why
  `queue-4.log` does not exist.
- `kubectl logs <pod>` works while a pod lives; `kubectl exec` is forbidden,
  and so is listing nodes.
