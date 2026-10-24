# Arm 3 — Blockers, gaps and workarounds

Append-only. Newest at the bottom. Every entry says what was needed, what was
tried, and what was done instead.

---

## 2026-10-09T02:12:50Z  SETUP — headless CLI authentication expired

`claude -p` returned:

```
Failed to authenticate. API Error: 401 OAuth access token has expired.
```

Credentials live in the macOS keychain (`Claude Code-credentials`, account
`aaylab`), not in a file. Re-authentication needs an interactive browser flow,
so it cannot be automated.

**Required before the run starts:** a human runs `claude login` in a terminal.

Until then every scheduled wake-up would produce nothing. `network/wake.sh` now
detects this specific failure, appends an entry here, commits and exits non-zero,
so it can never fail silently again.

---

## 2026-10-09T02:12:50Z  SETUP — gfortran library path broken (RESOLVED)

`igraph`, `Hmisc` and therefore `WGCNA` all failed to link:

```
ld: warning: search path '/usr/local/gfortran/lib' not found
ld: library 'gfortran' not found
```

R 4.2's `Makeconf` hard-codes `FLIBS` to CRAN's gfortran 8.2.0 at
`/usr/local/gfortran`, which is not installed. This machine has Homebrew
GCC 14.2.0 instead.

**Fix:** `~/.R/Makevars` overrides `FLIBS` to
`-L/usr/local/opt/gcc/lib/gcc/current -lgfortran -lquadmath -lm`.
The `opt` path is Homebrew's version-stable symlink, so a gcc upgrade will not
re-break it. This is a user-level R config outside the repository; nothing in
`Implementation/` was modified.

---

## 2026-10-09T02:12:50Z  SCOPE — index rebuild restricted, not banned

`MASTER_PROMPT.md` Phase 4 step 4 says to run `python index/build_index.py`. A bare run
rewrites `INDEX.db` and `INDEX.csv` **at the repository root**, outside the `network/`
scope, so it is not permitted.

`build_index.py --query "..."` is read-only and IS permitted — `IMAC_README.md` §7 names
it as the first discovery step, so banning it outright (as an earlier draft of
`RESUME_PROMPT.md` did) would have broken the prescribed workflow. The distinction is now
explicit in the prompt. The index will go stale for files added under `network/`; a human
should rebuild it after the run.

---
## 2026-10-09T02:47:51Z  MALFORMED AUTH TOKEN
`/Users/Aaylab/.arm3_token` is 50 bytes; a real token is ~108 and starts sk-ant-.
Most likely a paste truncated at a line break. Re-store the FULL token.

## 2026-10-09T22:16:26Z  base-001  hit 10800s cap; partial results kept; +21 new checkpoints (attempt refunded)
```
     ..Working with set 1 as reference set
[2026-10-09T18:05:38] draw 16/20 done: Zsummary range [4.57, 20.3]
[2026-10-09T18:05:38] draw 17/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-09T18:10:01] draw 17/20 done: Zsummary range [5.58, 24.37]
[2026-10-09T18:10:01] draw 18/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-09T18:14:22] draw 18/20 done: Zsummary range [3.53, 21.6]
[2026-10-09T18:14:22] draw 19/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
CAP EXCEEDED after 10800s
```

## 2026-10-10T01:36:01Z  PULL --REBASE FAILED (rc=1), ABORTED
Local commits are intact on main and unpushed. A human must reconcile.
```
fatal: unable to access 'https://github.com/AliAbdullah12347/Project-CompBio.git/': Could not resolve host: github.com
```

## 2026-10-14T22:14:58Z  base-001  hit 10800s cap; partial results kept; +5 new checkpoints (attempt refunded)
```
     ..Working with set 1 as reference set
[2026-10-14T17:38:20] draw 1/20 done: Zsummary range [5.43, 34.07]
[2026-10-14T17:38:20] draw 2/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-14T17:50:31] draw 2/20 done: Zsummary range [5.48, 30.02]
[2026-10-14T17:50:31] draw 3/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-14T18:02:43] draw 3/20 done: Zsummary range [4.76, 31.38]
[2026-10-14T18:02:44] draw 4/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
CAP EXCEEDED after 10800s
```

## 2026-10-14T22:15:00Z  PULL --REBASE FAILED (rc=1), ABORTED
Local commits are intact on main and unpushed. A human must reconcile.
```
fatal: unable to access 'https://github.com/AliAbdullah12347/Project-CompBio.git/': Could not resolve host: github.com
```

## 2026-10-15T09:08:10Z  NO PROGRESS IN 3 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 3 wake-ups (~15 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-15T14:08:03Z  NO PROGRESS IN 4 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 4 wake-ups (~20 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-15T19:08:07Z  NO PROGRESS IN 5 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 5 wake-ups (~25 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-15T20:40Z  ROOT CAUSE IDENTIFIED AND RESOLVED (session 6)

The "no progress" entries above were caused by two compounding failures:

1. **2026-10-14T22:15Z**: `git push` failed (`Could not resolve host: github.com`) after
   base-001 hit its cap and the driver tried to sync. The pull --rebase also failed.
   The queue was left with base-001=PARTIAL and 3 lith draws saved, but unpushed.

2. **2026-10-15T14:05Z and 19:05Z**: Daemon sessions started but the Claude CLI
   immediately exited with `API Error: Unable to connect to API (FailedToOpenSocket)`.
   Both sessions hit the same connectivity failure (~5-min window each time), so no
   work was done and no progress.jsonl entries were written.

3. **2026-10-15T09:08Z**: A session ran but hit the pull --rebase failure again
   (github.com still unreachable).

**Current status (20:40Z):** Network restored. Session 6 is running. base-001 is now
actively computing lith draws (draw 4 done at 20:37:56Z; draws 5–15 in progress).
Queue is RUNNING. The missing ~25 hours was a network outage, not a science or
infrastructure failure.

## 2026-10-15T20:57:34Z  NO PROGRESS IN 6 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 6 wake-ups (~30 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-16T00:08:05Z  NO PROGRESS IN 7 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 7 wake-ups (~35 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-16T04:08:03Z  NO PROGRESS IN 8 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 8 wake-ups (~40 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-16T09:08:01Z  NO PROGRESS IN 9 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 9 wake-ups (~45 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-16T14:08:00Z  NO PROGRESS IN 10 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 10 wake-ups (~50 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-16T19:08:06Z  NO PROGRESS IN 11 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 11 wake-ups (~55 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-17T00:08:11Z  NO PROGRESS IN 12 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 12 wake-ups (~60 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-17T04:08:06Z  NO PROGRESS IN 13 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 13 wake-ups (~65 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-17T09:08:08Z  NO PROGRESS IN 14 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 14 wake-ups (~70 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-17T14:08:00Z  NO PROGRESS IN 15 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 15 wake-ups (~75 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-17T19:08:04Z  NO PROGRESS IN 16 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 16 wake-ups (~80 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-18T00:08:08Z  NO PROGRESS IN 17 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 17 wake-ups (~85 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-18T04:08:04Z  NO PROGRESS IN 18 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 18 wake-ups (~90 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-18T09:08:03Z  NO PROGRESS IN 19 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 19 wake-ups (~95 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-18T14:08:01Z  NO PROGRESS IN 20 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 20 wake-ups (~100 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-18T19:07:57Z  NO PROGRESS IN 21 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 21 wake-ups (~105 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-18T21:38:32Z  base-001  exit 1 after 0s; see log.txt
```
[2026-10-18T17:38:31] === base-001 starting ===
Error in paste0(CACHE_KEY, EXPR_TAG) : object 'EXPR_TAG' not found
Execution halted
```

## 2026-10-19T00:16:58Z  NO PROGRESS IN 22 CONSECUTIVE SESSIONS
DONE count has been stuck at 1 across 22 wake-ups (~110 h).
The scheduler and sessions are running, but nothing is completing.
Check the newest entries above and `network/RUNLOG/progress.jsonl`.

## 2026-10-19T00:39:44Z  base-001  hit 10800s cap; partial results kept; +7 new checkpoints (attempt refunded)
```
     ..Working with set 1 as reference set
[2026-10-18T20:09:30] draw 2/20 done: Zsummary range [-0.69, 30.35]
[2026-10-18T20:09:30] draw 3/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-18T20:21:55] draw 3/20 done: Zsummary range [-0.55, 30.24]
[2026-10-18T20:21:55] draw 4/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-18T20:36:34] draw 4/20 done: Zsummary range [-1.14, 25.89]
[2026-10-18T20:36:34] draw 5/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
CAP EXCEEDED after 10800s
```

## 2026-10-19T04:33:20Z  base-001  hit 10800s cap; partial results kept; +15 new checkpoints (attempt refunded)
```
     ..Working with set 1 as reference set
[2026-10-19T00:06:42] draw 17/20 done: Zsummary range [-1.2, 26.59]
[2026-10-19T00:06:43] draw 18/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-19T00:18:27] draw 18/20 done: Zsummary range [-1.69, 27.1]
[2026-10-19T00:18:27] draw 19/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-19T00:30:12] draw 19/20 done: Zsummary range [-1.91, 22.72]
[2026-10-19T00:30:12] draw 20/20: subsampling bp_lith (n=74 from 152) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
CAP EXCEEDED after 10800s
```

## 2026-10-19T04:33:21Z  PULL --REBASE FAILED (rc=1), ABORTED
Local commits are intact on main and unpushed. A human must reconcile.
```
fatal: unable to access 'https://github.com/AliAbdullah12347/Project-CompBio.git/': Could not resolve host: github.com
```

## 2026-10-22T23:06:19Z  base-002-ceil  hit 10800s cap; partial results kept; +4 new checkpoints (attempt refunded)
```
       Calculating new MEs...
[2026-10-22T16:49:03] reference: 36 modules
           used  (Mb) gc trigger  (Mb) limit (Mb)  max used   (Mb)
Ncells  2431991 129.9    4380007 234.0         NA   4380007  234.0
Vcells 16068794 122.6   62814927 479.3      16384 480940418 3669.3
[2026-10-22T16:49:03] CEILING: 74 control preservation (500 perms) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-22T18:21:19] ceiling done: Zsummary range [9.44, 35.46]
           used  (Mb) gc trigger  (Mb) limit (Mb)  max used   (Mb)
Ncells  2455443 131.2    4380007 234.0         NA   4380007  234.0
Vcells 18111401 138.2   62814927 479.3      16384 480940418 3669.3
[2026-10-22T18:21:19] bp_nolith preservation (500 perms) ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
CAP EXCEEDED after 10800s
```

## 2026-10-22T23:06:21Z  PULL --REBASE FAILED (rc=1), ABORTED
Local commits are intact on main and unpushed. A human must reconcile.
```
fatal: unable to access 'https://github.com/AliAbdullah12347/Project-CompBio.git/': Could not resolve host: github.com
```

## 2026-10-23T02:06:24Z  base-002-ceil  hit 10800s cap; partial results kept; +10 new checkpoints (attempt refunded)
```
     ..Working with set 1 as reference set
[2026-10-22T21:39:27] draw 7 done: Zsummary range [2.92, 32.02]
[2026-10-22T21:39:28] draw 8/20: subsampling bp_lith ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-22T21:48:32] draw 8 done: Zsummary range [5.05, 29.64]
[2026-10-22T21:48:32] draw 9/20: subsampling bp_lith ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-22T21:57:38] draw 9 done: Zsummary range [6.64, 34.22]
[2026-10-22T21:57:38] draw 10/20: subsampling bp_lith ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
CAP EXCEEDED after 10800s
```

## 2026-10-23T02:06:25Z  PULL --REBASE FAILED (rc=1), ABORTED
Local commits are intact on main and unpushed. A human must reconcile.
```
fatal: unable to access 'https://github.com/AliAbdullah12347/Project-CompBio.git/': Could not resolve host: github.com
```

## 2026-10-23T03:46:46Z  PULL --REBASE FAILED (rc=1), ABORTED
Local commits are intact on main and unpushed. A human must reconcile.
```
fatal: unable to access 'https://github.com/AliAbdullah12347/Project-CompBio.git/': Could not resolve host: github.com
```

## 2026-10-23T06:46:50Z  inp-022  hit 10800s cap; partial results kept; +7 new checkpoints (attempt refunded)
```
     ..Working with set 1 as reference set
[2026-10-23T02:13:49] draw 2 done: Zsummary [10, 75]
[2026-10-23T02:13:49] draw 3/20: subsampling residualised bp_lith ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-23T02:25:04] draw 3 done: Zsummary [10.86, 81.51]
[2026-10-23T02:25:04] draw 4/20: subsampling residualised bp_lith ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-23T02:36:19] draw 4 done: Zsummary [13.32, 85.58]
[2026-10-23T02:36:19] draw 5/20: subsampling residualised bp_lith ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
CAP EXCEEDED after 10800s
```

## 2026-10-23T06:46:51Z  PULL --REBASE FAILED (rc=1), ABORTED
Local commits are intact on main and unpushed. A human must reconcile.
```
fatal: unable to access 'https://github.com/AliAbdullah12347/Project-CompBio.git/': Could not resolve host: github.com
```

## 2026-10-23T09:46:54Z  inp-022  hit 10800s cap; partial results kept; +15 new checkpoints (attempt refunded)
```
     ..Working with set 1 as reference set
[2026-10-23T05:14:53] draw 17 done: Zsummary [12.91, 61.2]
[2026-10-23T05:14:53] draw 18/20: subsampling residualised bp_lith ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-23T05:26:07] draw 18 done: Zsummary [7.8, 73.02]
[2026-10-23T05:26:08] draw 19/20: subsampling residualised bp_lith ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-23T05:37:23] draw 19 done: Zsummary [10.4, 71.38]
[2026-10-23T05:37:23] draw 20/20: subsampling residualised bp_lith ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
CAP EXCEEDED after 10800s
```

## 2026-10-23T09:46:54Z  PULL --REBASE FAILED (rc=1), ABORTED
Local commits are intact on main and unpushed. A human must reconcile.
```
fatal: unable to access 'https://github.com/AliAbdullah12347/Project-CompBio.git/': Could not resolve host: github.com
```

## 2026-10-24T07:08:58Z  inp-023  hit 10800s cap; partial results kept; +5 new checkpoints (attempt refunded)
```
     ..Working with set 1 as reference set
[2026-10-24T02:41:18] draw 1 done: Zsummary [7.75, 82.51]
[2026-10-24T02:41:18] draw 2/20 ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-24T02:53:15] draw 2 done: Zsummary [7.67, 78.86]
[2026-10-24T02:53:15] draw 3/20 ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
[2026-10-24T03:05:13] draw 3 done: Zsummary [8.99, 82.61]
[2026-10-24T03:05:13] draw 4/20 ...
      ..checking data for excessive amounts of missing data..
      ..calculating observed preservation values
      ..calculating permutation Z scores
     ..Working with set 1 as reference set
CAP EXCEEDED after 10800s
```

