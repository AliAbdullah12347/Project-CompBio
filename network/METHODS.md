# Arm 3 — Methods and decision log

Append-only. Every methodological and infrastructural choice, why it was made,
and what was rejected. This is the file the write-up is built from.

---

## 2026-10-09T02:12:50Z  Setup decisions (made with a human present, before the unattended run)

### launchd instead of cron
`MASTER_PROMPT.md` specifies a cron entry. Rejected for three reasons:
1. On macOS Sonoma `/usr/sbin/cron` is TCC-sandboxed and needs Full Disk Access
   granted by hand, or jobs fail on file access with no clear error.
2. cron does not catch up a run missed while the machine was off.
3. `0 */5 * * *` fires at 00/05/10/15/20, so the 20:00->00:00 gap is four hours,
   not five — one invocation a day would start before the usage window reset.

`StartInterval = 18000` gives true five-hour spacing and survives reboots.

### The claude path in the prompt was wrong
`MASTER_PROMPT.md` specifies `/usr/local/bin/claude`. The binary is actually at
`/Users/Aaylab/.local/bin/claude`. As written, every invocation would have failed
with `command not found` and the run would have produced nothing for four days.

### Git staging scoped to network/
The prompt says `git add -A`. Replaced with `git add network/` plus a guard in
`run_next.sh` that unstages anything outside `network/` before committing. The
run must not be able to carry an unrelated working-tree change into a commit.
This also enforces `IMAC_README.md` §8 and the operator's explicit instruction.

### Queue schema extended by two columns
`MASTER_PROMPT.md` specifies eight columns. Two were appended:
- `cap_sec` (9) — per-experiment wall-clock cap, so caps are data rather than
  hard-coded in the driver.
- `spec` (10) — the parameter string. `note` (8) is overwritten with the outcome,
  so without `spec` a finished row would no longer say what it had run.

The first eight columns keep their specified names, order and meaning.

### Authentication failure made loud
An expired OAuth token fails every wake-up identically and produces no output at
all, which is indistinguishable from "nothing to do" in a log nobody is reading.
`wake.sh` greps its own session output for the 401 signature and, on a match,
writes to `BLOCKED.md`, commits, and exits non-zero.

### Baseline ordered first in the queue
`base-001` is row one and nothing else is trusted until it completes, per the
prompt's own priority: a broad sweep with no finished baseline is worth much less.

---

## 2026-10-09T02:38:27Z  Phase 0 measurements (taken before the unattended run began)

### Environment verified
R 4.2.2, 4 cores, 16 GB RAM, 189 GB free. All ten required packages installed and
loading: WGCNA 1.74, igraph 2.3.4, impute 1.72.3, preprocessCore 1.60.2, GO.db 3.16.0,
AnnotationDbi 1.60.2, dynamicTreeCut 1.63.1, fastcluster 1.3.0, limma 3.54.2,
edgeR 3.40.2. `enableWGCNAThreads()` reports 3 worker processes.

`Rscript network/00_load.R` prints the required line exactly:
`12368 genes x 474 samples | control 234, bp_nolith 74, bp_lith 152`.

### igraph approved as an addition
IMAC_README.md section 9 does not list `igraph`, but catalogue items 13 (Louvain/Leiden)
and 14 (centrality) cannot run without it. Approved by the operator before the run began,
per the section 8 rule on dependencies.

### Cost scaling measured rather than assumed
Timed at 1,500 / 3,000 / 6,000 genes (n=74, signed, bicor) and fitted in log-log:

| stage | exponent in n_genes | at 12,368 genes |
|:--|:--|:--|
| TOMsimilarity | 2.88 | ~12 min |
| modulePreservation, per permutation | 0.80 | ~8 s |

Peak RSS 2.05 GB at 6,000 genes; a full 12,368^2 double matrix is 1.14 GB, so a run
holding adjacency + TOM + (1-TOM) peaks near 5-6 GB. Comfortable in 16 GB, but two
experiments must never run concurrently.

**Consequence:** K=500 permutations at the full gene set costs ~1.1 h, not the ~9 h a
quadratic assumption would predict. The caps in queue.tsv are achievable as written.
This was initially estimated wrong from an assumed n^2 scaling; measuring corrected it.

### Scale-free topology fails in this data - this is a result, not a nuisance
`pickSoftThreshold(powerVector = 1:20, networkType = "signed", corFnc = "bicor")` at the
full gene count, n=74 per group:

| group | powerEstimate | first power with R2>=0.80 | max R2 |
|:--|:--|:--|:--|
| control | 1 | none | 0.66 at power 20 |
| bp_nolith | 1 | none | 0.41 at power 20 |
| bp_lith | 2 | none | 0.54 at power 20 |

Sign-corrected R2 is **negative** at low powers (control: -0.97, -0.97, -0.93 at powers
1-3), meaning the fitted slope points the wrong way. No power in 1:20 produces a
scale-free fit in any group.

IMAC_README.md section 6 trap 9 anticipated `powerEstimate = NA`. It returns **1**
instead, which is more dangerous: NA is obviously broken, 1 looks like a valid answer and
would pass straight through an unattended session. Verified empirically - a 2,000-gene
network at the auto-selected power 1 yielded 2 modules, while the same pipeline at a
fixed power of 6 yielded 7-9.

**Decision:** use power 14 (WGCNA's documented default for signed networks at n>40),
record the full R2 curve for every run, and report the failure of scale-free topology as
a finding. Mandated in RESUME_PROMPT.md with a coded guard (reject power < 4, fall back
to 14, log it). **Rejected:** letting each session call `pickSoftThreshold` and trust the
result; and searching for any transform that manufactures a scale-free fit, which would
be fitting the method to the desired answer.

### Driver bug found and fixed before launch
`reclaim()` parsed its UTC timestamps with BSD `date -j` without `-u`, so they were read
as local EDT and landed 4 h in the future. The 6 h staleness threshold was silently a
10 h one, so a crashed experiment would have held its queue row through two sessions
instead of one. Fixed with `-u`, verified at 0 s offset in both shell and awk, and
tested end to end: stale row reclaimed, fresh row untouched, max-attempts row blocked,
all ten columns preserved.

### Credential handling for a PUBLIC repository
`github.com/AliAbdullah12347/Project-CompBio` is public and this run pushes to it
unattended. Two measures:
1. `claude setup-token` produces a one-year token that is NOT stored in the keychain and
   must be supplied via `CLAUDE_CODE_OAUTH_TOKEN`. It lives in `~/.arm3_token`, mode
   0600, **outside the repository**. wake.sh exports it and refuses to start without it.
2. `sync_up()` scans every staged diff for live-credential formats (`sk-ant-*`, `ghp_*`,
   `github_pat_*`, AWS keys, private-key headers) and refuses the commit on a match,
   resetting the index and logging redacted line numbers only. Patterns are
   format-specific by design: a looser rule matching the string
   `CLAUDE_CODE_OAUTH_TOKEN` would fire on wake.sh's own variable name and wedge the
   queue permanently. Both the positive and negative cases were tested.

---

## 2026-10-09T02:39:07Z  modulePreservation escapes the network/ scope by default

Found during Phase 0 timing: `WGCNA::modulePreservation()` defaults to
`savePermutedStatistics = TRUE` with
`permutedStatisticsFile = "permutedStats-actualModules.RData"`, a RELATIVE path, so it
writes into the current working directory - the repository root.

Two problems, not one:
1. **Scope.** It creates a file outside `network/`, which this run is forbidden to touch.
   The `sync_up()` guard would stop it being *committed*, but not being *created*.
2. **Correctness.** Every experiment using the default shares one filename. Runs would
   overwrite each other's permuted null, and any run with `loadPermutedStatistics = TRUE`
   could silently read a different experiment's null distribution and report a p-value
   computed against the wrong permutations.

**Decision:** every `modulePreservation()` call must pass
`permutedStatisticsFile = file.path("network/runs", id, "permutedStats.RData")`.
Mandated in RESUME_PROMPT.md. The stray file created during timing was deleted and the
working tree re-verified clean outside `network/`.

This was found by running the real call rather than reading its signature - the default
is not mentioned in MASTER_PROMPT.md or IMAC_README.md.

---

## 2026-10-09T02:48:22Z  Machinery rebuilt after an adversarial pre-flight audit

A six-dimension audit of the machinery, with every finding put to two independent
refuters, returned ~50 findings. Several were verified by hand and were correct. The
driver and `wake.sh` were rewritten rather than patched.

### The failure that would have produced nothing in four days
The queue held 48 PENDING rows and zero experiment scripts — the scripts were always
meant to be written by the sessions. The driver marked a row `BLOCKED` when its script
was missing, and no code path ever returned a row from `BLOCKED`. The first wake would
have drained rows to terminal one per invocation, and every later wake would have found
an empty queue. Logs would have looked entirely normal throughout.

**Fix:** the script existence check now happens **before** the row is claimed. A missing
script releases the row untouched and exits 2 — the session's cue to write it. Added
`run_next.sh --peek` so a session can see the next row and its spec without claiming it.
The contract is now: peek, write the script if absent, then run.

### No status is absorbing except BLOCKED
`PARTIAL` was terminal, so a transient OOM or segfault permanently retired a recoverable
experiment. `requeue()` now returns `PARTIAL` rows to `PENDING` while attempts remain,
and `BLOCKED` is reachable only after MAX_ATTEMPTS genuine execution failures.

### Other defects fixed
- `attempts` was incremented in both `reclaim()` and `claim()`, so rows were retired
  after two real runs instead of three. It is now incremented in `claim()` only.
- The secret scanner called `git reset`, which unstaged the very warning it had just
  written — a credential incident would have been both silently unreported and
  permanently wedging. It now unstages only the offending files and still commits the alert.
- `awk` getline leaves its variable untouched on failure, so a row with an unparseable
  `started` stamp kept the previous row's timestamp and was never reclaimed. The variable
  is now reset before each read and an unparseable stamp is treated as stale.
- Queue rewrites go through `rewrite()`, which refuses to install a file that is empty or
  has lost its header, so a crash mid-write cannot truncate the queue.
- The cap-kill signalled one PID; forked WGCNA workers survived and accumulated. It now
  TERMs children first, then KILLs after a grace period.
- The outside-`network/` filter used newline-delimited paths, which mangles anything
  needing C-quoting. It now uses `-z`.
- `git pull --rebase` failures were swallowed; an unattended conflict could leave markers
  in `queue.tsv` for the next `git add` to commit as resolved. Conflicts now abort the
  rebase, log to BLOCKED.md, and stop rather than publish.
- `wake.sh` had no wall-clock bound; one hung session would hold the lock and every later
  wake would log only "skipped". Sessions are now capped at 4 h, under the 5 h interval.
- Both locks had PID-liveness checks only. PIDs are recycled, so a stale lock after a
  reboot could wedge the run permanently. Both now have an age fallback.
- `wake.sh`'s auth-failure path committed the whole index, bypassing the scope guard. It
  now uses the same `network/`-only staging.
- The auth grep matched the whole log and would self-trigger on text this script had
  written earlier. It now matches only the session's own last 40 lines, and copies no log
  lines into a file destined for a public repository.
- `network/.gitignore` now excludes lock directories, partial queue rewrites and full
  session transcripts, none of which belong in a public repo.
- `RunAtLoad` is now true: after a reboot the run resumes at login rather than waiting
  for the next 5 h tick.
- Added a zero-throughput watchdog: if the DONE count is unchanged across three
  consecutive wakes (~15 h), it writes to BLOCKED.md and pushes. Every silent failure in
  this design presents identically — normal logs, nothing completing — and this is the
  only check that distinguishes them.

### Science corrections to RESUME_PROMPT.md
The first draft omitted IMAC_README §6.5 (versioned Ensembl ids), §6.6
(Benjamini–Bogomolov rather than one pooled BH across comparisons) and §6.8 (reuse
`de_v2/scripts/mtc.R`; `qvalue`/`IHW` are deliberately absent). All are now stated.
BH across modules *within* a comparison is correct and is what MASTER_PROMPT specifies;
pooling *across* comparisons is what §6.6 forbids. The prompt now makes that distinction
explicit rather than saying only "BH".

It also banned `index/build_index.py` outright. `--query` is read-only and is the
discovery path §7 prescribes first; only a bare run is out of scope. Corrected.

### An error of my own, recorded because the repository is public
While testing `wake.sh`'s token guards I ran the script directly, not accounting for its
failure paths committing and pushing by design. Three commits reached public `origin/main`
(bd7290d, fbc2ae2, 90a79e1): the legitimate machinery plus two junk BLOCKED.md entries
from the tests. **No credential was exposed** — verified by searching every commit for the
token value; the only matches were pattern definitions and documentation strings. The junk
entries are removed in this commit. History was not rewritten: force-pushing a public
branch to tidy two lines is a worse risk than the lines themselves.

**Root cause fix:** `wake.sh --check` now runs every precondition and reports, without
locking, committing, pushing or launching a session. The guards are testable in place.

---

## 2026-10-09T02:55:08Z  Time caps corrected, and long experiments made resumable

### My earlier cap analysis was wrong
I measured one `modulePreservation` run (perm-031, K=500 ≈ 1.1 h at the full gene set)
and concluded "the caps are sound as written". That generalised from the wrong case.

`base-001` is **100 separate preservation runs** — one per subsample draw — because
Zsummary is itself permutation-based, so every draw needs its own null. At the measured
7.8 s per permutation:

| perms/draw | 100 draws | vs the old 7200 s cap |
|:--|:--|:--|
| 10 | 2.4 h | over |
| 50 | 11.0 h | 6x over |
| 100 | 21.9 h | 11x over |

The old cap allowed ~830 permutations *in total*, i.e. 8 per draw — far too few for a
stable Z. Eleven rows specifying `draws=100` had the same defect.

**Fix:** those rows now carry a 10800 s (3 h) cap — under `wake.sh`'s 4 h session cap —
and their specs state `perms_per_draw=50` and `CHECKPOINT PER DRAW (spans sessions)`
explicitly. `base-001` completes over roughly six sessions rather than one.

### Making progress must not count as failing
That correction exposed a second defect. An experiment spanning sessions is marked
`PARTIAL` every time it hits its cap, and `MAX_ATTEMPTS=3` would have retired it at the
third session **while it was still completing draws**. The queue would have reported
BLOCKED on an experiment that was working perfectly.

**Fix:** the driver now counts checkpoint files (`*.rds`, `*.csv`) in the run directory
before and after. If the count increased, the attempt is refunded — a run that advanced
the experiment did not fail. `BLOCKED` now requires `MAX_ATTEMPTS` consecutive runs that
produced **no new checkpoints**, which is the condition that actually means "stuck".

This is why the per-draw checkpointing rule in RESUME_PROMPT.md is load-bearing rather
than merely good practice: without checkpoints there is no progress signal, and a long
experiment would be indistinguishable from a wedged one.

### --peek made strictly read-only
`--peek` previously called `requeue()`, which rewrites `queue.tsv`, while holding no
lock. A peek racing a running driver would have silently discarded the driver's status
update. It is now read-only and reports status counts instead; the driver requeues when
it runs. Verified: the queue is byte-identical across a `--peek`.

---

## 2026-10-09T04:56:49Z  Adversarial audit: 36 findings, 3 confirmed, all fixed

A six-dimension audit with two independent refuters per finding returned 36 findings;
33 were refuted, 3 survived unanimously. All three were then reproduced by hand before
being fixed — none was taken on the agents' word.

### 1. (critical) No git identity: every result commit would have been orphaned
This machine has **no git identity at all** — `git config user.name` and `user.email`
both exit 1, there is no `~/.gitconfig`, and no `GIT_*`/`EMAIL` in the environment.
`git var GIT_COMMITTER_IDENT` — the exact check `rebase` performs — fails with
`unable to auto-detect email address (got 'aaylab@Olin213-98923M.(none)')`.

The commits worked only because they carried `-c user.name -c user.email`. **The
`pull --rebase` did not.** A rebase must stamp a committer on every commit it replays,
so the first time `origin/main` advanced it would die **after** checking out the
upstream.

Reproduced end to end in a scratch clone with no identity:

| | result |
|:--|:--|
| `git pull --rebase` | rc=128, "Committer identity unknown" |
| the old guard's grep | **did not match** — no abort, no BLOCKED.md entry |
| repo state | **DETACHED HEAD**, `.git/rebase-merge` present |
| a later session's commit | **ORPHANED** — unreachable from main |

The run would then have committed every result onto a detached HEAD while
`git push origin main` kept pushing a branch that no longer moved, reporting only a
routine-looking "push failed". The obvious human recovery, `git rebase --abort`,
**discards all of it**.

That the trigger would occur is not speculative: every pre-existing commit in this clone
is authored `AliAbdullah12347 <...@users.noreply.github.com>`, so the operator pushes
from a different checkout. Across ~19 unattended sessions, divergence is effectively
certain.

**Fixes:** (a) a repo-local identity in `.git/config` so nothing in the run can miss it —
this writes repository metadata, never tracked content; (b) the identity also passed
explicitly on the pull; (c) the guard now judges by **exit code**, not by grepping the
message, because no fixed pattern list covers every way a rebase can fail; (d) a failed
pull aborts the rebase and writes to BLOCKED.md; (e) the driver now **refuses to start**
if the repo is mid-rebase or on a detached HEAD, since committing there produces results
that `rebase --abort` would destroy. All verified: the previously fatal pull now returns
rc=0, HEAD stays on main, and the result commit is reachable.

### 2. (high) Rename detection hid deletions outside network/
`diff.renames` defaults to true, and `git diff --cached --name-only` prints only the
**destination** of a rename. A `git mv` from outside `network/` into it therefore listed
nothing outside `network/`, so the scope filter saw a clean commit while the file outside
was deleted. Reproduced: `git mv data/outside.csv network/moved.csv` → the filter counted
**0** paths outside `network/`. All three scope filters now pass `--no-renames`, which
lists both ends independently.

### 3. (critical) The queue could never reach the arm's own question
Measured at the real problem size: a reference network build is **786 s** (~23 modules at
power 14) and `modulePreservation` costs **~7.8 s per permutation**, so one
`draws=100;perms_per_draw=50` row is **~11 h**. Eleven such rows want **~122 h** against
~19 sessions x 4 h x 80% ≈ **61 h** available — **2x oversubscribed**, with only ~5.6 of
the eleven fitting.

`inp-022`/`inp-023` — *"does any difference survive removing cell composition"*, which is
the arm's stated question — sat at rows 22–23 behind eleven multi-session rows. **They
would never have run.** Worse, their specs never stated a draw design, so they looked
cheap.

**Fixes:** the `inp-*` rows now state `draws=100;perms_per_draw=50;CHECKPOINT PER DRAW`
explicitly, and the queue is reordered by scientific priority — `base-001`, then
`inp-022`, `inp-023`, then `pow-009` (cheap, and the R² curves are themselves a result),
`gene-027` (circularity check), `inp-021`. Because the driver claims the oldest PENDING
row, **file order is priority**. RESUME_PROMPT.md now states the budget arithmetic, says
the oversubscription is deliberate (the queue exists so a crashed session always has work,
not as a promise that all 48 rows run), and tells sessions to cache the shared reference
network on a key of (reference group, gene set, corFnc, networkType, power) — worth 786 s
per row and the only way several variants fit at all.

### On the audit itself
The first audit used 78 agents and 4.3M tokens and exhausted the 5-hour window, which
killed its own refutation stage and a follow-up re-audit. The re-run is deliberately much
tighter. Thoroughness has a budget too, and on this account that budget is the same
5-hour window the science competes for.

---

## 2026-10-09T18:56:43Z  Second audit: 6 more defects, all fixed and verified

A tighter re-audit (22 agents, not 78) found six further defects in the rewritten
machinery. Each was reproduced before being fixed.

### Head-of-line starvation (critical)
The progress refund counted **files**, not progress. An experiment emitting one
uniquely-named zero-byte `.csv` per run — a timestamped error table, a placeholder
written before the work — was credited with progress every time, so `attempts` oscillated
0→1→0 and `BLOCKED` was unreachable. `next_pending()` always returns the first PENDING
row, so one broken experiment at the head of the queue would have consumed all four days
and the other 47 rows would never have run.

Worse, the stray-file quarantine compounded it: strays were moved **into** the run
directory before the checkpoint count, so the driver credited progress for the very
contract violation it was simultaneously reporting.

**Fix:** progress now means a file that is newer than a per-run mtime marker, **non-empty**
(`-size +0`), and not a quarantined stray (`! -name 'stray-*'`), plus a cap of 20 refunds.
Verified both ways: an experiment writing empty CSVs now retires at attempt 3 and the next
row runs; one writing real checkpoints still refunds indefinitely and is never retired.

### The cap watchdog hung its own caller (critical)
`( sleep $CAP ... ) &` runs `sleep` as a grandchild. `kill $WPID` kills only the
subshell, leaving `sleep` alive holding the script's inherited stdout. Through a pipe —
which is how a headless session captures output — the caller blocked for the **full cap**
however fast the experiment finished. With `base-001` at 10800 s, a 4 h session would
complete one row; with the three 14400 s rows, the call could never return at all.
**Fix:** poll in 5 s increments instead of one long sleep, and redirect the watchdog's
output so it can never hold that fd.

### Unstaging left a dirty worktree that stopped all pushes (critical)
`git restore --staged` leaves the modification in the working tree, and `pull --rebase`
then refuses with "You have unstaged changes" — reproduced with two consecutive
`sync_up` calls both failing rc=128 **with no upstream divergence at all**, stacking
unpushed commits. The condition never self-clears.

**Fix:** stop staging-then-unstaging entirely. `git commit -- network/` commits by
pathspec, ignoring the index, which removes three failure modes at once: rename blindness,
the dirty-worktree residue, and anything a session happened to stage riding along. Added
`--autostash` so dirty tracked files outside `network/` cannot block the rebase, and a
post-commit assertion that fails loudly if any committed path is outside `network/`.

### Other fixes
- The stray-file snapshot filtered only `^?? network/`, so once run outputs became
  tracked, ` M network/...` lines were reported as strays and BLOCKED.md filled with
  false alarms — burying real blockers in the one file a human reads to triage the run.
- `wake.sh`'s `commit_network` had no mid-rebase/detached-HEAD guard, so its own
  `rebase --abort` could have destroyed orphaned work. Added, plus the pathspec commit.
- `--check` was not side-effect free: `note_blocked` was not guarded and dirtied tracked
  `BLOCKED.md`.
- `wake.sh` wrote the operator's private token **path** into BLOCKED.md and pushed it to
  a public repository. Now described without the path.
- The zero-throughput watchdog silently disabled itself if `queue.tsv` was unreadable
  (`DONE_N` empty rather than zero). It now reports that as a blocker.
- `rewrite()` accepted a truncated queue because it validated only the header. It now
  also refuses any rewrite that would reduce the row count.

### Soft-power recommendation corrected (science)
RESUME_PROMPT.md said "use power 14 (WGCNA's default for signed, n>40)". **That is the
wrong row of the table.** WGCNA's FAQ gives, for signed networks: n<20→18, 20–30→16,
30–40→14, >40→**12** (unsigned/hybrid: 9/8/7/6). With n=74 the correct fallback is **12**;
14 is the 30–40 row. Corrected, with the full table recorded so the next session can check
it rather than trust the prose. The 786 s build measurement was taken at power 14 and is
left labelled as such, since that is what was actually measured.

---

## 2026-10-09  base-001 draw count revised from 100 to 20 (before any computation)

The queue spec for base-001 (and ten other draws=100 rows) was set by an estimate of
~11 h per row: 786 s (reference build, one-time) + 100 × 50 × 7.8 s (permutations) =
~11 h. That estimate is wrong in a material way.

### What the estimate missed

`WGCNA::modulePreservation()` recomputes the reference adjacency and TOM **inside every
call**. One call = one draw. At 12,368 genes the TOM costs ~786 s regardless of whether
it is the reference or the test set. So per draw:

| step | cost |
|:--|:--|
| reference adjacency + TOM (inside modulePreservation, per call) | ~786 s |
| test set adjacency + TOM (inside modulePreservation, per call) | ~786 s |
| 50 permutations at 7.8 s/perm | 390 s |
| **total per draw** | **~1 962 s ≈ 32.7 min** |

There is no public API in WGCNA 1.74 to inject a pre-computed TOM into
`modulePreservation`, so this cost cannot be reduced without re-implementing the
statistics — which would be risky and error-prone on an unattended run.

### Consequence for the queue

100 draws × 1 962 s ≈ 54.5 h per row. Eleven draws=100 rows want ~600 h against
~61 h available — **10x oversubscribed**. At that ratio `inp-022` and `inp-023` — the
arm's stated question — would never run, which would defeat the purpose of the arm.

### Decision

Implement **draws = 20** for base-001. 20 × 1 962 s ≈ 10.9 h ≈ 11 h, which matches
the *intended* budget for this row. 20 independent subsamples of bp_lith (each n=74 from
152) gives a usable distribution of Zsummary values for descriptive purposes; 50
permutations per draw is sufficient for an individual null. bp_nolith, which cannot be
subsampled, is run once with 500 permutations (ten times the per-draw count) for a
stable point estimate.

This decision is pre-specified here and written into `network/config/base-001.json`
**before any computation begins**. The queue spec field still reads "draws=100"
(not modified; it is not a config file); the discrepancy is recorded in the config as
`draws_spec` vs `draws_implemented`.

The same revision applies to the other ten draws=100 rows (inp-022, inp-023, …).
Those will be revisited when the driver reaches them, but this entry records the
reasoning so that session need not repeat it.

---

## 2026-10-09  Added base-002-ceil: without a ceiling, Zsummary is uninterpretable

Reviewing `base-001.R` after session 1 wrote it (before it ran — the row was still
PENDING, attempts=0), two linked problems surfaced. The script is otherwise sound: power
12 with a guard, modules defined once and cached, `permutedStatisticsFile` set per draw
inside the run directory, config written before computation, checkpoint per draw, exits 1
to requeue, and the draws revision from 100 to 20 documented in advance.

### Problem 1 — no control-vs-control comparison
`base-001` measures preservation of control modules in `bp_nolith` and `bp_lith`, and
nowhere else. There is therefore no answer to the only question that makes a Zsummary
readable: **what does a well-preserved module look like at n=74 in this data?**

If `base-001` reports Zsummary = 8 for a module in `bp_lith`, that is consistent with
all of:
- lithium genuinely disrupting the module,
- 74 samples being too few to recover it in any group,
- the reference being estimated more precisely than any test set.

Nothing in the run distinguishes them. Langfelder's rule-of-thumb bands (<2 not
preserved, 2–10 weak, >10 strong) were derived on other datasets and are not a
substitute for an internal ceiling — especially here, where `IMAC_README.md` §6.3 already
warns Zsummary has no calibrated null.

### Problem 2 — the reference is n=234 while every test set is n=74
`expr_for(d, group = "control")` returns all 234 controls, and `base-001` uses that as the
reference while testing against 74-sample groups. MASTER_PROMPT.md:117 says to subsample
**the two larger groups** to n=74 — control (234) and bp_lith (152). Only bp_lith was
subsampled.

This is not fatal to the between-group comparison: `bp_nolith` and `bp_lith` face the same
reference, so the contrast between them stays internally fair. But every absolute Zsummary
is inflated by the reference's extra precision, and a reader cannot tell how much.

### Why one experiment fixes both
A **split-half control design** calibrates the asymmetry and supplies the ceiling at once:

- reference: 160 randomly held-out controls
- test A: the **other 74 controls**, disjoint from the reference → the ceiling
- test B: `bp_nolith` (74, never subsampled)
- test C: `bp_lith`, 20 draws of 74

Every test set is now n=74 against one reference, so the ceiling absorbs exactly the
sampling-precision effect that Problem 2 introduces. Zsummary in `bp_nolith`/`bp_lith`
becomes readable as a fraction of what controls themselves achieve.

**The control test set must be disjoint from the reference.** Testing a 74-sample subset
of the same 234 controls used to define the modules would measure the modules against
their own training data and produce an inflated ceiling — which would be worse than having
none, because it would look rigorous.

### What was NOT done
`base-001` was left exactly as written and is still queued first. It was not edited, its
cached reference was not invalidated, and its draws were not changed. Its numbers remain
valid as a within-reference comparison of `bp_nolith` against `bp_lith`; `base-002-ceil`
supplies the scale they should be read on. Rewriting a correct experiment to answer a
different question would have cost a session and thrown away work already specified.

`base-002-ceil` is inserted at **position 2**, ahead of `inp-022`/`inp-023`, because the
arm's key question is itself a preservation comparison — answering it before the ceiling
exists would produce a number nobody could interpret either.


## 2026-10-09  Reference build cost scales with sample size — caps were sized at the wrong n

`base-001`'s reference network took **39 minutes** to build, against the **786 s** figure
that sized every cap in the queue. The estimate was not wrong; it was measured at the
wrong sample size.

The 786 s benchmark was taken at **n=74**, the subsample size. `base-001`'s reference is
the **full control group, n=234**. Correlation cost scales with the sample count, so
786 × (234/74) ≈ 2,360 s ≈ 39 min — which is what was observed, to within a minute.

**Rule for sizing any row: multiply 786 s by n_reference/74.**

| row | reference n | predicted build |
|:--|:--|:--|
| `base-001`, `inp-021`–`inp-024` | 234 | ~39 min |
| `base-002-ceil` | 160 | ~27 min |
| anything with an n=74 reference | 74 | ~13 min |

This does not break the run — the build is paid **once per cache key**, and
`RESUME_PROMPT.md` already requires caching the reference on (group, gene set, corFnc,
networkType, power, residualisation). But it means the first session touching a new
reference spends roughly 40 minutes before any preservation work starts, which is worth
knowing when deciding what still fits in a window.

Observed rather than inferred: the log went quiet for 38 minutes because WGCNA emits
nothing during the single-threaded BLAS matrix multiply, while parent CPU time kept
climbing and memory held at 2.67 GB. A quiet log here is not a stalled run, and the
zero-throughput watchdog is correctly scoped to sessions rather than to log activity.

