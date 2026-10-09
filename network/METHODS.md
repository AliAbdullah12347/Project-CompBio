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
