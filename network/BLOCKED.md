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

