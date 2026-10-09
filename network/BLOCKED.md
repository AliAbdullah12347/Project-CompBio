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

## 2026-10-09T02:12:50Z  SCOPE — index rebuild deliberately omitted

`MASTER_PROMPT.md` Phase 4 step 4 says to run `python index/build_index.py`.
That rewrites `INDEX.db` and `INDEX.csv` **at the repository root**, outside the
`network/` scope this run is restricted to. It is omitted from
`RESUME_PROMPT.md`. The index will go stale for files added under `network/`;
a human should rebuild it after the run.

---
## 2026-10-09T02:46:23Z  NO AUTH TOKEN
`/Users/Aaylab/.arm3_token` is missing or empty.
A human must run `claude setup-token` and store the FULL token there (chmod 600).

## 2026-10-09T02:46:26Z  MALFORMED AUTH TOKEN
`/Users/Aaylab/.arm3_token` is 50 bytes; a real token is ~108 and starts sk-ant-.
Most likely a paste truncated at a line break. Re-store the FULL token.

