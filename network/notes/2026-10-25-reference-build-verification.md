# Reference builds now verify the matrix, not the label

Added after a reference network built from raw `log2(TPM+1)` sat for a week in a cache
file whose key ended `_offsetrm`. The key records the mode the script was *supposed* to
use; nothing read the data.

## Three layers, build-time to load-time

**1. `assert_expr_mode(M, d, mode, what)` — 00_load.R.** Reconstructs a 400-gene random
probe of `M` via `expr_samples()` under the declared mode and compares values. On
mismatch it reports the max deviation *and* which mode the matrix actually looks like.

It compares values rather than checking "per-sample mean is zero" because the offset is
removed over all 12,368 genes — any gene-subsetted matrix has a non-zero mean and would
fail a naive check. Verified against a 2,000-gene subset.

**2. Called at all 30 `blockwiseModules()` sites** — 29 checks plus `bmind-047`, exempt
and documented: its input is a bMIND per-lineage profile, not a slice of `d$logtpm`.

**3. `ref_cache_ok(path, mode)` gates every load.** All 28 `file.exists(ref_cache)` sites
became `ref_cache_ok(ref_cache)`. It **returns FALSE and quarantines** rather than
raising, so the caller's existing "if not usable, build it" branch just rebuilds.
Raising would wedge the row — failing every attempt and reaching BLOCKED without ever
regenerating the thing it objected to.

A cache with no `expr_mode` stamp predates this and cannot be verified from its contents,
so it is refused rather than trusted. All three pre-existing caches are therefore
quarantined and will be rebuilt; that costs ~17 min each and is the conservative choice.

## Verified

| case | result |
|:--|:--|
| clean offset-removed matrix | passes, max deviation 0 |
| raw matrix claimed as offset-removed | caught, reported as "looks like: raw" |
| 2,000-gene subset | passes |
| unstamped cache | quarantined, FALSE |
| correctly stamped cache | passes, survives |
| cache stamped `raw` under an offset-removed run | quarantined, FALSE |

`network/lint_expr.sh` covers the complementary failure — a script bypassing `expr_for`
entirely. Together: lint catches the bypass, `assert_expr_mode` catches a bad build,
`ref_cache_ok` catches a bad reuse.
