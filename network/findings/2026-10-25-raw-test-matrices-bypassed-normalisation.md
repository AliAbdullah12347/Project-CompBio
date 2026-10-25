# 36 places scored a raw test matrix against an offset-removed reference

**Date:** 2026-10-25
**Analysis:** `network/lint_expr.sh`, `network/00_load.R`
**Evidence:** 24/25 · **Interest:** 9/10 · **Δ:** new

## Claim

The 2026-10-18 rebuild changed `expr_for()` so every experiment would see
offset-removed expression. **32 of 51 experiment scripts never called it for their test
sets.** They built test matrices directly:

```r
testExpr <- t(d$logtpm[, sel_sids])       # raw log2(TPM+1)
```

So the *reference* network came from offset-removed data while the *test* matrix was
raw. 36 such uses were found and fixed. This is what produced the incoherence reported
the same day in `2026-10-25-ceiling-vs-nolith-anomaly.md`, and it is a different fault
from the one I guessed there (I suspected the Z-block slot indexing at
`base-002-ceil.R:184-188`; that code is correct).

## Damage, exactly

| run | reference | test | verdict |
|:--|:--|:--|:--|
| `mod-000-resolution` | `expr_for` | — | **clean** — the 71-module result stands |
| `inp-022`, `inp-023` | `expr_for` | `expr_for` | **clean** |
| `base-001` | `expr_for` ✓ | bp_nolith `expr_for` ✓ · **20 lith draws raw** ✗ | half invalid |
| `base-002-ceil` | **raw** ✗ | **raw** ✗ | wholly invalid |

**`base-002-ceil` is the worst case and the most instructive.** Its reference was built
from raw expression at line 160, then cached under a key ending `_offsetrm`. The cache
key tagging I added on 2026-10-18 to make stale references impossible to reuse gave
*false assurance* here: the key records the mode the script was supposed to use, not the
data it actually read.

That fully explains the anomaly: cor(log10 module size, Zsummary) was **+0.911** in its
ceiling arm — the dense raw-graph signature, where Zsummary is essentially module size —
against **−0.036** in bp_nolith. It was not a biological result and it was not a slot
bug.

Retracted as a consequence: `base-001`'s "46 of 69 modules beat gold in all 20 draws"
(2026-10-22) used raw lith test matrices. **The bp_nolith half of that finding survives**
— gold Zsummary 14.01, medianRank 72 of 72, 63 of 69 modules above it — because both the
reference and the bp_nolith test went through `expr_for`.

## Fix

- `expr_samples(d, sids)` added to `00_load.R`: subsets by sample ID and applies the
  identical offset removal, computed over all 12,368 genes and the full sample set.
  Verified to reproduce `expr_for()` exactly on a group, and to leave per-sample sd at
  0.0000 against 0.4322 for raw.
- All 36 uses rewritten across 32 scripts; 51/51 still parse.
- `network/lint_expr.sh` fails any experiment reaching for `d$logtpm` as expression
  input. Historical `analyses/*/run.R` are reported as warnings — they ran before the
  rebuild and are already covered by `2026-10-15-library-offset-artifact.md`.
- Contaminated artifacts archived to `network/runs/_raw_test_contaminated_*/`, not
  deleted. `base-001`'s clean `nolith_preservation.rds` was kept, saving a 2 h recompute.

## The general lesson

A cache key records intent. It does not verify content. The `_offsetrm` tag made a raw
network *look* validated for a week. A check that reads the artifact — per-sample sd of
the matrix that built it should be 0 — would have caught this immediately, and is worth
adding to the reference-build path.

E=5 R=5 B=5 D=5 C=4
