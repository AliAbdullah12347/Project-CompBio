# base-001 complete: modules beat gold under subsampling — but Zsummary is not comparable across arms

> **CORRECTION 2026-10-25 — Part 1 is invalid; Part 2 stands.**
>
> The 20 bp_lith draws were computed with RAW test matrices. `base-001.R:344` built
> them as `t(d$logtpm[, sel_sids])`, bypassing `expr_for()` and so the per-sample
> offset removal, while the reference came from the offset-removed cache. Every
> bp_lith number below — gold 4.06, "46 of 69 beat gold in all 20 draws",
> mediumorchid at -1.27 — is therefore reference-vs-raw-test and must not be quoted.
>
> **The bp_nolith half survives**: reference and bp_nolith test both went through
> `expr_for`, so gold Zsummary 14.01 at medianRank 72 of 72, with 63 of 69 modules
> above it, is a valid result.
>
> **Part 2 also survives and is strengthened.** Zsummary genuinely is not comparable
> across these arms; the raw/offset-removed mismatch is a second, independent reason
> on top of the 50-vs-500 permutation difference.
>
> Fixed in 36 places across 32 scripts; see
> `2026-10-25-raw-test-matrices-bypassed-normalisation.md`. base-001's lith draws are
> archived and the row is re-queued.


**Date:** 2026-10-22
**Analysis:** `network/runs/base-001` (offset-removed, 71-module reference, 20 draws)
**Evidence:** 10/25 · **Interest:** 5/10  _(corrected; was 22/25 and 8/10)_ · **Δ:** new

## Part 1 — the gold reversal survives subsampling

The provisional result from the bp_nolith arm holds once all 20 bp_lith subsample draws
are in. Control reference, test = 74 bp_lith samples drawn from 152, 20 draws:

| | Zsummary | medianRank |
|:--|--:|--:|
| gold (random genes), mean over 20 draws | **4.06** (sd 0.93) | **68.5 of 72** (sd 1.2) |

| modules (69 non-grey) | count |
|:--|--:|
| beat gold in **all 20** draws | **46** |
| beat gold in ≥18 of 20 | 55 |
| beat gold in <10 of 20 | 7 |

Median per-module Zsummary sd across draws is 0.98 — the same order as gold's 0.93, so
the separation is not an artifact of unstable draws. Top: blue 26.73, brown 19.72,
yellow 19.23. This settles finding #5 against the version that was marked superseded on
2026-10-18 from the nolith arm alone, and it does so with the subsampling variance that
the arm's earlier 5-draw over-claim lacked.

**One module is genuinely not preserved: mediumorchid, Zsummary mean −1.27, beating gold
in 0 of 20 draws.** Six others beat gold in fewer than half the draws
(yellow4, lightcyan, skyblue2, darkslateblue, mediumpurple2). These are the candidates
for "wired differently in lithium-treated patients" and are where the arm's actual
question lives.

## Part 2 — the arms cannot be compared by Zsummary

This is the part that matters for the remaining queue. The two halves of base-001 were
run with different instruments:

| | permutations | test set |
|:--|--:|:--|
| bp_lith | **50** per draw | 74 subsampled from 152, ×20 |
| bp_nolith | **500** | 74, fixed, single run |

Zsummary is `(observed − mean_perm) / sd_perm`, so it moves with permutation count and
test-set construction. It does, badly, and **not by a constant factor**:

| | bp_lith | bp_nolith | ratio |
|:--|--:|--:|--:|
| gold | 4.06 | 14.01 | 3.45 |
| blue | 26.73 | 48.09 | 1.80 |
| brown | 19.72 | 38.59 | 1.96 |
| yellow | 19.23 | 58.63 | 3.05 |

Ratios from 1.80 to 3.45 mean this is not a rescale that a constant could undo.

**Consequence: "is module X less preserved in bp_lith than in bp_nolith?" cannot be
answered by comparing these two Zsummary columns** — and that comparison is the arm's
central question. Any such claim made from this run would be measuring the permutation
budget.

What remains valid:
- **within-arm** comparison against gold, which is what Part 1 does
- **medianRank**, which is ordinal within a run and so immune to the scaling
- cross-arm comparison **after** re-running with matched perms and matched test-set
  construction

## Recommended before the cross-group question is attempted

Match the instruments. Either raise bp_lith to 500 perms per draw (20× the cost), or run
bp_nolith as 20 subsampled draws at 50 perms the way bp_lith was. The second is far
cheaper and makes the arms symmetric. `base-002-ceil` — the split-half ceiling, next in
the queue — is the right place to calibrate what a *null* difference looks like under
whichever instrument is chosen.

E=5 R=5 B=4 D=5 C=3
