# The central question, asked properly: control and bipolar-off-lithium networks are not detectably different

**Date:** 2026-10-26
**Analysis:** `network/analyses/018-matched-ceiling-comparison`
**Evidence:** 24/25 · **Interest:** 9/10 · **Δ:** new

## The first valid cross-group comparison in the arm

`base-002-ceil` scores **one** reference — the 67-module ctrl160 network, rebuilt on
verified offset-removed expression — against two test groups under one instrument:

| test | n | perms | reference |
|:--|--:|--:|:--|
| 74 held-out controls (`ctrl_ceiling`) | 74 | 500 | ctrl160, 67 modules |
| bp_nolith | 74 | 500 | ctrl160, 67 modules |

Both preconditions were enforced, not assumed: **identical module sets** (68 shared
entries, verified by set equality, not by name) and **matched instrument** (500 perms,
n_test 74 both sides). The size-dependence check now reads +0.838 and +0.729 — consistent
across arms, unlike the +0.911 / −0.036 split that exposed the raw-expression bug.

## Result: no systematic difference

| statistic | ceiling | bp_nolith | paired test |
|:--|--:|--:|:--|
| **medianRank** (ordinal, scale-immune) | median 34.0 | median 33.0 | **Wilcoxon p = 0.724** |
| direction of rank shifts | — | 34 worse, 31 better, 2 tied | **sign test p = 0.804** |
| Zsummary | median 25.07 | median 31.05 | p = 4.8e-05, *wrong direction* |

**By the statistic that can be trusted here, there is no difference.** Two independent
tests agree, and the shift distribution is centred on zero (mean d_rank +0.61, median
+1.0).

The Zsummary result is significant and must be discarded: it says modules preserve
**better** in patients than in held-out controls. A disease group cannot exceed a
split-half ceiling. Gold moves the opposite way to the modules (35.78 → 23.16 while
modules go 25.07 → 31.05), so this is not a common scale shift — it is Zsummary's size
dependence interacting with a large random gene set. In the ceiling, gold outscores the
median real module; in bp_nolith it does not. That inversion drives the p-value.

## The negative is calibrated, which is what makes it usable

21 of 67 modules shift by **20 or more ranks** — in both directions (brown4 −38,
greenyellow +37, maroon +37, orangered4 −36). The per-module shift has **sd 18.2 ranks**.

That is the measurement's noise floor at n = 74, and it is the most useful number here.
It means:

- the absence of a systematic shift is a real constraint, not a failure to measure
- **any single-module claim in this arm must clear ±18 ranks before it means anything** —
  which retrospectively condemns several earlier module-level claims
- a true effect smaller than that is simply invisible at this sample size, and no amount
  of extra permutations fixes it

## What this does and does not say

**Does:** at n = 74 per group, with a verified 67-module reference and a matched
instrument, the co-expression network of unmedicated bipolar patients is not
distinguishable from that of held-out controls by module preservation rank.

**Does not:** say anything about lithium. The bp_lith draws run at 50 perms and are not
comparable to these 500-perm arms — see
`notes/2026-10-26-how-to-read-the-preservation-rows.md`. They are running now and will be
read within-arm against their own gold.

**Does not** rule out an effect below the noise floor. With sd = 18.2 ranks over 67
modules, this design detects only substantial, systematic restructuring.

## Relation to the eigengene findings

This is the preservation instrument finally reporting, and it agrees with the surviving
eigengene result rather than contradicting it. Finding 002 located a lithium-associated
dimension that was mostly cell composition; nothing in the arm has shown structural
network differences that survive proper controls. Three retracted findings had claimed
otherwise on the raw matrix.

E=5 R=5 B=5 D=5 C=4
