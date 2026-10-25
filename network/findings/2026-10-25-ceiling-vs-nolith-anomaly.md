# The one fully-matched comparison behaves incoherently — stop before interpreting it

**Date:** 2026-10-25
**Analysis:** `network/runs/base-002-ceil`
**Evidence:** 15/25 · **Interest:** 5/10 · **Δ:** -1.0

## The comparison finally exists

`base-002-ceil` scores **one** reference — the 160-control, 36-module ctrl160 network —
against three test sets. Two of them are fully matched:

| test set | n | perms | reference |
|:--|--:|--:|:--|
| 74 held-out controls (`ctrl_ceiling`) | 74 | 500 | ctrl160 |
| bp_nolith (`nolith`) | 74 | 500 | ctrl160 |

Same reference, same gene sets, same instrument. This is exactly the design the
2026-10-25 synthesis finding said did not exist — **that claim was wrong and is corrected
here.** No new experiment is needed to make the comparison; what is needed is to trust it.

## And it does not behave

| | ceiling (controls) | bp_nolith |
|:--|--:|--:|
| gold (random genes) Zsummary | 22.58 | **4.64** (−79%) |
| real modules, median Zsummary | 18.69 | **36.81** (+97%) |
| modules lower in bp_nolith | — | **4 of 36** |
| **cor(log10 module size, Zsummary)** | **+0.911** | **−0.036** |
| Spearman(Z_ceiling, Z_nolith) across modules | — | **+0.072** |
| medianRank, paired Wilcoxon | — | **p = 0.76** |

Three things here are not jointly believable:

1. **Real modules preserve *better* in patients than in held-out controls** (18.69 →
   36.81, 32 of 36 rising). A split-half ceiling is supposed to be the upper bound. A
   disease group exceeding it by ~2× is not a biological result, it is a symptom.

2. **Gold moves the opposite way to the modules** (−79% while modules go +97%) in the
   same comparison, on the same reference, with the same permutation budget.

3. **Zsummary is essentially module size in one arm and unrelated to it in the other**
   (+0.911 vs −0.036). This arm already knows Zsummary tracks size — it recorded
   cor(Zsummary, log10 size) = 0.879 earlier. Seeing that dependence vanish between two
   runs of the same statistic points at the computation, not the biology.

The only statistic that behaves consistently is the ordinal one: gold ranks 36 then 37
of 37 (worst in both), and across real modules the paired Wilcoxon on `medianRank` gives
**p = 0.76 — no detectable difference** between held-out controls and bipolar patients.
But Spearman between the two arms' rankings is only +0.312, so even the orderings
disagree more than they should.

## Extraction code correct — but expression inputs were raw (corrected 2026-10-25)

**2026-10-25 verification (this session):** The extraction code at base-002-ceil.R:182–192
was read end-to-end. Both the ceiling and nolith calls to `modulePreservation` pass
`multiData = list(Ref = ..., Test = ...)` — exact list names `"Ref"` and `"Test"` —
and `extract_mp()` uses the default `ref_name = "Ref"`, `test_name = "Test"`, producing
`"ref.Ref"` and `"inColumnsAlsoPresentIn.Test"`. The names match. Each call is separate
(no multi-test combined call), so there is no slot-swapping between ceiling and nolith.
The Z-block extraction was not the problem.

The correlation figures were re-computed from the checkpoint files and match:

```r
cor(log10(module_sizes), Zsummary)
# ctrl_ceiling:  0.911   ✓ matches finding
# bp_nolith:    -0.036   ✓ matches finding
```

**CORRECTION (commit 1177d8c, 15:40 EDT same date):** The initial conclusion "code correct,
anomaly genuine" was wrong. A concurrent session identified the real cause: `base-002-ceil.R`
line 160 built the reference network from raw expression (`t(d$logtpm[, sids])`) rather
than `expr_for()`, bypassing the per-sample library offset removal introduced 2026-10-18.
All test matrices were also raw. The cache key ending in `_offsetrm` recorded the INTENDED
mode, not the ACTUAL data.

`cor(log10 size, Zsummary) = +0.911` in the ceiling arm is the expected "dense raw-graph
signature": raw expression inflates inter-gene correlations, producing size-dependent
Zsummary. The anomaly (ceiling 2× better preserved than nolith; gold drops 79% while modules
rise 97%) is an artifact of different sample sizes (n=160 reference, n=74 test) interacting
with raw-graph density differences — not biology. The biological mechanism proposed in the
initial draft of this section was wrong and has been removed.

Results archived to `network/runs/_raw_test_contaminated_20261025T193947Z/`. Script
corrected in commit 1177d8c. base-002-ceil reset to PENDING for re-run with correct
offset-removed expression.

## Status of the claims this displaces

- The 2026-10-25 synthesis finding's "the experiment does not exist yet" is wrong;
  corrected there.
- The 2026-10-24 inp-022 "rescue" finding is separately corrected: its cross-arm
  comparison used two different reference networks.

E=4 R=2 B=3 D=3 C=3
