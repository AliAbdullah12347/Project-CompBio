# Arm 3 — Live Result Ranking

Two rankings, maintained side by side and rebuilt after **every** experiment. Each new
result is re-ranked against everything before it, so both tables are always a live
ranking of the whole run rather than a log.

- **Table A — Evidence** (0–25): how well the run actually *demonstrated* the thing.
- **Table B — Interest** (0–10): how striking, surprising or quotable the finding is.

They are kept separate on purpose. Interest is the honest answer to "did we find anything
cool", and it is worth tracking. But `IMAC_README.md` is explicit that this arm is
exploratory, that "a well-documented negative or inconclusive result is a complete
deliverable", and that "making something look significant is not". If one number mixed
both, a striking-but-thin result would outrank a clean null — which is the exact failure
mode the arm forbids.

**The gap is the useful part.** Interest far above Evidence means *be suspicious and go
get more evidence*. Evidence far above Interest means *solid, unexciting, and still worth
reporting*. The Δ column makes that visible instead of hiding it inside one score.

**Interest never decides what runs next.** Queue order is set by scientific priority in
`RESUME_PROMPT.md`. Do not reorder the queue, extend an analysis, or pick a variant
because it scored high on Interest. Rank it, note it, move on.

---

## Table A — Evidence (the one that governs)

Five components, 0–5 each, 25 total.

| # | Component | 0 | 3 | 5 |
|:--|:--|:--|:--|:--|
| **E** | Evidence | single draw, or descriptive only | Zsummary/medianRank across 100 draws | permutation null K≥500, BH across modules |
| **R** | Relevance | incidental | networks differ by group | bears directly on *does it survive removing composition* |
| **B** | Robustness | untested | holds across ≥2 parameter variants | holds across the filter sweep and/or tobacco imputations |
| **D** | Decisiveness | inconclusive | clear direction | would change what the paper claims — **including a clear null** |
| **C** | Completeness | PARTIAL fragment | complete run, headline numbers | complete + config.json + result.json + README |

Rules: score what the run demonstrated, not what it suggests. A **well-powered** null
scores high on D; a null that merely lacked permutations scores low on E *and* D — say
which it is. Never re-score an old result upward to make a new one look consistent; if a
score changes, say so in Notes and in `METHODS.md`. `PARTIAL` rows cap at C=2.

| Rank | Score | id | What it showed | E | R | B | D | C | Notes |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 1 | 20/25 | `base-001/gold` | **A random gene set scores Zsummary 20.26 — and 3 of 4 real modules score below it.** WGCNA's `gold` module is by definition a random sample of all network genes. In control-vs-bp_nolith (n=74, 500 perms) it scores 20.26, so Langfelder's bands (>10 = strong) put *noise* in the strong band. Only blue (25.72) beats random; turquoise (15.48), yellow (11.01) and brown (6.41) fall below it. | 5 | 5 | 3 | 5 | 2 | **R raised 4->5 at 15 draws:** no module reliably beats random — blue 9/15 draws, turquoise 3/15, yellow and brown 0/15. **Correction:** at 5 draws I reported gold sd 0.93 and called it stable; at 15 it is sd 1.86, range [15.86, 23.58]. Five non-independent draws do not estimate a spread — the same error this run had just warned about. The finding strengthened; only that sub-claim was wrong. Demonstrates empirically what IMAC_README §6.3 states as a caution. **Every Zsummary this arm reports must be quoted against gold, not the published bands** — "Zsummary 14, strongly preserved" would be below this data's random baseline. B=2: one comparison, one group, run still PARTIAL. Does not replace `base-002-ceil`, which isolates the group effect rather than module informativeness. |
| 2 | 18/25 | `base-001/modules` | **18 modules were found, then merged down to 4.** `cutreeDynamic` detected 18 modules on 234 controls at power 12; `mergeCloseModules` at `mergeCutHeight = 0.25` collapsed them to 4, leaving turquoise (6,621) and blue (5,171) holding 95.3% of the transcriptome. That threshold merges any eigengene pair correlating >0.75, and at mean connectivity 8.3% of the genome nearly every pair does. | 5 | 4 | 3 | 4 | 2 | **Supersedes the earlier reading of this run**, which took the 4-module summary line at face value and concluded the data has no modular structure. It does — 18 modules' worth. E raised to 5: the pre-merge count is stated directly in the log, not inferred. `mod-000-resolution` respecified to scan mergeCutHeight — the right axis, and ~10x cheaper since the power-12 TOM is cached. |
| 3 | 15/25 | `base-001/sft` | **Scale-free topology fails at n=234 too — not a small-sample artifact.** Signed R² never reaches 0.80 at any power 1–20 on the full control group and all 12,368 genes; max 0.547 at power 20; `powerEstimate` returns 1. At the power actually used (12) signed R² = **0.162**, slope −0.41, and **mean connectivity = 1,031** — each gene tied to ~8% of the genome. Replicates the n=74 Phase 0 scan in all three groups. | 4 | 2 | 3 | 4 | 2 | From `base-001`'s tracked `sft_r2_curve.csv`, config written before computation. Two sample sizes (74, 234) and two gene sets now agree. C=2: the parent run is still PARTIAL and has no README yet. Promote when `pow-009` adds the other two groups at n=234. |
| 4 | 11/25 | `phase0-power` | Scale-free scan at n=74 per group: R² never reaches 0.80, negative at low powers (control −0.97/−0.97/−0.93 at p1–3); `powerEstimate` 1 / 1 / 2. | 3 | 2 | 2 | 4 | 0 | Superseded in scope by `base-001/sft` but retained: it is the only scan covering **all three groups**, which base-001 does not. |
| 5 | 6/25 | `phase0-cost` | Reference build 786 s (~23 modules); preservation ~7.8 s/perm; TOM scales n^2.88, preservation n^0.80. | 4 | 0 | 1 | 1 | 0 | Infrastructure, not science. Listed because it set every cap in the queue. |

---

## Table B — Interest (tracked, never acted on)

Single 0–10 score. Anchors:

| Score | Means |
|:--|:--|
| 0–2 | Expected. Infrastructure, calibration, or a confirmation nobody would query. |
| 3–4 | Mildly notable. Worth a sentence in the write-up. |
| 5–6 | Genuinely interesting — a real pattern, or a clean null on something people assumed. |
| 7–8 | Surprising. Overturns an expectation, or lands on the arm's central question either way. |
| 9–10 | Would change how someone reads the lithium/bipolar literature. Reserve this. |

| Rank | Interest | Evidence | Δ | id | Why it is interesting |
|:--|:--|:--|:--|:--|:--|
| 1 | 9/10 | 20/25 | **+1.0** | `base-001/gold` | The statistic the entire arm rests on puts random genes in its top interpretive band, in this dataset, by a factor of two — and WGCNA hands you the control that proves it, in the same output table, labelled `gold`. Any paper reporting "Zsummary = 14, strongly preserved" on a network this dense is reporting something below its own noise floor. **Raised to 9 at 15 draws:** the best module in the decomposition beats random in 9 of 15 subsamples and the second-best in 3 of 15 — the decomposition carries almost no information the random control does not. |
| 2 | 8/10 | 18/25 | +0.8 | `base-001/modules` | A single default parameter — `mergeCutHeight = 0.25`, which nobody would think to question — turned 18 modules into 2 covering 95% of the genome, and the pipeline reported nothing wrong. Everything downstream would have been computed on that, and the natural interpretation ("4 modules ≈ 4 blood lineages") is a story that fits an artifact. The near-miss is the interesting part: the summary line and the detail three lines above it support opposite conclusions. |
| 3 | 7/10 | 15/25 | +1.0 | `base-001/sft` | WGCNA's central modelling assumption does not hold in this blood data at **any** sample size or power — and the function reports `powerEstimate = 1`, a plausible-looking number, rather than failing. Anyone running a stock pipeline here builds a degenerate network and never learns. The mean-connectivity figure sharpens it: at the working power each gene is tied to ~8% of the genome, which is not what a co-expression module structure is supposed to look like. |
| 4 | 6/10 | 11/25 | +1.6 | `phase0-power` | Same finding at n=74, across all three groups. |
| 5 | 1/10 | 6/25 | −1.4 | `phase0-cost` | Timing numbers. Interesting only in that they are why the run is survivable at all. |

_Δ = (Interest/10 − Evidence/25) × 10. Positive means the story is ahead of the proof._

_`base-001` is still PARTIAL — its R² curve is a completed, tracked artifact and is ranked
on that basis, capped at C=2 until the parent run finishes. Ranked here rather than at run
completion so the table stays live; the session will re-rank it when `base-001` lands._
