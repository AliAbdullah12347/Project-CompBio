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
| 1 | 11/25 | `phase0-power` | **Scale-free topology fails in this data.** R² never reaches 0.80 at any power 1–20 in any group; negative at low powers (control −0.97/−0.97/−0.93 at p1–3). `pickSoftThreshold` returns 1 / 1 / 2 — the function failing, not a threshold. | 3 | 2 | 2 | 4 | 0 | Measured at the real gene count, n=74/group, before any experiment. Catalogue item 4 says the R² curves are themselves a result. No config/result.json yet, hence C=0; promote when `pow-009` lands. |
| 2 | 6/25 | `phase0-cost` | Reference build 786 s (~23 modules); preservation ~7.8 s/perm; TOM scales n^2.88, preservation n^0.80. | 4 | 0 | 1 | 1 | 0 | Infrastructure, not science. Listed because it set every cap in the queue. |

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
| 1 | 6/10 | 11/25 | **+1.6** | `phase0-power` | WGCNA's central modelling assumption simply does not hold in this blood data, in any group — and the function reports a plausible-looking `powerEstimate = 1` rather than failing loudly. Anyone running a stock WGCNA pipeline here would silently build a degenerate network. |
| 2 | 1/10 | 6/25 | −1.4 | `phase0-cost` | Timing numbers. Interesting only in that they are why the run is survivable at all. |

_Δ = (Interest/10 − Evidence/25) × 10. Positive means the story is ahead of the proof._

_No experiment has completed yet. The rows above are Phase 0 measurements, included so the
first real result has something to be ranked against._
