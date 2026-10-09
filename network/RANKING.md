# Arm 3 — Live Result Ranking

Rebuilt after **every** experiment. Newest scores fold into the same table, so this is
always a ranking of everything the run has produced so far, not a log.

## How a result is scored — and why it is scored this way

Five components, 0–5 each, **25 total**. The rubric deliberately rewards *evidential
strength*, not *excitement*. `IMAC_README.md` is explicit that this arm is exploratory
and that "a well-documented negative or inconclusive result is a complete deliverable",
while "making something look significant is not". A ranking that scored results by how
striking they looked would reward exactly the behaviour the arm forbids — so a clean,
well-powered **null can score 25/25**, and a dramatic finding resting on one draw cannot.

| # | Component | 0 | 3 | 5 |
|:--|:--|:--|:--|:--|
| **E** | Evidence | single draw, or descriptive only | Zsummary/medianRank across 100 draws | permutation null K≥500, BH across modules |
| **R** | Relevance | incidental | networks differ by group | bears directly on *does it survive removing composition* |
| **B** | Robustness | untested | holds across ≥2 parameter variants | holds across the filter sweep and/or tobacco imputations |
| **D** | Decisiveness | inconclusive | clear direction | would change what the paper claims — **including a clear null** |
| **C** | Completeness | PARTIAL fragment | complete run, headline numbers | complete + config.json + result.json + README |

**Scoring rules**
- Score what the run *demonstrated*, not what it suggests. Unsupported extrapolation is 0 on E.
- A null that is **well-powered** scores high on D. A null that merely failed to reach
  significance with too few permutations is low on E and D both — say which it is.
- Never re-score an old result upward to make a new one look consistent. If new work
  changes an old score, say so in the Notes column and in `METHODS.md`.
- `PARTIAL` rows are ranked on what they actually produced, capped at C=2.

## Ranking

| Rank | Score | id | What it showed | E | R | B | D | C | Notes |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| 1 | 11/25 | `phase0-power` | **Scale-free topology fails in this data.** R² never reaches 0.80 at any power 1–20 in any of the three groups; negative at low powers (control −0.97/−0.97/−0.93 at p1–3). `pickSoftThreshold` returns powerEstimate 1 (control, bp_nolith) and 2 (bp_lith) — the function failing, not a threshold. | 3 | 2 | 2 | 4 | 0 | Measured at the real gene count, n=74/group, before any experiment ran. Catalogue item 4 says the R² curves are themselves a result. Not yet a formal run — no config/result.json — hence C=0. Promote when `pow-009` lands. |
| 2 | 6/25 | `phase0-cost` | Reference network build 786 s (~23 modules at power 14); modulePreservation ~7.8 s/permutation; TOM scales n^2.88, preservation n^0.80. | 4 | 0 | 1 | 1 | 0 | Infrastructure, not science. Listed because it set every time cap in the queue and is the reason the run is survivable. |

_No experiment has completed yet. Rows above are Phase 0 measurements, included so the
ranking is never empty and so the first real result has something to be compared against._
