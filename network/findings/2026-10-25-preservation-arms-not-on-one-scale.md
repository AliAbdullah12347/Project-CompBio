# Five experiments done; the decomposition is real, but the cross-group question is still unanswerable

**Date:** 2026-10-25
**Analysis:** `network/analyses/017-preservation-synthesis`
**Evidence:** 14/25 · **Interest:** 5/10 · **Δ:** -1.0

## What is settled

Three preservation arms have now run to 20 draws each on offset-removed expression.
Comparing each arm's modules against **its own** gold (a random gene set drawn from that
same network under that same instrument):

| arm | modules | gold Zsummary | gold rank pct | beat gold in all 20 draws |
|:--|--:|--:|--:|--:|
| ceiling — ctrl160 ref, test = 74 held-out **controls** | 36 | 20.98 (sd 1.66) | 95th | 3 |
| observed — ctrl234 ref, test = 74 **bp_lith** | 71 | 4.06 (sd 0.93) | 95th | 48 |
| comp-adj — composition-residualised ref | 44 | 15.55 (**sd 8.13**) | 100th | 16 |

In every arm gold sits at or near the bottom of the rank distribution. **The
decomposition is real** — which is precisely what finding #5 denied on the raw matrix,
and this now holds across three independently built networks.

Worth flagging: the composition-residualised arm's gold has **sd 8.13**, against 1.66 and
0.93 elsewhere. Residualising for the ILR balances appears to destabilise the network
substantially. That is a result about the method, and it bears on `inp-023` (running now)
and on `tob-045`, which both build on the same residualisation.

## What is NOT settled, and why

**Whether control-vs-bp_lith preservation is lower than control-vs-control.** That is the
arm's central question, and these three runs cannot answer it.

Each arm built its **own reference network** — four distinct caches now exist:

```
ref_control_all_..._offsetrm.rds                       71 modules   base-001
ref_ctrl160_all_..._offsetrm.rds                       67 modules   base-002-ceil  (36 was from the contaminated raw-expression reference)
ref_control_all_..._composition_ilr_b1b4_..._offsetrm  44 modules   inp-022
ref_control_all_..._technical_ilr_age_sex_rin_b1b4_... —            inp-021/023
```

WGCNA names modules by size rank *within a run*, so `blue` in one arm and `blue` in
another are unrelated gene sets. Zsummary values from different reference networks are
not on one scale, and the gold figures above show it directly: 20.98, 4.06 and 15.55 for
the same random-gene control.

**I made exactly this error first.** The initial version of
`analyses/017-preservation-synthesis/run.R` matched modules across arms by colour,
reported "27 of 36 modules below ceiling, paired Wilcoxon p = 1.3e-4" and "composition
adjustment +21.22, p = 2.9e-11", and both are meaningless. The script now refuses the
cross-arm comparison and says why. The numbers are recorded here so they are not quoted
from some later copy.

> **CORRECTED same day.** The section below says the matched test does not exist.
> It does: `base-002-ceil` already scores the ctrl160 reference against BOTH the 74
> held-out controls and bp_nolith at 500 perms each, and against 20 bp_lith draws at
> 50. I had read only its lith draws. The real problem is not a missing experiment
> but that the matched comparison behaves incoherently - see
> `2026-10-25-ceiling-vs-nolith-anomaly.md`.

## The experiment that would settle it

One reference, two test sets, identical instrument:

- reference: the **ctrl160** network (already built and cached)
- test A: the 74 held-out controls — this is `base-002-ceil`, already done
- test B: 74 bp_lith subsampled, 20 draws, 50 perms — **does not exist yet**

With both tests scored against the same 67-module reference (clean expr_for() run; the
36-module count was from the contaminated reference), Zsummary is on one scale and the
question becomes a paired test over identical gene sets. `base-001` cannot serve as
test B because it scores against the 234-control network, not ctrl160.

This is a single ~3 h row and it is the only thing standing between this arm and its
stated question.

## CORRECTION 2026-10-25 — two of three cited arms used raw expression matrices

Commit 1177d8c found that base-001's 20 lith draws and all of base-002-ceil used raw
`t(d$logtpm[...])` rather than `expr_for()` for test matrices. This affects the synthesis
table above:

| arm | status |
|:--|:--|
| ceiling — ctrl160 ref, test=74 controls | **INVALID** (reference also raw; both arms raw) |
| observed — ctrl234, test=74 bp_lith, 20 draws | **INVALID** (lith draws had raw test matrices) |
| comp-adj (inp-022) | **clean** — uses expr_for() throughout |

The claim "The decomposition is real" still holds, supported by:
- inp-022 (ILR-residualised): gold Z mean=15.55, always last in all 20 draws
- inp-023 (technical+ILR): gold always last in 6/6 draws so far (replicating)
- base-001 nolith arm: gold Z=14.01, ranked last (this arm was clean)

The ceiling arm gold Z=20.98 and base-001 lith gold Z=4.06 are from invalid runs and
cannot be quoted. base-002-ceil and base-001 reset to PENDING; will re-run with corrected
scripts.

## CORRECTION 2026-10-26 — base-001 re-run complete; updated synthesis table

base-001 completed 2026-10-26T00:17Z with offset-removed expression for both reference
and test matrices. The synthesis table is partially updated with verified numbers:

| arm | modules | gold Zsummary (lith draws) | gold rank pct | beat gold in all 20 draws |
|:--|--:|:--|:--|--:|
| ceiling — ctrl160 ref, test = 74 bp_lith | **67** | **PENDING** (lith draws starting ~04:40Z) | TBD | TBD |
| observed — ctrl234 ref, test = 74 bp_lith | **71** | **21.74 (sd 7.05)** | **100th (72/72)** | **17/71** |
| comp-adj (inp-022) — ILR residualised | 44 | 15.55 (sd 8.13) | 100th | 16/44 |

Key change for the "sd" comparison: the old base-001 lith gold sd was **0.93** (from the
contaminated run). The verified value is **7.05**. Comparing 7.05 (base-001) to 8.13
(inp-022) shows *both arms have similarly high draw-to-draw variance* — the suggestion
that ILR residualisation "destabilises the network substantially" was based on the
contaminated sd and is retracted.

The claim "the decomposition is real" is now supported by:
- inp-022 (ILR-residualised): gold Z mean=15.55, always last in all 20 draws ✓
- inp-023 (technical+ILR): gold always last in 7/7 draws (PARTIAL; resuming)
- base-001 nolith: gold Z=14.01, ranked 72/72 (last) ✓
- **base-001 lith: gold Z=21.74 (sd 7.05), 72/72 last in all 20 draws** ✓ (new)

The ceiling finding (base-002-ceil lith rows) will update this table when available.
"The experiment that would settle it" (base-002-ceil lith draws) is running; nolith
checkpoint expected ~04:40Z, then lith draws at ~8 min each, cap at 06:17Z.

E=3 R=3 B=3 D=3 C=2
