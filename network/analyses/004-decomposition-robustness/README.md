# 004 — Does the lithium dimension replicate across decompositions?

**Why.** Finding 002 rested on one decomposition (`deepSplit=4, mergeCutHeight=0`), which
is what held its robustness score down. All ten cells of the `mod-000-resolution` grid were
rebuilt and asked the same question, with the same a-priori rule each time (Li & Ji
effective dimensions, BH across them, PC chosen by smallest BH — never post hoc).

## Result

| deepSplit | mergeCut | modules | lithium PC | BH | sig | attenuation | PC r² composition |
|--:|--:|--:|:--|--:|:--|--:|--:|
| 2 | 0.00 | 18 | PC2 | 0.0049 | **YES** | 65.7% | 62.6% |
| 4 | 0.00 | **33** | PC2 | 0.0029 | **YES** | 68.1% | **73.3%** |
| 2 | 0.05 | 14 | PC2 | 0.0096 | **YES** | 65.1% | 59.0% |
| 4 | 0.05 | 14 | PC2 | 0.0214 | **YES** | 59.6% | 47.5% |
| 2 | 0.10 | 8 | PC2 | 0.0165 | **YES** | 60.7% | 46.3% |
| 4 | 0.10 | 8 | PC2 | 0.0436 | **YES** | 54.7% | 36.5% |
| 2 | 0.15 | 4 | PC2 | 0.0627 | no | 49.4% | 27.9% |
| 4 | 0.15 | 4 | PC1 | 0.1210 | no | 50.8% | 8.2% |
| 2 | 0.25 | 3 | PC2 | 0.0638 | no | 49.3% | 27.1% |
| 4 | 0.25 | 3 | PC2 | 0.0640 | no | 49.4% | 27.1% |

**Six of ten** decompositions show a lithium-associated dimension at BH < 0.05, and
**all four usable cells** (≥10 modules, largest <25% of genes) are among them. It is
always PC2 except in one degenerate cell.

**None of the six survives composition adjustment.** Applying BH across those six, every
adjusted p is 0.166. The best single cell reaches p = 0.0398 uncorrected, which is why
finding 002 reported "attenuated to p = 0.081" under full adjustment rather than
"abolished" — the same conclusion holds across the grid.

## The unexpected part

**correlation(n_modules, composition share) = +0.879.**

The *better* the decomposition, the *more* of the lithium dimension is cell composition:
73.3% at 33 modules, 27.1% at 3. The degenerate decompositions do not strengthen the
compositional explanation — they dilute it, because collapsing everything into two blocks
blurs the composition-linked structure along with everything else.

This matters for how the finding is read. If the compositional explanation were an
artifact of a particular parameter choice, it should be strongest somewhere arbitrary.
Instead it strengthens monotonically as the decomposition improves, which is what a real
effect being resolved more clearly looks like.

It also means the degenerate `mergeCutHeight=0.25` structure `base-001` used would have
*understated* the composition story, not invented it.

## What this does not establish
Still one dataset, one cohort, one composition estimator. Replication across
decompositions tests parameter robustness, not external validity. And it remains unable to
separate a compositional mediator from a confounder — Arm 2's question.

**Files.** `config.json`, `results/robustness_grid.csv`, `results/summary.rds`.
