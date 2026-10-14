# 002 — Effective dimensionality, and the arm's question answered

**Why this exists.** Analysis 001 found the 33 module eigengenes correlate at median
|r| = 0.822. If they are not independent then testing 33 correlated modules is both
mis-corrected *and underpowered* — the real question is whether lithium associates with
the few dimensions they actually span.

## Effective dimensionality

| decomposition | modules | median abs r | Meff (Li & Ji) | PCs for 80% | PC1 var |
|:--|--:|--:|--:|--:|--:|
| `deepSplit=4, mergeCutHeight=0` | 33 | 0.822 | **6** | **2** | 79.3% |
| `mergeCutHeight=0.25` (base-001) | 4 | 0.569 | 3 | 2 | 64.2% |

Eigenvalue spectrum of the 33 eigengenes: **26.18, 3.43, 1.55, 0.82, 0.59, 0.22, 0.08…**

Thirty-three modules span roughly **six** effective dimensions, and **two** principal
components carry 80% of the variance. This is not an artifact of the merge threshold — the
degenerate 4-module structure gives Meff 3, so *both* decompositions are low-dimensional.
The dimensionality is a property of the data, not of the parameter choice.

## The question, asked of the dimensions rather than the modules

Testing lithium against the 6 effective dimensions (number fixed by Li & Ji before looking
at any association, BH across the 6):

| PC | var | β marginal | p | BH | β comp-adj | p | attenuation |
|:--|--:|--:|--:|--:|--:|--:|--:|
| PC1 | 79.3% | −0.750 | 0.302 | 0.528 | −0.534 | 0.477 | 28.8% |
| **PC2** | **10.4%** | **−0.908** | **0.00048** | **0.0029** | **−0.290** | **0.0398** | **68.1%** |
| PC3 | 4.7% | 0.143 | 0.419 | 0.528 | −0.108 | 0.502 | 24.7% |

**PC2 is associated with lithium.** Permutation null (2,000 label shuffles):
**p = 0.0015** — 2 of 2,000 permutations reached the observed |β|. The association is real.

## And it is cell composition

| nested model | β | p | attenuation |
|:--|--:|--:|--:|
| marginal | −0.9081 | 0.00048 | — |
| + composition (ILR b1–b4) | −0.2896 | 0.0398 | **68.1%** |
| + age only | −0.8561 | 0.00123 | 5.7% |
| + age, sex, RIN, seqPC1–3 | −0.8266 | 0.00162 | 9.0% |
| **FULL: composition + technical** | **−0.2387** | **0.0811** | **73.7%** |

| PC2 variance explained by | |
|:--|--:|
| ILR balances b1–b4 | **73.3%** |
| the 5 lineage fractions | 71.3% |
| lithium alone | 5.3% |

**Technical covariates explain 9% of the association. Age explains 5.7%. Composition
explains 68%.** PC2 is three-quarters cell composition, and once that is removed the
lithium association falls to p = 0.081.

## What this says, stated carefully

At the module-eigengene level, in BP1 (n = 226): **a lithium-associated dimension exists
in the co-expression structure, and it does not survive removing cell composition.**

That is the project's founding hypothesis — that a blood "bipolar signature" may be a
lithium signature acting through neutrophil proportion rather than through any cell
changing behaviour — measured directly on the network side.

## What it does not say

- "Does not survive" means **attenuated to p = 0.081**, not abolished. A residual effect
  of this size cannot be excluded at n = 226.
- One decomposition, one parameter set. `base-003-remerge` asks the *preservation*
  question at the same parameters and may behave differently.
- PC1 (79.3% of variance) is unrelated to lithium; it tracks seqPC1 (r = 0.399) and age
  (r = 0.257) — a technical/demographic axis.
- Composition and lithium are not independent by construction: lithium raises neutrophil
  proportion, so "adjusting for composition" removes part of lithium's own causal path.
  This analysis cannot separate a compositional mediator from a compositional confounder.
  That separation is Arm 2's mediation question, and the right next step is to read this
  against `mediation/results/`.

**Files.** `config.json`, `results/effective_dimensionality.csv`,
`results/lithium_vs_effective_dimensions.csv`, `results/summary.rds`.
