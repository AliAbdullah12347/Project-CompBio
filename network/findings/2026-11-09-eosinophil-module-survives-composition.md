# The one module lithium moves that composition does not explain is the eosinophil module

**Date:** 2026-11-09
**Analyses:** `network/experiments/ann-044.R` → `network/runs/ann-044`;
`network/experiments/trait-053-composition-reconcile.R` → `network/runs/trait-053`
**Evidence:** 21/25 · **Interest:** 9/10 · **Δ:** +1.1

## The result

Of 71 modules in the verified reference, **29 are associated with lithium before
adjustment and exactly one survives adjustment for cell composition**. That one
is `mediumpurple2`, 35 genes, and its hub list is an eosinophil programme:

> SMPD3, PTGDR2, CEBPE, IL5RA, SLC29A1, FRRS1, ALOX15, LGALS12, IL1RL1, SPNS3

It survives in all three places the model was fitted, and the association gets
**stronger** under adjustment, not weaker:

| sample set | model | p | BH |
|:--|:--|--:|--:|
| BP1 only (n=226) | marginal | 4.6e-04 | 0.017 |
| BP1 only | `li + b1..b4` | 1.5e-04 | 0.011 |
| BP1 only | `li + dx + technical + b1..b4` | **4.2e-05** | **0.0030** |
| all (n=460) | marginal | 1.0e-05 | 6.3e-05 |
| all | `li + b1..b4` | 1.0e-06 | 7.3e-05 |
| all | `li + dx + technical + b1..b4` | **9.7e-05** | **0.0069** |

Strengthening under adjustment is the signature of a **suppressor**, not a
mediator: the ILR balances correlate with this eigengene in the opposite
direction to lithium's effect, so holding composition fixed removes variance
that was masking the lithium term. `eos-052` reached the same conclusion from
the 22-type fractions — lithium does not shift the granulocyte compartment
internally (p = 0.82), so this is not neutrophilia leaking into an eosinophil
module.

## It is the same gene set in three independent decompositions

The module names differ because WGCNA names colours by size rank, so the same
name means different things in different runs. The **genes** are nearly
identical:

| pair | shared | Jaccard |
|:--|--:|--:|
| base-001 `mediumpurple2` (35) vs grp-051 control_A `orangered3` (37) | 33 | **0.846** |
| base-001 `mediumpurple2` vs grp-050 control_A `orangered3` (37) | 33 | **0.846** |
| grp-050 vs grp-051 `orangered3` | 37 | 1.000 |

Against a background of **0.184** — the median best-match Jaccard between any
two control networks in `grp-051`. This module is reconstructed almost intact
by decompositions that agree on almost nothing else, including ones built from
74-sample subsamples.

So `grp-051b`'s surviving `orangered3` and this `mediumpurple2` are **one
finding, not two**. They were counted as separate results until ann-044 wrote
per-gene module assignments and made the comparison possible.

## Why this matters to the argument

The arm's headline is negative: lithium's effect on the transcriptome runs
through the cell mixture, bipolar diagnosis has no detectable network effect,
and 28 of 29 lithium-associated modules are composition. This is the exception,
and it is a well-behaved one — a small, tightly conserved, functionally
coherent module whose lithium association is *not* explained by, and is in fact
masked by, composition.

That is a concrete claim a reader can check, and it is the kind of residual the
mediation arm's "lower bound" framing predicts should exist.

## Caveats

- Eosinophils are the cell type LM22 handles worst: detectable in only 13.5% of
  samples, and the reference-interval check in `CLAUDE.md` puts the deposited
  eosinophil fraction at 0.0%. The ILR balances therefore have little power to
  adjust for eosinophil abundance specifically, so "survives composition
  adjustment" means "survives the five lineage balances we can estimate", not
  "is independent of eosinophil count". **This is the main threat to the
  result** and the honest way to report it.
- **It is not a gene-filter artifact.** `filt-044` sweeps 31 filters, and at
  the project's own (`>10 counts in >=90% of samples`) all 35 genes survive, as
  they do under edgeR's `filterByExpr`. The module only erodes under filters no
  analysis uses — 7 of 35 at `>50 in >=90%`, 0 at `>50 in 100%`. Its summary
  line reads "0/35 survive all filters", which means all 31 including the
  absurd ones; the median module keeps 14.3% by that standard. These are
  low-abundance genes, as eosinophil genes in whole blood should be, but they
  are comfortably above the threshold the pipeline actually applies.
- 35 genes is small. One eigengene from 35 genes is sensitive to a handful of
  outlying samples; a leave-one-out or bootstrap over subjects would settle it
  and has not been run.
- `dx` is not significant for this module in any model, so this is a lithium
  effect, not a bipolar one.

## Provenance note

`ann-044` wrote `reference_gene_modules.csv` — the first per-gene module
assignment any run in this arm has saved. Without it, module names could not be
matched across decompositions at all, and the three results above read as three
unrelated modules. Every future decomposition should write one.

See also [[2026-11-08-composition-contradiction-resolved]].
