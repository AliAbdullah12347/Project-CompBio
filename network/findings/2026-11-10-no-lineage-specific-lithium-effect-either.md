# Deconvolving to five lineages does not reveal the lithium effect that whole blood hides

**Date:** 2026-11-10
**Analysis:** `network/experiments/bmind-047.R` -> `network/runs/bmind-047`
**Evidence:** 20/25 · **Interest:** 8/10 · **Δ:** +0.8

## The result

The whole-blood network arm is negative, and the obvious objection is that
bulk expression averages five cell types together, so a lithium effect
confined to one of them would be invisible. `bmind-047` tests that directly:
take the bMIND per-lineage expression estimates, build a separate WGCNA
network per lineage on controls, and run the same preservation test.

**Nothing in any lineage.**

| lineage | modules | bp_lith medianRank | bp_nolith medianRank | p | BH |
|:--|--:|--:|--:|--:|--:|
| gran | 14 | 8.0 | 7.5 | 0.784 | 1.000 |
| mono | 20 | 10.9 | 10.0 | 0.729 | 1.000 |
| T | 16 | 8.5 | 8.5 | 0.821 | 1.000 |
| NK | 16 | 8.2 | 8.5 | 1.000 | 1.000 |
| B | 18 | 9.6 | 9.0 | 0.983 | 1.000 |

Paired Wilcoxon on medianRank, bp_lith averaged over its 20 subsample draws
against the full bp_nolith run, over the modules each lineage shares. **0 of
5 lineages**, and the smallest p is 0.729 - further from significance than
any of the 25 whole-blood parameter variants, whose range starts at 0.441.

Preservation is strong throughout: median lith Zsummary 16.5 to 19.8 per
lineage. One module of the 84 across all five falls below Zsummary 2 -
`lightyellow` in mono, at 1.68 - and it is at 2.04 in bp_nolith, so it is a
weakly preserved module in both patient groups rather than something lithium
takes apart.

## Why this is the test that closes the argument

Every other negative in this arm is a whole-blood negative, and all of them
share one escape route: bulk co-expression is dominated by cell composition,
so a within-cell-type rewiring could be diluted past detection. That is a
serious objection - it is the same reasoning the whole project is built on.

This row removes it. The deconvolution is the DE arm's own, run without the
phenotype so the cell-type estimates are not shaped by the comparison, and
the networks are built per lineage on controls before lithium status is ever
consulted. If lithium restructured co-expression inside granulocytes,
monocytes, T, NK or B cells, this is where it would show, and it does not.

## What it does not rule out

- **Five lineages, not 22.** bMIND was run on the five-way lineage collapse,
  not the full LM22. A rewiring confined to a rare subtype - the eosinophil
  programme `mediumpurple2` points at, for instance - would be inside `gran`
  here and diluted again. The eosinophil result stands on its own evidence,
  not on this row.
- **Deconvolution error.** Per-lineage estimates carry posterior uncertainty
  that WGCNA does not see; `bmind-047` uses the posterior mean and discards
  the SE. A noisier matrix is harder to find structure in, so this test is
  conservative - it can hide an effect but is unlikely to invent one.
- **The global thinning is still there.** This tests module preservation, not
  overall graph density, and the 11-22% sparser lithium graph reported in
  [[2026-11-10-lithium-graph-is-sparser-in-three-independent-draws]] is not
  addressed by it.

## Two things had to be fixed before the row would run at all

- **52 of 100 draws were being rejected.** `modulePreservation()` refuses an
  entire call if any gene in any set fails `goodGenes()`, and a bMIND
  per-lineage profile has a few near-constant genes in most 74-sample
  subsamples - the Gibbs posterior for a gene the sampler cannot resolve in
  that cell type collapses, and which genes those are depends on the draw.
  gran flagged none and finished clean; mono and B lost 16 draws each; T lost
  its nolith run and all 20. Fixed by trimming the flagged genes from the
  test matrix only - 1 to 4 genes of 12,368, always the same handful
  (ENSG00000112053 and ENSG00000198216 in T, ENSG00000154898 and
  ENSG00000147036 in B). Every control matrix flags zero, so the reference
  modules never moved.
- **The row would otherwise have reported COMPLETE having done nothing.**
  `normalise_bmind()` tested `is.list()` first and returned the profile
  object unchanged, so the five lineage names were never found, the loop body
  never ran, and `all(logical(0))` is `TRUE`. A silent success, not a failure.

## Provenance

The profiles are the 2026-10-01 bMIND run, not a rebuild. `bmind-046` was
queued to redraw them for 3.1 h; bMIND is a Gibbs sampler, so a second run
returns a different posterior draw, and every DE-arm result downstream was
computed from the first. It adopts them instead and records that it did.

See also [[2026-11-09-perm031-at-157-draws-is-a-clean-negative]] and
[[2026-11-09-eosinophil-module-survives-composition]].
