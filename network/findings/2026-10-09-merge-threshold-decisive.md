# The merge threshold is decisive: 33 modules vs 3

**Date:** 2026-10-09
**Analysis:** `network/runs/mod-000-resolution`
**Evidence:** 22/25 · **Interest:** 7/10 · **Δ:** -1.8

## Claim
A 2x5 grid at n=234 shows `cutreeDynamic` finds 18 modules at deepSplit 2 and 33 at
deepSplit 4; `mergeCutHeight` then determines everything. At 0.25 (base-001's setting):
3 modules, largest 97.5% of genes. At 0: 33 modules, largest 11.68%, median size 248.
Four of ten cells meet the usability target.

## Mechanism
Median eigengene pair correlation is 0.739-0.765 against a merge threshold of 0.75, so
over half of all module pairs qualify on the first pass and the merge cascades.

## Consequence
**Recommended for all downstream rows: deepSplit=4, mergeCutHeight=0.** Everything
base-001 measured was computed on the degenerate cell and must be re-run.

E=5 R=5 B=3 D=5 C=4. Interest 7: confirms a hypothesis rather than overturning one.
