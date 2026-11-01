# ann-042: module overlap with LM22 markers and cell-type lithium DEGs

**Date:** 2026-11-01
**Analysis:** `network/runs/ann-042` — Fisher exact test of 71 modules against 5 gene sets (LM22 signature genes, and cell-type lithium DEG sets T/mono/gran/NK), BH across modules within each set
**Evidence:** 17/25 (E=3 R=5 B=3 D=3 C=3) · **Interest:** 5/10 · **Δ:** −1.4

## Claim
Do the modules that carry lithium-DEG signal coincide with cell-type marker modules?

## Results (background N=12,368; LM22 350 of 547 symbols matched)
14 module-set overlaps significant at BH 5%.
- **LM22_all (350 genes):** 8 modules. Strongest: thistle1 (22 overlap vs 1.9 expected, OR 17), lavenderblush3 (20 vs 1.7, OR 18), mediumpurple2 (11 vs 1.0, OR 16), plus brown4, floralwhite, thistle2, antiquewhite4, lightyellow.
- **CT_LI_gran (109):** yellow, darkred, lightyellow. **CT_LI_T (28):** mediumpurple2 (4 vs 0.08), yellow. **CT_LI_mono (20):** yellow. **CT_LI_NK (4):** untestable.
- **Key overlap:** mediumpurple2 and lightyellow are both LM22-marker modules and cell-type lithium-DEG modules; yellow and darkred carry the lithium DEG signal in several cell-type sets and in the ann-041 whole-blood list, but only lightyellow/mediumpurple2 are also LM22-enriched.
Full table: `results_geneset_overlap.csv`.

## Interpretation and limits
- Lithium-DEG-enriched modules partly coincide with cell-type marker modules, which fits the composition account (finding #2). It is only partial: yellow and darkred are not LM22-enriched here, so those are not explained by this test.
- The CT_LI sets are small (4–109 genes), so ORs for them are unstable (e.g. mediumpurple2 OR 66 from 4 genes).
- BH is within set, not pooled across the five sets. The CT_LI sets are nested subsets of related lists, so they are not independent.
- gene-027 already showed the preservation results do not depend on LM22 genes; this does not contradict it.

<!-- SCORES: E=3 R=5 B=3 D=3 C=3 -->
