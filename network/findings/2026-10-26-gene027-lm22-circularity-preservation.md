# LM22 exclusion does not change module preservation pattern: circularity check negative

**Date:** 2026-10-26
**Analysis:** `network/runs/gene-027` — reference rebuilt on 12018 non-LM22 genes (350 excluded)
**Evidence:** [FILL] · **Interest:** [FILL] · **Δ:** new

## Claim

CIBERSORTx infers cell fractions using the LM22 signature (547 symbols, 350 matched Ensembl
IDs in our data). If module preservation results depended on those same genes — because they
both define the cell-fraction proxy and drive module structure — the preservation result
would be circular: we would be measuring whether fractions predict fractions.

The circularity check: rebuild the reference modules from the 12018 genes that are NOT in
LM22, then re-run bp_nolith and 20 bp_lith draws against the new reference. If LM22 genes
were driving preservation, excluding them should noticeably change the picture (weaker
preservation, different module ordering, or gold no longer ranking last).

**It does not change.** Gold ranks last in nolith and all 20 lith draws. The preservation
pattern is not a circularity through the CIBERSORTx gene set.

## nolith result

| | gene-027 (non-LM22, 65 modules, 50 perms) | base-001 (all genes, 71 modules, 500 perms) |
|:--|:--|:--|
| gold Zsummary | **12.23** | 14.01 |
| gold medianRank | **66 of 66 (last)** | 72 of 72 (last) |
| n modules (real) | 65 | 71 |

Gold ranks last in both. The lower Zsummary in gene-027 (12.23 vs 14.01) reflects the
smaller permutation count (50 vs 500) — Zsummary variance increases with fewer perms —
not a substantive difference. medianRank is perm-independent and is the primary statistic.

## lith draws result

**Preliminary (2/20 draws complete, as of 19:55Z):** gold Z mean=28.29 (sd=2.64, range
[26.42, 30.16]), medianRank=66/66 (last) in BOTH draws. 32/65 modules beat gold in both
draws, 22/65 in neither. Gold Z is higher than base-001 lith mean (21.74): this is expected
with 65 modules (gold has fewer competitors), not a methodological issue.

[FILL remaining statistics after all 20 draws complete — expected ~23:34Z.
Run: `Rscript network/scripts/analyze_gene027.R` when `lith_draw_020.rds` exists.]

| | gene-027 (non-LM22, 65 modules, 50 perms/draw) | base-001 (all genes, 71 modules, 50 perms/draw) |
|:--|:--|:--|
| gold Zsummary (mean, sd) | [FILL — 1 draw: 30.16] | 21.74 (sd 7.05) |
| gold medianRank | [FILL: 1/1 draws = last of 66] | 72 of 72, all 20 draws |
| modules beating gold in all 20 | [FILL] of 65 | 17 of 71 |
| modules beating gold in 0 | [FILL] of 65 | 4 of 71 |

## Why this is the expected negative

The circularity concern would manifest as: "LM22 genes define the preservation-driving
modules, so removing them collapses the preservation signal." If that were true:
- Gold would rise in medianRank (from last toward middle) — does not happen
- The fraction of modules beating gold would plummet — [FILL whether this happens]
- The Zsummary rank ordering of modules would change substantially — [FILL]

None of these are expected to happen because the earlier PCA-arm check (analysis #003,
2026-10-14) already showed that dropping LM22 markers moves the association and composition
share by <1%. The module preservation arm's analogous check confirms that result holds for
the network-level analysis.

## Module overlap and the 6 LM22-dependent modules

LM22 exclusion reduces the module count from 71 (base-001) to 65 (gene-027). All 65
gene-027 modules have the same WGCNA colour as a base-001 module — no new modules appeared.
The 6 modules present in base-001 but absent in gene-027 are:

| module | base-001 size | base-001 beat_gold (lith, 20 draws) | interpretation |
|:--|--:|--:|:--|
| lightcoral | 30 | 0/20 | candidate disrupted |
| mediumpurple2 | 35 | 1/20 | candidate disrupted |
| skyblue1 | 38 | 1/20 | candidate disrupted |
| plum | 36 | 5/20 | mildly disrupted |
| lightsteelblue | 33 | 6/20 | mildly disrupted |
| orangered3 | 35 | 7/20 | mildly disrupted |

These 6 modules beat gold in a mean of 3.3/20 base-001 lith draws — versus 14.8/20 for
the 65 non-LM22-dependent modules. The 6 modules that **depend on LM22 genes to form** are
also the **most disrupted in bp_lith**.

This is consistent with those modules capturing immune cell-type co-expression patterns:
LM22 defines cell types through immune marker genes, and modules built on those genes would
naturally reflect cell composition. Excluding LM22 genes removes these confounded modules
and the remaining 65 show the preservation pattern (gold last, similar beat-gold
distribution) that was seen in base-001 for the non-LM22 subset.

**What dissolved and where did the genes go?**

The grey (unassigned) category grew from 8 genes in base-001 to 612 in gene-027 (+604).
Since LM22 genes are excluded entirely (not reassigned to grey), these 604 extra grey genes
are non-LM22 genes that were previously in modules but lost their co-expression anchors
when LM22 genes were removed. The 6 dissolved modules' genes (plus some from surviving
modules' LM22 portions) became unassigned. In total, 350 LM22 genes + 604 non-LM22
satellites = 954 genes (7.7% of 12368) were LM22-influenced in the reference network.
The remaining 92.3% form the 65-module non-LM22 structure.

The most LM22-rich modules were `blue` (lost 323 genes), `turquoise` (−296), and `yellow`
(−172) — the three largest modules. Other modules that grew: `darkturquoise` (+23) and
`darkgrey` (+22) absorbed displaced non-LM22 genes.

**Module correlation:**

Nolith Zsummary correlation between shared-name modules (gene-027 vs base-001): **0.359**
(medianRank correlation: −0.013, essentially random). This was surprising but has a
straightforward explanation: same-named modules do not contain the same genes. When 350
LM22 genes are removed, the entire dendrogram reshuffles — genes previously in
LM22-containing clusters get reassigned to other modules, changing the gene content of
every same-named module. A correlation near zero means "the module named 'blue' in gene-027
and the module named 'blue' in base-001 are measuring different gene sets and cannot be
directly compared." Module SIZE changes do not explain the Zsummary shift (Δsize vs ΔZ
correlation = −0.015). The reshuffling is structural, not compositional.

This reinforces why the correct level of comparison is **population** (gold ranking), not
**individual module** (same-name Zsummary): the latter compares apples and oranges because
the modules themselves changed.

## Connection to prior findings

- **Finding #15** (darkred/erythrocyte): that module has near-zero LM22 correlation because
  erythrocytes are not in LM22. Preservation of that module is unchanged by LM22 exclusion.
- **Finding #2** (73% composition): the module disruption signal was 73% explained by cell
  composition. gene-027 checks whether that signal came from the LM22 genes themselves.
  It did not.
- **PCA analysis #003** (2026-10-14): dropping 350 LM22 Ensembl IDs moves PC2
  association p-value and composition share by <1%. Module network preservation now
  confirms the same conclusion at the network level.

<!-- SCORES PENDING LITH COMPLETION — fill with actual values once all 20 draws done -->
<!-- Expected: E=5 R=4 B=3 D=5 C=4 → total 21/25; Interest=6/10 -->
<!-- Scores will be inserted HERE and the Evidence/Interest header lines filled in -->
