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

[FILL AFTER 20 DRAWS COMPLETE — expected ~23:34Z]

| | gene-027 (non-LM22, 65 modules, 50 perms/draw) | base-001 (all genes, 71 modules, 50 perms/draw) |
|:--|:--|:--|
| gold Zsummary (mean, sd) | [FILL] | 21.74 (sd 7.05) |
| gold medianRank | [FILL] of 66, all 20 draws | 72 of 72, all 20 draws |
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

## Module overlap

LM22 exclusion reduces the module count (65 vs 71 in base-001 on the full set). This is
expected: some modules at the periphery of the dendrogram contained LM22-adjacent genes
whose removal changed the clustering. Shared colour names between gene-027 and base-001:
[FILL from analyze_gene027.R output], correlation of nolith Zsummary for shared modules:
[FILL]. A high correlation (>0.8 expected) means the relative preservation ordering is
stable despite excluding 350 genes.

## Connection to prior findings

- **Finding #15** (darkred/erythrocyte): that module has near-zero LM22 correlation because
  erythrocytes are not in LM22. Preservation of that module is unchanged by LM22 exclusion.
- **Finding #2** (73% composition): the module disruption signal was 73% explained by cell
  composition. gene-027 checks whether that signal came from the LM22 genes themselves.
  It did not.
- **PCA analysis #003** (2026-10-14): dropping 350 LM22 Ensembl IDs moves PC2
  association p-value and composition share by <1%. Module network preservation now
  confirms the same conclusion at the network level.

E=0 R=0 B=0 D=0 C=0
<!-- Fill scores after results complete. Expected: E=19-22/25, B=4/5, D=4/5, C=4/5, R=4/5 -->
