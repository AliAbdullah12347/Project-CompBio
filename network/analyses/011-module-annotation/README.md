# 005 — Module annotation, sample structure, hub genes

Three catalogue items at the usable decomposition (33 modules), all from cached data.

## 16. Sample structure — the important one

Eigengenes on all 474 labelled samples, PCA, one-way ANOVA of each PC on group:

| PC | var | F | p (BH) | control | bp_nolith | bp_lith |
|:--|--:|--:|--:|--:|--:|--:|
| **PC1** | 77.5% | 21.41 | **5.2e-09** | −1.43 | +0.98 | +1.73 |
| **PC2** | 10.1% | 9.54 | **1.8e-04** | +0.22 | +0.36 | **−0.52** |
| PC3 | 6.1% | 2.46 | 0.115 | | | |

**PC1 is a diagnosis axis** — controls on one side, both bipolar groups on the other, with
an ordering control < bp_nolith < bp_lith.

**PC2 is a lithium axis** — control and bp_nolith sit together (+0.22, +0.36) and bp_lith
separates (−0.52). That is the *same component* analyses 002 and 004 identified as
lithium-associated, found independently here on the **full 474-sample cohort** rather than
BP1 only. Two different sample sets, same axis.

This matters because 002/004 tested lithium within BP1, where a sceptic could ask whether
PC2 was simply the within-BP1 contrast by construction. Here PC2 is defined over all three
groups and still places bp_lith apart from both control and bp_nolith.

## 15. Module / DEG overlap (Fisher exact, BH across 33 modules)

**Cell-type lithium DEGs (120 genes mapped): 6 of 33 modules enriched.**

| module | genes | overlap | expected | OR | BH |
|:--|--:|--:|--:|--:|--:|
| greenyellow | 438 | 15 | 4.2 | 3.99 | 5.8e-04 |
| saddlebrown | 64 | 6 | 0.6 | **11.05** | 5.8e-04 |
| **green** | 835 | **0** | 8.1 | **0.00** | 4.0e-03 |
| tan | 407 | 12 | 3.9 | 3.33 | 4.8e-03 |
| darkorange | 90 | 5 | 0.9 | 6.22 | 9.9e-03 |
| magenta | 460 | 12 | 4.5 | 2.93 | 9.9e-03 |

The lithium DE signal is **not** spread uniformly across the network — it concentrates in
specific modules (saddlebrown OR 11.05) and is **significantly depleted** from one
(green: 0 observed against 8.1 expected across 835 genes).

**LM22 markers (350 genes): 6 of 33 enriched**, led by darkgreen (OR 4.70) and turquoise
(OR 1.90); purple (OR 0.15) and cyan (OR 0.11) are depleted. So cell-type marker genes do
concentrate in particular modules — consistent with analysis 001's finding that only one
module correlates strongly with a lineage fraction, since concentration of markers and
correlation with an estimated fraction are different things.

## 11. Hub genes (kME, control reference)

Top hub per module, largest modules first: turquoise *ARFGAP2* (0.961), blue *RBMX*
(0.961), brown *TRIM4* (0.959), yellow *SIKE1* (0.979), green *TRMT44* (0.953), red
*TTBK2* (0.969), black *NRD1* (0.976), pink *CUL9* (0.970).

Every top hub has kME between 0.95 and 0.98 — a suspiciously uniform range. In a sparser
network hub kME varies much more by module. This is the density showing up again: when
mean connectivity is 8.3% of the genome, every module has a gene almost perfectly
correlated with its eigengene, so "hub" carries less information here than it normally
would. Hub identity should not be over-interpreted on this data.

**Files.** `config.json`, `results/module_vs_degs_fisher.csv`,
`results/sample_structure_anova.csv`, `results/hub_genes_kme.csv`, `results/summary.rds`.
