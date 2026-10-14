# 003 — LM22 circularity check

**The concern.** Analysis 002 found PC2 of the module eigengenes is lithium-associated and
73.3% explained by the ILR composition balances. Those balances derive from CIBERSORTx
fractions, which are estimated **from the LM22 marker genes in this same expression
matrix**. The composition estimate and the expression PC therefore share genes, so the
73.3% could be partly circular. `IMAC_README` catalogue item 7 names exactly this check.

**The test.** Drop the 350 LM22 markers present among the 12,368 genes (2.83%, spread
across 30 of 33 modules — not concentrated). Recompute eigengenes from the *same* module
assignments, redo the PCA, and re-ask both questions. The lithium-associated dimension is
identified by the same a-priori rule as 002 (test the Li & Ji effective dimensions, BH
across them) — no post-hoc PC selection.

**Result.**

| gene set | genes | Meff | PC | β | p | BH | attenuation | PC r² composition |
|:--|--:|--:|:--|--:|--:|--:|--:|--:|
| all genes | 12,368 | 6 | PC2 | −0.9081 | 4.8e-4 | 0.0029 | 68.1% | 73.3% |
| **LM22 excluded** | **12,018** | 6 | PC2 | **−0.8996** | 5.0e-4 | **0.0030** | **67.4%** | **72.2%** |

Correlation between the two gene sets' components: **PC1–PC4 all 1.000**.

**Verdict: finding 002 is not explained by the LM22 marker genes.** The lithium
association survives (BH 0.0030) and the composition share holds (72.2% vs 73.3%).
Removing every gene used to estimate composition moves the result by about one percent.

**A second observation.** That the PCs are *identical* (r = 1.000) after removing 2.83% of
genes — selected to be the most composition-informative ones available — says the
structure is extremely diffuse. No small gene subset is pivotal, which is consistent with
a network whose mean connectivity is 8.3% of the genome. It also means gene-level
exclusion experiments on this data will generally be uninformative: the signal is not
carried by identifiable gene subsets.

**What this does not test.** Circularity through a different route: CIBERSORTx fractions
are estimated from the whole expression profile, not only from marker genes, so some
shared variance is unavoidable by construction. This check removes the most direct path,
not every path.

**Files.** `config.json`, `results/pc_lithium_by_geneset.csv`, `results/side_by_side.csv`,
`results/summary.rds`.
