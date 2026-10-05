# Superseded — do not cite these numbers

Everything in this folder came from `bmind_de()`, which was abandoned after two
defects were verified by direct test. The reasoning is in `../METHODS.md` §9.

**The defects**

1. **A hard p-value floor at 2×10⁻⁴.** Only ~600 distinct p-values appeared
   across 61,840 tests. Benjamini–Hochberg at FDR 0.05 over 12,368 genes needs
   the smallest p-value below 4.0×10⁻⁶. The floor is 50× too large, so *no* gene
   can reach significance however strong the true effect. Raising MCMC iterations
   tenfold moved the floor from 0.086 to 0.073 — it is structural, not a
   sampling-resolution problem.

2. **Covariates silently discarded.** `bmind_de(covariate = NULL)`,
   `bmind_de(covariate = CV)` and `bmind_de(covariate_bulk = CV)` return
   byte-identical p-values. Every result here is unadjusted for age, sex,
   tobacco, RIN, plate, sequencing PCs and the cell-composition balances.

Any "zero differentially expressed genes" in these files is therefore
**arithmetically forced, not observed**, and means nothing.

**What replaced them**

| superseded | replacement |
|:--|:--|
| `02_bmind.R` | `../scripts/06_bmind_profiles.R` — `bMIND()` without the phenotype |
| `05_celltype_de.R` | `../scripts/07_celltype_limma.R` — limma in a design we control |
| `celltype_*.csv`, `celltype_pvalues.rds` | `../results/celltype_limma_*.csv` |
| `bmind_*.rds` | `../results/bmind_profiles.rds` |

Kept rather than deleted because the defects are a reportable finding in their
own right: Boltz et al. used `bmind_de` and reported zero cell-type lithium
genes across 1,730 samples. These files are the evidence for that claim.
