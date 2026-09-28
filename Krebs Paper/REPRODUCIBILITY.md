# What can and cannot be reproduced from public data

Krebs et al. (2020), *Psychological Medicine* 50(15), 2575–2586.
GEO accession **GSE124326**.

The paper deposited its **gene expression data only**. Several analyses in it draw on inputs
that were never deposited and are not obtainable from GEO. This file records exactly which,
so that nothing in `src/` pretends to a reproduction it cannot perform. No code is shipped
for a blocked analysis — a script that cannot run is worse than an honest gap.

---

## Reproducible, and implemented here

| Paper result | Script | Status |
|:--|:--|:--|
| Sample QC to 444 subjects | `00_prepare_data.py` | reproduced |
| Gene filter (>10 counts in ≥90%) | `00_prepare_data.py` | reproduced — **12,353 genes, matching their DE table exactly** |
| Table 1 demographics | `00_prepare_data.py` | reproduced (one-subject deposit/paper discrepancy, below) |
| BD case/control DE — 6 DEGs | `02_differential_expression.R` | reproduced — same 6 genes |
| Lithium DE | `02_differential_expression.R` | see "Unresolved" below |
| Lithium DE + cell-type covariates | `02_differential_expression.R` | reproduced — r(logFC) = 0.989 vs their table |
| WGCNA network and modules | `03_wgcna.R` | see results |
| CIBERSORT proportions → lithium (stepAIC) | `04_celltype.R` | see results |
| Cell-type vs reference intervals | `04_celltype.R` | reproduced |

The deposited CIBERSORT LM22 fractions are in the series matrix, so CIBERSORT itself does
**not** need re-running — we use the authors' own output, which is what they analysed.

---

## Blocked: input never deposited

### 1. Polygenic risk score analysis
**Needs:** individual-level genotypes.
The supplementary states genotyping was done on 234 cases and 187 controls with the Illumina
Infinium OmniExpressExome, imputed to 6,828,668 SNPs. **None of this is in GEO** — the
accession contains expression only, and no dbGaP/EGA accession is given for the genotypes.
Without individual genotypes a PRS cannot be computed at all, so the PRS-vs-diagnosis t-test
(t = −3.42, p = 6.88 × 10⁻⁴) and the PRS-as-DE-trait analysis are **not reproducible by anyone
outside the original group**.

### 2. MAGMA gene-set analysis
**Needs:** GWAS summary statistics for three traits, plus a 1000 Genomes LD reference panel.
- BD (Stahl et al. 2019) — publicly available from PGC.
- Schizophrenia (PGC 2014) — publicly available.
- **Self-reported depression (Hyde et al. 2016) — 23andMe. Not publicly available**; it requires
  a data transfer agreement with 23andMe.

So at most two of the three traits could ever be reproduced, and only after ~2–3 GB of
downloads. Figure S9 and Table S8 are therefore **partially blocked**. The same applies to
the secondary sleep-trait analysis (Jones et al. 2016) and to the LD Score Regression
heritability/genetic-correlation estimates in Table S7.

### 3. BMI sensitivity analysis
**Needs:** per-subject BMI.
The supplementary reports BMI was available for 380 of 480 subjects and was used for an ANOVA,
a stepAIC on cell types, and a controls-only DE. **BMI is not among the 35 characteristics
fields in the series matrix.** Figures S2 and S3 and Table S2 cannot be regenerated.

This is worth noting for our own project: it means the E-value approach to bounding body-fat
confounding is not a stylistic choice but the only option available.

### 4. LM22 signature gene list
**Needs:** the LM22 matrix distributed with CIBERSORT.
The 547-gene LM22 signature is available only after registering at cibersort.stanford.edu; it
is not in the GEO deposit and is not redistributable. This blocks the exact form of two checks:
- "16 of 60 neutrophil signature genes were also lithium DEGs" (OR 4.64, p = 4.45 × 10⁻⁶)
- the module × cell-type hypergeometric and linear-model enrichment (Figure 3), which needs the
  331 LM22 signature genes expressed in this sample.

The *deposited cell fractions* are unaffected — only the gene-level signature membership is.

### 5. DAVID functional annotation
**Needs:** the DAVID web service (Figure 1c, module annotation in Table 2).
DAVID is an interactive web tool whose results depend on its database version at query time,
so the 2019 enrichment scores are not reproducible even with access. Any modern re-run would
use a different annotation release and give different numbers.

### 6. Upstream processing
Reads were aligned with TopHat2 to hg19 and quantified with HTSeq union mode. **FASTQ files
were not deposited** — GEO holds the count matrix only. Everything upstream of counts
(alignment rate 96.0%, 33.9% duplicates, the 18 Picard sequencing metrics) is taken as given.
The three sequencing-metric PCs *are* deposited, so the covariate model is unaffected.

---

## Resolved: the deposited lithium sheet is mislabelled, but the paper's analysis is sound

The paper reports **976** lithium DEGs in its Results, **897** in the MAGMA methods section,
while the deposited `Li_DEGs_PreCellTypeCorrection` sheet contains **3,031** at FDR < 0.05.
We can now account for this.

**The deposited sheet is not the lithium main effect.** It is the contrast *bipolar-on-lithium
versus control* — the case and lithium coefficients summed. Independently confirmed twice, by
`contrasts.fit` and by an equivalent three-level reparameterisation agreeing to 2 × 10⁻¹⁴:

```
r(logFC) = 0.9993     sign concordance 99.0%
n significant 3,001 vs their 3,031     overlap 2,941     Jaccard 0.951
```

Regressing the deposited logFC on all 15 model coefficients with free weights returns
R² = 0.9990 with weights +1.014 on `case` and +1.003 on `lithium`, everything else below
0.026 — the sheet *is* that sum, not something resembling it.

**But the paper's own text describes a correctly specified lithium main effect**, and our model
reproduces it on every statistic the paper reports:

| | Paper's 976 | Our lithium main effect | Deposited sheet |
|:--|--:|--:|--:|
| n at FDR < 0.05 | 976 | 1,023 | 3,031 |
| \|FC\| mean | 0.20 | **0.195** | 0.165 |
| \|FC\| max | 0.82 | **0.815** | 1.026 |
| \|FC\| sd | 0.10 | **0.100** | 0.086 |
| up-regulated | 77.3% | 71.2% | 48.9% |
| carry-over to post-correction (paper: 194/233, 83.2%) | 83.2% | **197 (84.5%)** | 181 (77.7%) |
| Anand overlap genes significant | — | **20/20** | 16/20 |
| Breen overlap genes significant | — | **130/134** | 91/134 |

The max \|FC\| is decisive on its own: the deposited sheet reaches 1.026, so the paper's
976-gene set — capped at 0.82 — cannot be any top-N subset of it.

**Conclusion.** The error is confined to the deposited supplementary sheet, which holds
"BP-on-lithium vs control" where it should hold the lithium main effect. The paper's reported
analysis is reproducible; its deposited file is not the analysis it describes. The residual
1,023-vs-976 gap is consistent with the one-subject deposit/publication discrepancy below.

The 897 in the MAGMA section remains unexplained; no threshold on either sheet yields exactly
976 or 897. Anyone reusing "the 976 lithium DEGs" must say which list is meant — and should
not use the deposited pre-correction sheet for that purpose.

Full search in `src/02c_presheet_resolution.R` (81 specifications,
`results/presheet_model_sweep.csv`); independent confirmation in `src/verify_presheet_claim.R`.

> Note on `02c`: its sweep scores contrasts with `contrasts.fit`, which is approximate for
> multi-coefficient contrasts under voom's gene-specific weights. logFC and all correlations
> are unaffected; the significant-gene count is slightly conservative (2,966 rather than the
> exact 3,001). `verify_presheet_claim.R` carries the exact figures.

---

## Deposit vs publication: one-subject discrepancy

| | Paper | Deposit |
|:--|--:|--:|
| Cases post-QC | 240 | 239 |
| Controls post-QC | 204 | 205 |
| BP1 / BP2 | 227 / 13 | 226 / 13 |
| Non-lithium cases | 88 | 87 |
| Lithium users | 152 | 152 |

One subject is classified as a case in the paper and as a control in the deposit. The lithium
user count is unaffected. We report the deposit, because that is what the code can analyse.
