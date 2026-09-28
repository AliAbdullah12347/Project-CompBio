# What is and is not reproducible

Boltz et al. analysed **1,730 whole-blood samples with genotypes** (dbGaP **phs002856.v1**,
controlled-access; approval needs a faculty-PI Data Access Request). None of that is public. This
recreation applies their methods to **GSE124326** (Krebs et al. 2020): 444 post-QC samples from the
same Utrecht cohort, sequenced by the same lab, with deposited LM22 cell fractions. The two cohorts
overlap in origin but not in size, depth or processing, so **the target is the same method and the
same qualitative conclusions, not the same counts.** `results/validation_scorecard.csv` scores
every claim on that basis.

---

## 1. Not reproducible — skipped

| Paper analysis | Needs | Why it cannot be done here |
|:--|:--|:--|
| cis-eQTL mapping, QTLtools, 8 cell types + bulk (Fig 2, Table 1 eGene counts) | genotypes | GSE124326 has none; phs002856 is controlled |
| Decon-eQTL comparison (Table S1) | genotypes | as above |
| Storey π1 replication vs bulk / eQTL Catalogue / OneK1K (Table S2, Fig S5) | eQTLs | as above |
| FUSION TWAS, heritability, colocalisation, HTR6 conditioning (Tables 2–4, S4–S12, Fig 3) | genotypes + GWAS | as above; Tables S4–S12 are kept in `data/supplementary/` for reference |
| lithium × SNP interaction eQTL (Tables S13–S14, Fig 4A–B) | genotypes | as above |
| CIBERSORTx vs clinical blood counts, n = 143 (Fig S1) | CBC records | not deposited |
| UCLA ATLAS neutrophil validation (Fig S7) | EHR | private biobank |
| Table S3 **BLUEPRINT** rows | BLUEPRINT median TPM | not downloaded (large); DICE/Schmiedel rows are done |

## 2. Reproduced — with these deviations

| # | Paper | Here | Effect |
|:--|:--|:--|:--|
| D1 | kallisto TPM, GENCODE v33lift37, 5.9 M reads/sample | TPM rebuilt from HTSeq union counts ÷ GENCODE v19 union-exon length, ~8 M reads | Gene-level approximation of kallisto TPM. Small RNAs (RN7SL1/2, RN7SK) take ~60% of TPM mass in these libraries, as they would under kallisto, which lowers every other gene's TPM |
| D2 | CIBERSORTx with batch correction | CIBERSORT LM22 fractions deposited in GSE124326 | No CIBERSORTx account needed. Monocytes are ~4× higher here (0.23 vs 0.05) and memory B is ~0: LM22 fits whole blood poorly without batch correction |
| D3 | "8 cell types with proportions > 0.02" — naive B, memory B, CD8, naive CD4, memory resting CD4, resting NK, monocytes, neutrophils | The same rule selects a different eight: drops **memory B** (mean 0.00001) and **CD8** (0.007); adds activated memory CD4 and resting mast | bMIND cannot estimate a cell type whose fraction is ~0. Regressions report Boltz's eight *and* the rule's eight; bMIND uses the rule's |
| D4 | covariates age, sex, **RNA concentration**, RIN | age, sex, RIN | RNA concentration is not in the GEO deposit |
| D5 | cases = BP + SCZ | cases = BP-I + BP-II (239) | cohort has no SCZ |
| D6 | "at least 1 TPM in at least 436 individuals (about 25%)" | ≥ 1 TPM in ≥ 25% (111 of 444) → 18,006 genes | rule carried as a proportion |
| D7 | "log2-transformed" | log2(TPM + 1) | pseudocount unstated; +1 is what bMIND itself uses |
| D8 | "first 50 expression PCs" | `prcomp` on genes centred, not scaled, all 444 samples, computed once | scaling and sample set unstated; sentence order implies once, before subsetting |
| D9 | bMIND "with flag np = TRUE" | bMIND() with no prior supplied (its data-driven-prior mode); `np = TRUE` in bmind_de() | MIND 0.3.3, the only public release, has `np` only in `bmind_de()`. See §3 |
| D10 | bMIND on all genes | real bMIND MCMC on a fixed random 2,000 genes (seed 20240201) | full run ≈ 13 h on this laptop; every use is a cross-gene correlation or PCA |
| D11 | bmind_de via MCMC | exact posterior for all 18,006 genes; real bmind_de() on 400 genes to confirm agreement | see §3 |
| D12 | bmind_de gene set: 19,099 genes (from Table S15; not described) | the 18,006 bulk-filtered genes | the paper does not say how 19,099 were chosen |
| D13 | preprint prediction: metric not stated; Supplementary Table 2 is an image | test-set AUC; qualitative comparison only | the table could not be retrieved (bioRxiv bot wall) |

## 3. bmind_de: why the exact posterior is the same test

`bmind_de(np = TRUE)` fits, gene by gene, the fixed-effects model

```
x = PCs·γ + Σ_k W_k·co·β_k,co + Σ_k W_k·ca·β_k,ca + ε
```

by MCMCglmm with its default priors: flat on coefficients (`V = 1e10·I`) and inverse-Wishart
`(V = 1, ν = 0)` — i.e. p(σ²) ∝ 1/σ² — on the residual. Under exactly those priors the posterior of
each contrast d_k = β_k,ca − β_k,co is Student-t with n − p degrees of freedom centred on the OLS
estimate, so bmind_de's p-value 2·min(P(d>0), P(d<0)) equals the classical two-sided t-test p. The
MCMC approximates that limit with Monte Carlo error and a 1/nsamp floor, which bmind_de lifts by
re-running strong genes with up to 10⁶ draws. `src/04_bmind.R` computes the limit exactly for all
genes in seconds, and runs the real `bmind_de()` on 400 genes; `results/bmind_de_mcmc_vs_exact.csv`
records their agreement.

## 4. Problems found in the paper itself

1. **Case/control logFC range is misreported.** Text: −0.126 to 0.104. Deposited Table S15: −0.217
   to 0.104; three DEGs lie below −0.126 (lowest PDK4, −0.217).
2. **The Krebs comparison list is unstated, and it is not a deposited one.** Boltz's 100 lithium DEGs
   overlap Krebs' *documented* lithium model (recreated in `../Krebs Paper`, ~1,000 DEGs; Krebs'
   text says 976) at 32 genes, OR 5.92, p 1.1e-12 — matching the reported 33, OR 6.43,
   p 4.7e-14. Against Krebs' deposited pre-cell-type sheet the overlap is 35 but OR only 1.8;
   against the post-correction sheet, 26 at OR 22.6. So Boltz compared against the list in Krebs'
   text, which Krebs never deposited.
3. **The bMIND DE gene universe (19,099) is not the bulk DE universe (17,194)** and is not described.
4. **The lithium-prediction analysis was removed between preprint and publication** with its
   conclusion ("gene expression does not provide additional predictive value over the cell type
   proportions"). It rested on one random 70/30 split, with gene selection and PC residualisation
   done on all samples before splitting. It is recreated in `src/06_preprint_prediction.R`,
   together with the repeated-split check it lacked.

## 5. Findings specific to this cohort

- **Case/control DE is dominated by collection batch.** Every case is in assessment group A; the
  controls span A and B. Boltz's model adjusts only for 50 PCs. Here it gives 236 DEGs; adding
  assessment group drops that to 80, and group A alone gives 33. Boltz's cohort may not share this
  structure, but any case/control analysis of GSE124326 must adjust for it.
- **Lithium non-users vs controls:** none of Boltz's eight types differ (as they report), but
  activated memory CD4 T cells — a type only the > 0.02 rule selects here — do (Bonferroni
  p = 0.019). The controls are QC-filtered non-randomly (35 of 240 dropped vs 0 of 226 BP-I), so this
  is flagged, not interpreted.
