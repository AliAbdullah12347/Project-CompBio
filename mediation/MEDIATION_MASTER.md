# Mediation arm — master document

Project 7, BIOL/MATH 481, Colgate · Ali Abdullah & Ava Gillespie · generated 19 October 2026 by `scripts/15_master_doc.R` from the saved results — every number below is read from a results file.

**How to read this.** The first section is a brief overview of everything. Each later section goes into detail, in the order the work was done. Section 16 lists every flagged gene by name. Section 17 lists the files, so any number here can be traced to the script that produced it.

## Overview — everything we did, in brief

### The question and the design

- **Question:** how much of lithium's effect on whole-blood gene expression runs *through the cell mixture* rather than directly?
- **Sample:** the 226 bipolar I subjects of GSE124326 (152 lithium users, 74 non-users). Lithium varies only within bipolar I.
- **Exposure** lithium (0/1) · **mediator** the cell mix, as ILR balances · **outcome** gene expression, one model per gene (12,018 genes; the 350 LM22 signature genes held out) · **covariates** age, sex, tobacco, RIN, plate, sequencing PCs 1–3.
- **Two versions of the mediator, both run end to end:** the pipeline's tree (5 lineages → 4 balances) and the proposal's tree (6 cell types → 5 balances).
- **Protocol fixed in writing before any estimate** (`PRESPEC.md`), including two later addenda written before their analyses ran.

### How it was estimated and checked

- **Estimator:** natural direct and indirect effects (NDE, NIE) for a vector mediator, allowing a lithium × cell-mix interaction (VanderWeele & Vansteelandt 2014) — chosen over product-of-coefficients because Prof. Ay's feedback required the interaction to be tested, not assumed.
- **Verification:** 27 checks, 0 failures, against brute-force Monte-Carlo counterfactuals, independent `lm()` fits, a 4,000-replicate bootstrap and 500-dataset interval coverage.
- **External validation:** the mediation *total* effect reproduces the separate differential-expression arm (per-gene correlation r = 0.9985; 1,224 of the DE arm's 1,376 lithium genes recovered; Jaccard 0.860).
- **Independent adversarial review** of the statistics found four real defects after the checks had passed; all were fixed and everything was re-run. A fifth claim was tested and rejected.

### Main results

- **Lithium changes the expression of 1,272 genes** (total effect, BH 5%) — 888 up, 384 down.
- **Lithium shifts the cell mix:** granulocytes rise relative to monocytes (p = 0.012) and myeloid relative to lymphoid (p = 0.047); nothing moves inside the lymphoid compartment.
- **About 40% of lithium's effect runs through the cell mix** — median proportion mediated 0.396 (95% CI 0.152–0.567); 0.449 on the proposal tree (0.236–0.630).
- **Correcting for the mediator being estimated** raises it to 0.587 (lineage tree); the uncorrected figure is a lower bound.
- **245 genes keep a significant direct effect** after accounting for the cell mix — candidates for lithium acting inside cells (204 of them up-regulated).
- **The lithium × cell-mix interaction is negligible:** 1 gene of 12,018.

### What the results do and do not show

- **The proportion mediated is not itself evidence of mediation.** With lithium labels shuffled, the same genes still give a median PM of 0.383 (permutation p = 0.49; z = 0.03). It is a valid *description* of how a real total effect splits, not a test.
- **Per-gene mediated effects cannot be declared significant.** The count swings from 0 to 3,647 between trees and the permutation null alone reaches 7,150.
- **The limit is the design, not the data:** 79.5% of the uncertainty comes from the lithium → cell-mix link, so even perfect expression measurement caps the best gene at |z| = 2.92 against the 4.60 needed. Roughly 563 subjects would be required.
- **Unmeasured confounding is the largest threat:** an unmeasured factor associated with both lithium and expression at a risk ratio of 1.30 would explain the typical effect away; age, which we do adjust for, sits at 2.33.

### Added since — circularity sensitivity analysis

- The mediation was re-run on **only the LM22 signature genes** and on **the full filtered set including them**, from one common fit so only the gene set differs.
- Of the 547 signature genes, 6 are absent from the annotation, 191 are essentially unexpressed in blood, and **350 are expressed** — those are the ones analysed.
- **Circularity is real but does not inflate the answer.** The cell mix explains 1.57× more of a signature gene's expression than of an ordinary gene's (partial R² 0.284 vs 0.181), yet signature genes do not show a larger proportion mediated (0.337 vs 0.391; difference not significant), and including them moves the overall figure by -0.0008.

### Added since — alpha and beta diversity

- Transcriptome-wide summaries replace the 12,018 per-gene outcomes: **alpha** diversity (how evenly one person's expression is spread across genes; Shannon) and **beta** diversity (how different two people's profiles are; Aitchison distance).
- **Lithium lowers transcriptome diversity through the cell mix:** indirect effect -0.0341 nats (95% CI -0.0715 to -0.00695; permutation p = 0.020), ≈ 3.4% fewer effective genes. The total effect alone is not significant (p = 0.080).
- **A model-free decomposition agrees:** 95% of lithium's effect on bulk diversity comes from the blood being *reweighted* toward granulocytes, a lineage with a narrower transcriptome; changes inside lineages are not detectable.
- **Lithium makes the cell mix itself less even** (p = 0.003), and does not make people's cell mixes or transcriptomes more variable (dispersion p = 0.91 and 0.50).
- **Beta diversity:** lithium shifts whole-transcriptome structure (p = 0.0003); 67% of that shift is shared with the cell mix. The formal joint-significance test fails on the lineage tree because lithium's effect on the four balances as a whole is weak (p = 0.080), but passes on the proposal tree (p = 0.041). At the pathway level the shift is almost entirely shared with the cell mix and no direct part survives (p = 0.35).
- **Confirmatory family (Holm, 2 tests):** alpha mediation adjusted p = 0.041 (**significant**); beta mediation adjusted p = 0.080 (not significant).
- **Exploratory — bipolar off lithium vs controls:** no change in diversity, cell mix or average transcriptome, but patients may be *more heterogeneous* among themselves on abundance-weighted distances (dispersion p = 0.011, 0.013). A lead, not a finding.

### What changes in the proposal and paper

- Report the proportion mediated with its interval, never a count of mediated genes.
- Correct the proposal's claim that ILR ratios were more consistent across deconvolutions — aggregation from 22 types to 5 lineages did that, not ILR.
- State which cell-type tree is used (both were run; they agree).
- Drop "only" from the myeloid-vs-lymphoid prediction, or keep it with the refinement that granulocytes rise relative to monocytes as well.
- Extract the obesity → neutrophil figure from Furuncuoğlu et al. (2016) for the E-value benchmark; that paper was not in the project library.

---

## 1. The question and the design

Lithium raises neutrophil proportions. Blood is a mixture, so changing the mixture changes bulk gene expression arithmetically, without any cell behaving differently. The differential-expression arm had already found lithium's whole-blood signature to be large and largely compositional. This arm puts a number on it.

| Element | Choice |
|:--|:--|
| Sample | 226 bipolar I subjects: 152 lithium users, 74 non-users |
| Exposure | lithium use, 0/1 |
| Mediator | the cell mix as ILR balances, treated jointly as one vector |
| Outcome | voom log-CPM expression, one model per gene, 12,018 non-signature genes |
| Covariates | age, sex, tobacco, RIN, sequencing plate, sequencing PCs 1–3 |
| Dropped | assessment group — constant within bipolar I, so including it would make the design singular |

The protocol (`PRESPEC.md`) was written before any mediation estimate existed, because the pre-proposal feedback warned that trying exposures, mediators and outcomes until something is significant is a multiple-testing problem. It declares one primary model, every threshold, and a closed list of sensitivity and exploratory analyses.

---

## 2. The mediator: the cell mix as ILR balances

Cell fractions sum to one, so they cannot all enter a regression, and a rise in one forces a fall in another. Isometric log-ratio (ILR) balances on a haematopoietic tree solve both problems: each balance is a log-ratio between two groups of cell types. The written proposal and the implemented pipeline used different trees, so **both were run**.

|   | Lineage tree (pipeline) | Proposal tree |
|:--|:--|:--|
| parts | granulocyte, monocyte, T, NK, B | neutrophil, monocyte, CD4 naive, CD4 memory, NK, B |
| balances | 4 | 5 |
| cell mass used | 100% | 92.4% |
| samples needing zero replacement | 0 | 71 of 474 (15%) |
| median reliability (ICC) | 0.614 | 0.474 |

*Reliability = agreement between our CIBERSORTx deconvolution and Krebs et al.'s deposited CIBERSORT fractions (ICC, absolute agreement).*

| Tree | Balance | Contrast | Reliability (ICC) |
|:--|:--|:--|:--|
| lineage | b1 | gran+mono vs T+NK+B | 0.720 |
| lineage | b2 | gran vs mono | 0.585 |
| lineage | b3 | T vs NK+B | 0.642 |
| lineage | b4 | NK vs B | 0.516 |
| proposal | b1 | neut+mono vs cd4n+cd4m+NK+B | 0.551 |
| proposal | b2 | neut vs mono | 0.654 |
| proposal | b3 | cd4n+cd4m vs NK+B | 0.474 |
| proposal | b4 | cd4n vs cd4m | 0.200 |
| proposal | b5 | NK vs B | 0.448 |

The proposal tree's CD4-naive-vs-memory balance has a reliability of only 0.200 — two deconvolutions of the same blood barely agree on it. This is the main reason the lineage tree is primary.

---

## 3. The estimator, and how it was verified

```
mediator models   M_k = b0_k + b1_k·lithium + b2_k'·C + e_k                      (k = 1..K)
outcome model     Y   = t0 + t1·lithium + Σ t2_k·M_k + Σ t3_k·(lithium·M_k) + t4'·C + e

NDE = t1 + Σ t3_k·b0_k            natural direct effect
NIE = Σ (t2_k + t3_k)·b1_k         natural indirect effect (through the cell mix)
TE  = NDE + NIE                    total effect
PM  = NIE / TE                     proportion mediated (only where TE is significant)
```

These are Pearl's (2001) natural effects in the regression form of VanderWeele & Vansteelandt (2009, 2010), extended to a vector mediator (VanderWeele & Vansteelandt 2014). Because the mediator is a vector, NIE is the *joint* indirect effect through the whole cell mix; per-balance pieces are descriptive, not separate causal effects. Standard errors use the delta method with an exactly block-diagonal covariance, and limma's empirical-Bayes variance moderation for consistency with the DE arm.

| Check | Result |
|:--|:--|
| Monte-Carlo counterfactuals (4 million draws) | matches the closed form to 1.4e-4 |
| independent lm() computation | agrees to 1e-9 |
| four-way decomposition sums to TE | exact |
| centring changes nothing | identical to 1e-8 |
| delta-method SE vs 4,000-replicate bootstrap | within 1.9–6.3% (delta slightly low) |
| 95% interval coverage, 500 datasets | 93.0% (NDE), 94.4% (NIE) |
| product-of-coefficients is a genuine refit, not PIE relabelled | verified |
| total suite | 27 checks, 0 failures |

**The delta method is mildly anticonservative** (~5% low in simulation, ~8% low against the bootstrap on the real data), so bootstrap intervals are quoted throughout.

---

## 4. Stage one: what lithium does to the cell mix

| Tree | Balance | Contrast | Effect | t | p | BH q |
|:--|:--|:--|:--|:--|:--|:--|
| lineage | b1 | gran+mono vs T+NK+B | +0.1558 | +2.00 | 0.0471 | 0.094 |
| lineage | b2 | gran vs mono | +0.1026 | +2.55 | 0.0116 | 0.046 |
| lineage | b3 | T vs NK+B | +0.0354 | +0.56 | 0.5780 | 0.588 |
| lineage | b4 | NK vs B | +0.0629 | +0.54 | 0.5876 | 0.588 |
| proposal | b1 | neut+mono vs cd4n+cd4m+NK+B | +0.1922 | +1.70 | 0.0914 | 0.229 |
| proposal | b2 | neut vs mono | +0.1020 | +2.38 | 0.0183 | 0.092 |
| proposal | b3 | cd4n+cd4m vs NK+B | +0.0245 | +0.15 | 0.8806 | 0.980 |
| proposal | b4 | cd4n vs cd4m | +0.0028 | +0.03 | 0.9799 | 0.980 |
| proposal | b5 | NK vs B | +0.0919 | +0.42 | 0.6732 | 0.980 |

Both the myeloid-vs-lymphoid balance and the granulocyte-vs-monocyte balance move **up**, as the neutrophilia literature predicts (Kast 2008; Krebs et al. 2020). The lymphoid-internal balances do not move. After correction across the four balances only granulocytes-vs-monocytes survives (q = 0.046). The proposal's "myeloid relative to lymphoid **only**" holds in the sense that nothing shifts inside the lymphoid compartment, with the refinement that granulocytes also rise relative to monocytes.

Tested as one multivariate whole (MANOVA, Pillai), lithium's effect on the four balances together gives p = 0.080 (lineage) — weaker than the single strongest balance, because three of the four balances do not move and dilute the test. This matters for the beta-diversity mediation test in section 12.

---

## 5. The primary mediation result

|   | Lineage (primary) | Proposal |
|:--|:--|:--|
| genes tested | 12,018 | 12,018 |
| significant total effect | 1,272 | 1,075 |
| significant direct effect | 282 | 193 |
| significant indirect effect (per gene, delta method) | 0 | 3,647 |
| lithium × cell-mix interaction | 1 | 0 |
| median proportion mediated | 0.396 | 0.449 |
| 95% bootstrap CI | 0.152 – 0.567 | 0.236 – 0.630 |
| share of genes with PM > 0.50 | 30.3% | 38.7% |

### Four-way decomposition (VanderWeele 2014), lineage tree

| Component | Meaning | Share |
|:--|:--|:--|
| CDE | needs neither mediation nor interaction | 53.6% |
| INT_ref | interaction only | 3.8% |
| INT_med | mediation and interaction | 5.7% |
| PIE | pure mediation | 36.9% |

Mediation (PIE + INT_med) accounts for 42.6% of the decomposed effect; the direct part (CDE) is the largest single component.

### Was allowing the interaction worth it?

Only 1 gene shows a significant interaction, so the proposal's original product-of-coefficients estimator was defensible in practice. Refitting it properly — without the interaction columns — it differs from ours by 8.8% of the median effect and disagrees in sign for 2.30% of genes. (An earlier version of this comparison was tautological; see section 14.)

### External validation

The mediation total effect is, by construction, what the DE arm estimated with different code, a different model and a different gene set. They agree: per-gene correlation r = 0.9985; 1,224 of the DE arm's 1,376 lithium genes recovered; Jaccard 0.860. The decomposition splits an independently established effect; it does not invent one.

---

## 6. Can the proportion mediated be trusted?

### Bootstrap — the interval

*B = 1,000 subject resamples; TMM, voom, mediator and outcome models refitted in every replicate.*

| Tree | Genes | Median PM | 95% CI | Resamples above 0 | Resamples above 0.5 |
|:--|:--|:--|:--|:--|:--|
| lineage | 1272 | 0.396 | 0.152 – 0.567 | 99.7% | 14.6% |
| proposal | 1075 | 0.449 | 0.236 – 0.630 | 100.0% | 31.3% |

The pre-declared "substantial" threshold (PM > 0.5) is **not met**: 0.5 lies inside both intervals.

### Permutation — the null

*B = 1,000 shuffles of the lithium labels. The 'median proportion mediated' rows are superseded by the fixed-gene-set null below.*

| Tree | Statistic | Observed | Null median | Null max | Permutation p |
|:--|:--|:--|:--|:--|:--|
| lineage | genes with NIE at BH 5% |     0 |     0 |  7150 | 1.0000 |
| lineage | genes with TE at BH 5% |  1272 |     0 |  1266 | 0.0010 |
| lineage | max \|z\| for NIE | 2.874 | 1.749 |  4.17 | 0.0410 |
| lineage | median proportion mediated | 0.3957 | 0.1729 | 0.7939 | 0.1818 |
| proposal | genes with NIE at BH 5% |  3647 |     0 |  6717 | 0.0060 |
| proposal | genes with TE at BH 5% |  1075 |     0 |   993 | 0.0010 |
| proposal | max \|z\| for NIE | 3.305 | 1.945 |  3.97 | 0.0130 |
| proposal | median proportion mediated | 0.449 | 0.1411 | 0.7058 | 0.1724 |

**The total effect is unambiguous** (p = 0.001). **The mediated-gene count is meaningless**: with lithium doing nothing by construction, a single shuffle produced 7,150 "significant" mediated genes.

### The proportion mediated against its own null

*Gene set held fixed at the genes significant in the real data; 300 shuffles.*

| Tree | Observed | Null median | Null 95% range | Permutation p | z |
|:--|:--|:--|:--|:--|:--|
| lineage | 0.396 | 0.383 | -0.083 to 0.897 | 0.492 | 0.03 |
| proposal | 0.449 | 0.380 | -0.238 to 1.051 | 0.409 | 0.21 |

> **The most important caveat in the arm.** The proportion mediated sits at the same place under pure noise (z = 0.03). PM = NIE / (NDE + NIE) has its numerator inside its denominator, so it lands near 0.4 even when nothing is happening. It is a legitimate description of how a *real* total effect splits — and the total effect is real — but it is **not evidence that mediation exists**. That evidence comes from lithium moving the cell mix and from the indirect-effect tests in sections 12 and 4.

---

## 7. Power — why per-gene answers are out of reach

|   | Value |
|:--|:--|
| share of Var(NIE) from the lithium → cell-mix link | 79.5% |
| median observed \|NIE\| | 0.0393 |
| minimum detectable NIE, 80% power, uncorrected | 0.0645 |
| minimum detectable NIE, 80% power, Bonferroni | 0.1254 |
| largest attainable \|z\| with perfect expression measurement | 2.92 |
| \|z\| needed (Bonferroni, 12,018 genes) | 4.60 |
| subjects needed | ~563 (vs 226) |

You cannot be more certain about a mediated effect than about the lithium → cell-mix link that carries it. That part of the variance shrinks only with more subjects, never with deeper sequencing. This was knowable before any data were analysed, which is exactly why the diversity analysis (section 12) asks the question once, globally.

---

## 8. Sensitivity analyses (declared in advance)

| # | What changed | Median PM | Verdict |
|:--|:--|:--|:--|
| — | primary | 0.396 | — |
| S1 | LM22 genes analysed separately (superseded by section 11) | 0.355 | see section 11 |
| S2 | Krebs' independent deconvolution as mediator | 0.389 | holds; per-gene NIE correlate 0.995 |
| S3 | regression calibration for the estimated mediator | 0.587 | raises it; see section 9 |
| S5 | 20 tobacco imputations, Rubin-pooled | 0.396 | exact no-op: none of the 226 had tobacco imputed |
| S6 | zero-replacement constant 0.10–0.80 | 0.448 – 0.457 | not load-bearing |
| S9 | no covariates / age only | 0.479 / 0.437 | adjustment lowers PM; age does real work |
| S11 | proposal's 6-part tree | 0.449 | agrees |
| S12 | unweighted outcome models | 0.400 | immaterial (NIE correlate 0.9984) |

---

## 9. Correcting for the mediator being estimated

Classical measurement error in a mediator attenuates its coefficient and understates mediation. The proposal commits to regression calibration (Valeri, Lin & VanderWeele 2014) with the uncorrected figure reported as a lower bound.

| Tree | Median reliability | PM uncorrected | PM corrected |
|:--|:--|:--|:--|
| lineage | 0.614 | 0.396 | 0.587 |
| proposal | 0.474 | 0.449 | 0.767 |

Two problems had to be fixed to do this honestly. The full multivariate correction is **inadmissible** — the reliabilities imply more error than the observed covariance can hold — and at the admissibility boundary it is **numerically unstable**:

| Correction strength c | Condition number | Median PM |
|:--|:--|:--|
| 0.000 | 1.00 | 0.396 |
| 0.100 | 1.20 | 0.414 |
| 0.200 | 1.50 | 0.438 |
| 0.300 | 2.03 | 0.464 |
| 0.400 | 3.21 | 0.497 |
| 0.500 | 7.91 | 0.566 |
| 0.561 | 91.86 | 0.877 |

So the stable per-balance correction is reported instead of a point on a cliff edge. The *direction* is certain (every version raises PM); the *size* is not. The proposal tree corrects further only because its mediator is measured worse. The total effect is held exactly fixed throughout, because measurement error in a mediator cannot change it.

---

## 10. Unmeasured confounding — the E-value

| Quantity | Median E-value | At the top gene |
|:--|:--|:--|
| total effect | 1.305 | 2.407 |
| indirect (mediated) effect | 1.233 | 1.531 |
| total effect, CI limit | 1.000 | 2.146 |

An unmeasured confounder needs a risk ratio of about 1.30 with both lithium use and expression to explain the typical effect away. Age — which we adjust for — reaches 2.33 on this scale, so a confounder no stronger than age would suffice. Illness severity, mood state and other psychotropics are plausible at that strength and are not recorded. BMI is not recorded either; the proposal's benchmark (Furuncuoğlu et al. 2016) is not in the project library and was not guessed.

---

## 11. Circularity sensitivity analysis

**Why.** The cell fractions were estimated *from* the expression matrix using the LM22 signature genes, so for those genes the mediator is close to a function of the outcome. The question is how much that inflates the answer.

### Which of the 547 signature genes exist in this data

*Measured against the raw count matrix, not assumed.*

| Status | Freq |
|:--|:--|
| absent from annotation | 6 |
| all-zero | 2 |
| expressed (C2) | 350 |
| present, fails filter | 189 |

The 189 unexpressed genes mostly mark cell types that are rare or absent in blood (regulatory T cells, plasma cells, macrophages, mast cells, dendritic cells). They carried essentially no weight in the deconvolution either.

### Design

One common fit — normalisation, voom weights and the variance prior computed once on all 12,368 genes — so the analyses differ **only** in which genes are summarised. The earlier marker check (S1) had normalised the markers on their own counts, which puts cell composition into the denominator; it is superseded.

### Results

*Total / Direct = genes with a significant total / direct effect (BH within each set).*

| Tree | Set | Genes | Total | Direct | Median PM | Cell-mix partial R² |
|:--|:--|:--|:--|:--|:--|:--|
| lineage | C1 non-markers | 12,018 | 1,273 | 277 | 0.391 | 0.181 |
| lineage | C2 LM22 markers (expressed) | 350 | 48 | 23 | 0.337 | 0.284 |
| lineage | C3 all filtered genes | 12,368 | 1,329 | 300 | 0.391 | 0.183 |
| lineage | C4 every LM22 gene with counts | 539 | 48 | 26 | 0.303 | — |
| proposal | C1 non-markers | 12,018 | 1,083 | 196 | 0.445 | 0.201 |
| proposal | C2 LM22 markers (expressed) | 350 | 47 | 21 | 0.352 | 0.291 |
| proposal | C3 all filtered genes | 12,368 | 1,129 | 213 | 0.443 | 0.204 |
| proposal | C4 every LM22 gene with counts | 539 | 41 | 21 | 0.298 | — |

*C1 non-markers, C2 expressed LM22 markers, C3 everything. 1,000 bootstrap replicates; 300 permutations.*

| Tree | Set | Median PM | Bootstrap 95% CI | Permutation p |
|:--|:--|:--|:--|:--|
| lineage | C1 | 0.391 | 0.148 – 0.562 | 0.488 |
| lineage | C2 | 0.337 | 0.118 – 0.500 | 0.522 |
| lineage | C3 | 0.391 | 0.147 – 0.559 | 0.488 |
| proposal | C1 | 0.445 | 0.227 – 0.620 | 0.392 |
| proposal | C2 | 0.352 | 0.193 – 0.573 | 0.495 |
| proposal | C3 | 0.443 | 0.225 – 0.620 | 0.392 |

*Paired: both sets come from the same bootstrap replicates.*

| Tree | Contrast | Difference | Paired 95% CI |
|:--|:--|:--|:--|
| lineage | C2 markers - C1 non-markers | -0.0542 | -0.1410 to 0.0194 |
| lineage | C3 all - C1 non-markers | -0.0008 | -0.0038 to 0.0008 |
| proposal | C2 markers - C1 non-markers | -0.0927 | -0.1350 to 0.0344 |
| proposal | C3 all - C1 non-markers | -0.0019 | -0.0040 to 0.0014 |

*Normalising on 12,018 genes (original) vs 12,368 (common fit).*

| Tree | NIE correlation | Median PM, original | Median PM, common fit |
|:--|:--|:--|:--|
| lineage | 1.0000 | 0.3957 | 0.3914 |
| proposal | 1.0000 | 0.4490 | 0.4452 |

### What it means

- **The circularity signature is real.** The cell mix explains 1.57× more of a signature gene's expression (partial R² 0.284 vs 0.181; Wilcoxon p < 1e-15). Exactly what circularity predicts.
- **It does not inflate the proportion mediated.** Signature genes show a median PM of 0.337 against 0.391 for the rest — a paired difference of -0.0542 (95% CI -0.1410 to 0.0194). The difference is not significant, and it points the *opposite* way from inflation. Pre-declared criterion 1 (stronger cell-mix link for signature genes) holds; criterion 2 (larger share mediated) does not — exactly the case the protocol named in advance: circularity strengthens the mediator–outcome link without inflating the share mediated.
- **Including them changes nothing.** The full set's median PM differs from the non-marker set's by -0.0008.
- **Why would it lean lower at all? A likely reason, not tested directly:** signature genes carry proportionally more *direct* effect (48% of their significant genes vs 22%). They mark sub-types *within* lineages — eosinophils inside granulocytes, plasma cells inside B cells — and shifts between sub-types are invisible to a lineage-level mediator, so they surface as "direct". The "direct" effect is only direct relative to the resolution of the mediator, which makes the headline PM conservative on this count too.
- **The most specific signature genes are not the most circular** (Spearman ρ = 0.036, p = 0.50): circularity works through the deconvolution as a whole, not gene by gene.
- **The normalisation basis is immaterial** (per-gene NIE correlation 1.0000).
- **Verdict:** excluding the signature genes was the right call, but it was not load-bearing — the ~0.39 headline is not a circularity artefact.

---

## 12. Alpha and beta diversity, mediated through the cell mix

Requested at the weekly check-in. Research and every decision are in `DIVERSITY_RESEARCH.md` and `PRESPEC.md` §8; the decisions were fixed before any of this ran.

### The decisions, briefly

| Question | Decision | Why |
|:--|:--|:--|
| Which meaning of alpha/beta? | All three: per-person / between-person (primary); within- vs between-cell-type (secondary); cell mix's own diversity (secondary) | standard omics usage answers the question with one test; the other two explain the mechanism |
| Location or dispersion? | Both; dispersion mediated only if dispersion differs | PERMANOVA is unreliable under unequal dispersion in unbalanced designs (Anderson & Walsh 2013); ours is 152:74 |
| Exposure | lithium (primary); bipolar vs control (exploratory) | lithium is the arm's pre-specified primary; exposure-swapping is a multiple-testing problem |
| Pathway database | MSigDB Hallmark, 50 sets (primary); KEGG legacy (sensitivity) | built to remove redundancy; 50 interpretable programmes; KEGG legacy is flat, not nested |
| Pathway representation | aggregated abundance by fractional counting, not GSVA | diversity needs non-negative abundances; GSVA scores can be negative and have no entropy precedent |
| Alpha measure | Shannon (primary); Hill q = 0, 2 | the only measure that partitions exactly with unequal weights (Jost 2007) |
| Gene universe | protein-coding, filtered, non-LM22: 10,288 genes | structural non-coding RNAs (RN7SL1 6.9%, RN7SL2 4.8%, RN7SK 3.5%) dominate all-gene counts |
| Beta distance | Aitchison (primary); Jensen–Shannon, Bray–Curtis | same geometry as the ILR mediator |
| Packages | vegan, msigdbr (22.8 MB + 84.6 MB gene-set cache) | GSVA and energy (MODIMA) were not needed |

### Alpha diversity — per person

*Outcomes on the log scale of Hill numbers (nats): an effect d changes the effective number of genes by exp(d). 10,000 bootstrap replicates and 10,000 permutations.*

| Tree | Outcome | Total effect | Total p | Indirect (NIE) | NIE 95% CI | NIE perm. p |
|:--|:--|:--|:--|:--|:--|:--|
| lineage | PRIMARY: Shannon, protein-coding non-LM22, TPM | -0.04 | 0.080 | -0.0341 | -0.0715 to -0.00695 | 0.0204 |
| lineage | richness (all protein-coding non-LM22 genes) | -0.00785 | 0.461 | 0.00301 | -0.00279 to 0.011 | 0.0627 |
| lineage | inverse Simpson | -0.0896 | 0.161 | -0.0701 | -0.158 to -0.0119 | 0.0253 |
| lineage | sens: LM22 markers included | -0.0416 | 0.067 | -0.0363 | -0.075 to -0.00779 | 0.0185 |
| lineage | sens: mitochondrial genes excluded | -0.0407 | 0.078 | -0.0345 | -0.0724 to -0.00706 | 0.0205 |
| lineage | sens: read-count proportions | -0.0331 | 0.006 | -0.0251 | -0.0481 to -0.0056 | 0.0128 |
| lineage | sens: all filtered genes incl. ncRNA (expected to be technical) | -0.092 | 0.110 | -0.015 | -0.0543 to 0.0262 | 0.1122 |
| lineage | PATHWAY primary: Hallmark, fractional counting | -0.00469 | 0.532 | -0.00544 | -0.0132 to -0.000203 | 0.0484 |
| lineage | pathway sens: Hallmark without LM22 markers | -0.00461 | 0.558 | -0.00568 | -0.0137 to -0.000365 | 0.0401 |
| lineage | pathway sens: KEGG legacy | 0.00968 | 0.228 | 0.0104 | 0.00067 to 0.0194 | 0.0222 |
| proposal | PRIMARY: Shannon, protein-coding non-LM22, TPM | -0.0401 | 0.080 | -0.0438 | -0.0811 to -0.0125 | 0.0064 |
| proposal | richness (all protein-coding non-LM22 genes) | -0.00732 | 0.461 | 0.00351 | -0.00426 to 0.013 | 0.1043 |
| proposal | inverse Simpson | -0.0884 | 0.161 | -0.0854 | -0.179 to -0.00392 | 0.0255 |
| proposal | sens: LM22 markers included | -0.0417 | 0.067 | -0.0461 | -0.0844 to -0.0142 | 0.0052 |
| proposal | sens: mitochondrial genes excluded | -0.0408 | 0.078 | -0.0443 | -0.0818 to -0.0127 | 0.0063 |
| proposal | sens: read-count proportions | -0.0331 | 0.006 | -0.0304 | -0.0525 to -0.0122 | 0.0016 |
| proposal | sens: all filtered genes incl. ncRNA (expected to be technical) | -0.0813 | 0.110 | -0.0195 | -0.0656 to 0.0276 | 0.1117 |
| proposal | PATHWAY primary: Hallmark, fractional counting | -0.00463 | 0.532 | -0.00837 | -0.0169 to -0.000759 | 0.0175 |
| proposal | pathway sens: Hallmark without LM22 markers | -0.00453 | 0.558 | -0.00834 | -0.017 to -0.00065 | 0.0174 |
| proposal | pathway sens: KEGG legacy | 0.00853 | 0.228 | 0.0122 | 0.0026 to 0.0209 | 0.0180 |

- **Primary:** the indirect effect is -0.0341 (95% CI -0.0715 to -0.00695; permutation p = 0.0204) — lithium lowers transcriptome diversity by about 3.4% through the cell mix. Interaction p = 0.74.
- **The total effect alone is not significant** (p = 0.080). An indirect effect can be detectable when the total is not: the total also carries a small, noisy direct part, and the indirect test is far more targeted (Zhao, Lynch & Chen 2010).
- **Robust:** with signature genes included, with mitochondrial genes removed, on read-count proportions, on inverse Simpson and on the proposal tree (table). Also after a **post-hoc** adjustment for the structural-ncRNA share — added because primary diversity correlates −0.50 with that library-quality marker: indirect effect -0.0349 (95% CI -0.0722 to -0.00776).
- **Not robust where it should not be:** richness shows nothing (it is depth-driven), and the all-gene version shows nothing because structural non-coding RNAs drown it out.
- **Pathways:** Hallmark diversity falls through the cell mix (barely: CI just excludes zero), KEGG diversity *rises* — pathway-level diversity has no consistent direction, and is reported as descriptive.

### Within- versus between-cell-type decomposition

*H(bulk) = Σ f_k·H(lineage k) + JSD between lineages + residual (Jost 2007; Lin 1991). The five effects add up to the total exactly.*

| Component | Meaning | Lithium effect | 95% CI | Share of total |
|:--|:--|:--|:--|:--|
| H_bulk | bulk diversity (total) | -0.0405 | -0.0874 to 0.00453 | 100% |
| reweighting | MIXTURE: lineage weights change, lineage diversity held at average | -0.0386 | -0.0635 to -0.0153 | 95% |
| within_lineage | CELLS: each lineage's own diversity changes, weights held at average | -0.0117 | -0.0441 to 0.0204 | 29% |
| cross | weights and lineage diversity changing together | -0.00164 | -0.0122 to 0.00822 | 4% |
| between_JSD | MIXTURE: divergence between lineages, weighted by the mix | 0.000461 | -0.00507 to 0.0057 | -1% |
| residual | what the bMIND mixture does not reproduce | 0.0109 | -0.00638 to 0.0299 | -27% |

**95% of lithium's effect on bulk diversity comes from reweighting** — blood shifting toward granulocytes, a lineage with a narrower transcriptome — and that component's interval excludes zero. Changes *inside* lineages are not detectable. This is a model-free confirmation of the mediation result, with the caveat that it rests on bMIND's per-sample lineage estimates (the mixture it builds tracks real bulk diversity at r = 0.94).

### The cell mix's own diversity

| Configuration | Effect on cell-mix Shannon | p | NIE via cell-mix diversity | 95% CI |
|:--|:--|:--|:--|:--|
| lithium \| lineage tree (PRIMARY) | -0.0311 | 0.0027 | -0.0393 | -0.0703 to -0.0138 |
| lithium \| proposal tree | -0.0448 | 0.0030 | -0.0462 | -0.0799 to -0.0161 |
| bipolar vs control \| lineage tree (EXPLORATORY) | -0.00128 | 0.9054 | -0.00196 | -0.0384 to 0.0326 |

Lithium makes the cell mix **less even** — a single number for its evenness carries the whole indirect effect on transcriptome diversity. A shift toward granulocyte dominance also lowers the ceiling on bulk diversity, because the between-lineage term can never exceed the diversity of the mix (JSD ≤ H(f), verified in every sample).

### Beta diversity — between people

*9,999 permutations; 1,000 bootstrap replicates for the shared fraction (Euclidean layers only; Jensen–Shannon and Bray–Curtis have no interval). Like the proportion mediated, the shared fraction is shown only where the shift itself (B1) is significant. Values slightly above 100% occur because an adjusted R² can dip just below zero once the cell mix is conditioned on; read them as 'all of it'.*

| Tree | Layer | Shift p (B1) | Dispersion p (B2) | Direct p (B3i) | Shared with cell mix | Shared 95% CI | Path a p | Path b p |
|:--|:--|:--|:--|:--|:--|:--|:--|:--|
| lineage | gene_aitchison | 0.0003 | 0.497 | 0.0057 | 67% | 0% to 78% | 0.080 | 0.0001 |
| lineage | gene_aitchison_markers | 0.0003 | 0.446 | 0.0043 | 66% | -1% to 78% | 0.080 | 0.0001 |
| lineage | pw_hallmark_aitchison | 0.0088 | 0.318 | 0.3518 | 99% | -10% to 100% | 0.080 | 0.0001 |
| lineage | pw_hallmark_noLM22_aitchison | 0.0062 | 0.262 | 0.3304 | 98% | 7% to 101% | 0.080 | 0.0001 |
| lineage | pw_kegg_aitchison | 0.0011 | 0.163 | 0.0610 | 74% | -7% to 85% | 0.080 | 0.0001 |
| lineage | gene_jsd | 0.0011 | 0.418 | 0.0854 | 80% | — | 0.080 | 0.0001 |
| lineage | gene_bray | 0.0006 | 0.299 | 0.0575 | 81% | — | 0.080 | 0.0001 |
| proposal | gene_aitchison | 0.0005 | 0.497 | 0.0047 | 71% | 5% to 79% | 0.041 | 0.0001 |
| proposal | pw_hallmark_aitchison | 0.0088 | 0.318 | 0.3926 | 101% | -14% to 100% | 0.041 | 0.0001 |

- **Lithium shifts whole-transcriptome structure** (p = 0.0003), and **does not change its heterogeneity** (p = 0.50), so the shift is a genuine location change rather than a dispersion artefact of the unbalanced design.
- **About 67% of that shift is shared with the cell mix**, but the bootstrap interval reaches down to 0%, and lithium's effect survives conditioning on the cell mix (p = 0.0057).
- **Joint significance depends on the tree.** The cell mix strongly structures expression in both (path b p = 0.0001). Lithium's effect on the balances as a whole is p = 0.080 on the lineage tree (fails) and p = 0.041 on the proposal tree (passes). The step that limits per-gene mediation limits this too.
- **At the level of biological programmes, the shift is essentially compositional.** On Hallmark pathways lithium shifts structure (p = 0.0088), but once the cell mix is conditioned on no direct effect remains (p = 0.35); 99% is shared with the cell mix (without signature genes: 98%, 95% CI 7% to 101%).
- **On abundance-weighted distances** (Jensen–Shannon, Bray–Curtis), lithium's shift is again significant (p = 0.0011, 0.0006) and its direct part is borderline once the cell mix is conditioned on (p = 0.085, 0.058).
- **B4 (dispersion mediation) was not run**, exactly as declared, because dispersion did not differ on the primary layer — nor on any other layer for lithium.

*The five leading principal axes of the primary layer, each mediated through med.R (lineage tree). Descriptive.*

| Axis | Variance explained | Total effect p | Indirect p | PM | Correlation with granulocyte fraction |
|:--|:--|:--|:--|:--|:--|
| PC1 | 20.5% | 0.009 | 0.024 | 0.69 | 0.54 |
| PC2 | 11.7% | 0.014 | 0.007 | 0.79 | -0.59 |
| PC3 | 7.6% | 0.299 | 0.773 | — | 0.14 |
| PC4 | 6.2% | 0.030 | 0.033 | 0.59 | 0.38 |
| PC5 | 4.7% | 0.309 | 0.792 | — | 0.04 |

The main axes of transcriptome variation that lithium moves (PC1, PC2, PC4) are 59–79% mediated and track the granulocyte fraction — the beta-diversity shift and the cell-mix shift are largely the same direction in expression space.

### The confirmatory family

*Two pre-declared confirmatory tests, Holm-corrected. Joint significance's p-value is the larger of its two paths (intersection–union test).*

| Test | Raw p | Holm-adjusted p | Verdict |
|:--|:--|:--|:--|
| alpha mediation (NIE, primary outcome) | 0.0204 | 0.0408 | significant |
| beta mediation (joint significance) | 0.0803 | 0.0803 | not significant |

Why can the alpha indirect effect be significant while path a is not? They are different tests. Path a asks whether lithium moves the four balances *anywhere*, spending four degrees of freedom; three balances do not move and dilute it. The indirect effect weights each balance by how much it matters for diversity, which concentrates on the one balance lithium does move (granulocytes vs monocytes). The cell-mix evenness result says the same thing in one number.

### Exploratory: bipolar (off lithium) vs control

- **No change in alpha diversity** (total effect p = 0.623; indirect 0.0111, CI -0.0277 to 0.0493) and **no change in the cell mix's evenness** (p = 0.905).
- **No clear location shift** on the primary distance (p = 0.061; with signature genes p = 0.049) — consistent with the DE arm finding almost nothing gene by gene.
- **But patients are more heterogeneous** on the abundance-weighted distances: dispersion p = 0.011 (Jensen–Shannon), 0.013 (Bray–Curtis), 0.015 (KEGG pathways), with bipolar subjects further from their centroid than controls (Jensen–Shannon 0.1096 vs 0.1006). Not on the primary Aitchison distance (p = 0.369).
- **How to read it:** bipolar disorder off lithium does not move the average blood transcriptome, but patients may differ from one another more than controls do — the "Anna Karenina" pattern (Zaneveld et al. 2017), here confined to the most abundant genes. This is **exploratory**: several distances were examined, it does not hold on the primary one, and the comparison carries the QC-filtered-controls caveat. Worth stating as a lead, not a finding.

---

## 13. Exploratory analyses declared in the protocol

- **T1 — age as exposure (per 10 years).** Age barely moves the cell mix (all balances |t| < 1.4); 90 genes show a total effect of age and their median proportion mediated is only 0.084.
- **T2 — does age's mediated effect differ by lithium status?** No gene differs significantly, and the comparison is underpowered by a factor of 5.78 (74 non-users). A null here is uninformative.

---

## 14. Mistakes found and corrected

Recorded because a methods record that reports only successes is not a methods record.

1. **A validation test was wrong, not the code** — it expected two estimators to agree exactly where they should differ by sampling noise. Replaced with an exact identity plus a z-test.
2. **The total-effect identity appeared to fail** — traced to the deliberately unweighted mediator models; proved to machine precision once the projections match (`02c_identity.R`).
3. **The product-of-coefficients comparison was tautological** (found by independent review): it reused the interaction model's coefficients. Refitted properly; 25.3% → 8.8%.
4. **The regression calibration was inadmissible and the script carried on regardless** (review). Now shrinks to the admissible range and reports a sweep.
5. **The calibration let the total effect drift 18.5%** (review), which is impossible in principle. Now held exactly fixed.
6. **The bootstrap never produced an interval for the headline number** (review). Now computes the pre-declared estimand per replicate.
7. **One review claim was tested and rejected** — a variance 'fix' that is right only conditional on the sample means; the bootstrap showed it would make the standard errors worse.
8. **The calibration was quoted for both trees when it had run on one.** Run on the proposal tree, with that tree's own reliabilities.
9. **The first marker-only check normalised markers on their own counts**, putting composition into the denominator. Superseded by the circularity analysis (section 11).
10. **Path a was briefly coded as RDA** where the protocol declared MANOVA; they test different things. The declared test is used and RDA reported beside it.
11. **Two shell-escaping errors** (a regular expression and a backreference) were caught by assertions before producing results.

---

## 15. What to carry into the paper

1. **Headline:** roughly 0.396 of lithium's whole-blood effect is mediated by cell composition (95% CI 0.152–0.567), rising to 0.587 after correcting for measurement error. Both are lower bounds; they bracket the 0.5 the proposal expected.
2. **Do not present PM as proof of mediation** — present it as a decomposition of a real total effect, with the permutation result beside it.
3. **Do not report a count of mediated genes.**
4. **The cleanest single finding:** lithium lowers transcriptome diversity through the cell mix (NIE p = 0.020, Holm-adjusted 0.041), and 95% of the change is the blood being reweighted toward granulocytes.
5. **State the power ceiling:** ~563 subjects would be needed for per-gene significance.
6. **State the confounding threat plainly:** a confounder as strong as age would be enough, and illness severity and other drugs are unrecorded.
7. **Correct the proposal** on the ILR-consistency claim and on which tree is used.

---

## 16. The genes

### What "significant" means here

The model is fitted once per gene, so every gene has its own total effect (TE), direct effect (NDE) and indirect effect through the cell mix (NIE). Significance is Benjamini–Hochberg at 5% across the gene family. Each gene is placed in one category by asking, in order:

| Category | Rule | Plain meaning | Caveat |
|:--|:--|:--|:--|
| TOTAL | TE significant | lithium changes this gene in whole blood, by any route — the solid list | reproduces the DE arm (r = 0.9985) |
| DIRECT | TOTAL and NDE significant | the change survives accounting for the cell mix: candidates for lithium acting inside cells | 'direct' = not explained by this mediator; sub-type shifts within a lineage would land here |
| VIA_MIX | TOTAL, NDE not significant, NIE bootstrap CI excludes 0 | the change is consistent with running through the cell mix, and no direct part can be detected | a non-significant NDE is absence of evidence; the NIE interval is per gene, uncorrected |
| UNRESOLVED | TOTAL, neither of the above | lithium changes it, and this study cannot say through which route | — |
| OPPOSING | NDE significant but TE not | a direct effect that the indirect effect cancels — the net change is hidden, not absent | rare by nature |
| INTERACTION | lithium × cell-mix interaction significant | the effect of the cell mix on this gene depends on lithium | one gene; the interaction is negligible overall |

*NIE_delta / NIE_boot / NIE_perm: per-gene significance of the indirect effect by three methods. NIE_consensus: significant by all three.*

| Category | Lineage tree | Proposal tree | In both |
|:--|:--|:--|:--|
| TOTAL | 1,272 | 1,075 | 1,074 |
| DIRECT | 245 | 167 | 164 |
| VIA_MIX | 784 | 761 | 624 |
| UNRESOLVED | 243 | 147 | 114 |
| OPPOSING | 37 | 26 | 23 |
| INTERACTION | 1 | 0 | 0 |
| NIE_delta | 0 | 3,647 | 0 |
| NIE_boot | 0 | 3,883 | 0 |
| NIE_perm | 0 | 5,888 | 0 |
| NIE_consensus | 0 | 3,385 | 0 |

> **Why per-gene NIE significance is not a category.** It is 0 genes on the lineage tree and thousands on the proposal tree, from estimates that correlate at 0.98, and the permutation null alone reaches 7,150 genes. Those lists are given in section 16.6 for completeness, flagged, and should not be reported as discoveries.

### 16.1 Direct-effect genes — lineage tree (primary)

*245 genes, ordered by significance of the direct effect. Effects are log2 fold changes in expression (users vs non-users).*

| Gene | TE | NDE | NDE q | NIE | PM | Also DIRECT on proposal tree |
|:--|:--|:--|:--|:--|:--|:--|
| TSPAN2 | 0.786 | 0.629 | 1.2e-21 | 0.157 | 0.20 | yes |
| PIGB | 0.622 | 0.486 | 1.4e-18 | 0.137 | 0.22 | yes |
| RFX2 | 0.67 | 0.524 | 5.3e-14 | 0.146 | 0.22 | yes |
| LINC00877 | 0.576 | 0.459 | 1.6e-13 | 0.118 | 0.20 | yes |
| FAR2 | 0.506 | 0.389 | 7.5e-13 | 0.117 | 0.23 | yes |
| MIR24-2 | 0.661 | 0.522 | 1.9e-12 | 0.139 | 0.21 | yes |
| SLC31A2 | 0.522 | 0.377 | 3.0e-12 | 0.144 | 0.28 | yes |
| MYADM | 0.712 | 0.551 | 4.2e-12 | 0.162 | 0.23 | yes |
| FAM126B | 0.504 | 0.336 | 2.6e-11 | 0.168 | 0.33 | yes |
| FAM63A | 0.461 | 0.296 | 7.0e-11 | 0.166 | 0.36 | yes |
| RBM47 | 0.457 | 0.291 | 7.0e-11 | 0.166 | 0.36 | yes |
| ARMC2 | 0.413 | 0.351 | 4.6e-10 | 0.0622 | 0.15 | yes |
| CCDC19 | 0.763 | 0.576 | 5.9e-10 | 0.187 | 0.24 | yes |
| RAF1 | 0.316 | 0.165 | 2.1e-09 | 0.151 | 0.48 | yes |
| ADAM17 | 0.279 | 0.215 | 4.3e-09 | 0.0637 | 0.23 | yes |
| TAGLN2 | 0.454 | 0.302 | 6.7e-09 | 0.152 | 0.33 | yes |
| ARHGEF40 | 0.522 | 0.36 | 6.7e-09 | 0.162 | 0.31 | yes |
| IL8 | -0.648 | -0.799 | 1.9e-08 | 0.151 | -0.23 | yes |
| KBTBD7 | 0.431 |  0.3 | 2.1e-08 | 0.131 | 0.30 | yes |
| DMXL2 | 0.369 | 0.287 | 3.6e-08 | 0.0822 | 0.22 | yes |
| PRRG4 | 0.677 | 0.467 | 4.9e-08 | 0.211 | 0.31 | yes |
| RAB11FIP1 | 0.443 | 0.274 | 5.9e-08 | 0.169 | 0.38 | yes |
| FLVCR1 | 0.391 | 0.353 | 6.4e-08 | 0.0384 | 0.10 | yes |
| ZNF487 | 0.406 |  0.3 | 1.2e-07 | 0.106 | 0.26 | yes |
| B4GALT5 | 0.47 | 0.275 | 2.0e-07 | 0.195 | 0.41 | yes |
| CYSLTR2 | 0.429 | 0.437 | 2.7e-07 | -0.00824 | -0.02 | yes |
| KLF6 | 0.286 | 0.181 | 3.3e-07 | 0.105 | 0.37 | yes |
| ECT2 | 0.417 | 0.396 | 3.5e-07 | 0.0215 | 0.05 | yes |
| VIM-AS1 | 0.31 | 0.244 | 3.5e-07 | 0.0657 | 0.21 | yes |
| TAGLN2P1 | 0.454 | 0.303 | 3.5e-07 | 0.151 | 0.33 | yes |
| USP4 | 0.251 | 0.141 | 3.6e-07 | 0.11 | 0.44 | yes |
| CKAP4 | 0.453 | 0.29 | 5.9e-07 | 0.163 | 0.36 | yes |
| NHSL2 | 0.443 | 0.291 | 7.2e-07 | 0.152 | 0.34 | yes |
| C12orf55 | 0.655 | 0.491 | 1.0e-06 | 0.164 | 0.25 | yes |
| LRRFIP1 | 0.254 | 0.169 | 1.1e-06 | 0.0846 | 0.33 | yes |
| C3orf62 | 0.324 | 0.176 | 1.7e-06 | 0.148 | 0.46 | yes |
| STX11 | 0.41 | 0.267 | 2.0e-06 | 0.144 | 0.35 | yes |
| TIPARP | 0.252 | 0.196 | 2.2e-06 | 0.0557 | 0.22 | yes |
| CNTNAP3 | 0.906 | 0.628 | 2.4e-06 | 0.278 | 0.31 | yes |
| ANKRD13A | 0.282 | 0.161 | 4.0e-06 | 0.121 | 0.43 | yes |
| MPP7 | 0.269 | 0.222 | 4.8e-06 | 0.0466 | 0.17 | yes |
| IDS | 0.35 | 0.251 | 5.6e-06 | 0.0984 | 0.28 | yes |
| C16orf54 | 0.317 | 0.199 | 6.2e-06 | 0.118 | 0.37 | yes |
| MSL1 | 0.291 | 0.128 | 7.1e-06 | 0.163 | 0.56 | yes |
| PPP1R15A | 0.376 | 0.242 | 9.4e-06 | 0.134 | 0.36 | yes |
| DNAJC3 | 0.375 | 0.243 | 9.4e-06 | 0.132 | 0.35 | yes |
| CPPED1 | 0.454 | 0.277 | 9.4e-06 | 0.177 | 0.39 | yes |
| SLC31A1 | 0.306 | 0.247 | 9.8e-06 | 0.0589 | 0.19 | yes |
| EFHC1 | 0.337 | 0.299 | 9.9e-06 | 0.0379 | 0.11 | yes |
| ARHGAP9 | 0.315 | 0.188 | 1.5e-05 | 0.127 | 0.40 | yes |
| BAZ1A | 0.323 | 0.214 | 1.8e-05 | 0.109 | 0.34 | yes |
| GAS7 | 0.37 | 0.291 | 1.9e-05 | 0.0787 | 0.21 | yes |
| CAMKK2 | 0.29 | 0.182 | 2.0e-05 | 0.108 | 0.37 | yes |
| CAMK1D | 0.355 | 0.256 | 2.0e-05 | 0.0986 | 0.28 | yes |
| ARHGAP19 | 0.293 | 0.185 | 3.0e-05 | 0.108 | 0.37 | yes |
| MED25 | 0.425 | 0.296 | 3.2e-05 | 0.129 | 0.30 | yes |
| FOS | 0.456 | 0.272 | 4.7e-05 | 0.184 | 0.40 | yes |
| CFLAR | 0.252 | 0.129 | 5.7e-05 | 0.123 | 0.49 | yes |
| HIPK3 | 0.257 | 0.174 | 6.5e-05 | 0.0827 | 0.32 | yes |
| VPS37C | 0.439 | 0.346 | 6.5e-05 | 0.093 | 0.21 | yes |
| TATDN3 | 0.267 | 0.214 | 7.6e-05 | 0.0536 | 0.20 | yes |
| FKBP15 | 0.182 | 0.127 | 1.0e-04 | 0.0551 | 0.30 | yes |
| AC116366.5 | 0.398 | 0.208 | 0.0001 | 0.19 | 0.48 | yes |
| CD97 | 0.407 | 0.262 | 0.0001 | 0.145 | 0.36 | yes |
| ZBTB34 | 0.275 | 0.137 | 0.0001 | 0.137 | 0.50 | yes |
| SLCO3A1 | 0.284 | 0.175 | 0.0001 | 0.109 | 0.38 | yes |
| RPGR | 0.31 | 0.189 | 0.0002 | 0.12 | 0.39 | yes |
| GPR157 | 0.434 | 0.331 | 0.0002 | 0.103 | 0.24 | yes |
| TAGAP | 0.261 | 0.167 | 0.0002 | 0.0939 | 0.36 | yes |
| SPAG1 | 0.43 | 0.31 | 0.0003 | 0.12 | 0.28 | yes |
| RP11-73E17.2 | 0.313 | 0.216 | 0.0003 | 0.0972 | 0.31 | yes |
| ANKRD28 | 0.236 | 0.218 | 0.0003 | 0.0179 | 0.08 | yes |
| NFKB2 | 0.24 | 0.179 | 0.0003 | 0.0617 | 0.26 | yes |
| KATNBL1 | 0.329 | 0.221 | 0.0004 | 0.107 | 0.33 | yes |
| KIAA0556 | 0.293 | 0.213 | 0.0004 | 0.0805 | 0.27 | yes |
| MTX1 | 0.346 | 0.196 | 0.0004 | 0.15 | 0.43 | yes |
| TMEM216 | 0.265 | 0.197 | 0.0004 | 0.0679 | 0.26 | yes |
| POLD3 | 0.192 | 0.139 | 0.0004 | 0.0538 | 0.28 | yes |
| MTMR12 | 0.161 | 0.136 | 0.0004 | 0.025 | 0.16 | yes |
| RN7SL600P | 0.42 | 0.272 | 0.0006 | 0.148 | 0.35 | yes |
| RC3H1 | 0.204 | 0.115 | 0.0006 | 0.0886 | 0.43 | yes |
| ASB13 | -0.187 | -0.173 | 0.0006 | -0.0145 | 0.08 | yes |
| TSC22D3 | 0.321 | 0.202 | 0.0006 | 0.119 | 0.37 | yes |
| DNAJC3-AS1 | 0.455 | 0.34 | 0.0007 | 0.115 | 0.25 | yes |
| PIK3R6 | 0.433 | 0.426 | 0.0008 | 0.00669 | 0.02 | yes |
| LINC00963 | 0.317 | 0.188 | 0.0009 | 0.129 | 0.41 | yes |
| TPPP3 | 0.43 | 0.456 | 0.0009 | -0.0262 | -0.06 | yes |
| WDFY2 | 0.21 | 0.163 | 0.0010 | 0.0468 | 0.22 | yes |
| RELL1 | 0.314 | 0.215 | 0.0010 | 0.0994 | 0.32 | yes |
| CASZ1 | 0.301 | 0.287 | 0.0011 | 0.0144 | 0.05 | yes |
| KLHL11 | -0.16 | -0.12 | 0.0012 | -0.0407 | 0.25 | yes |
| DNAJB5 | 0.611 | 0.465 | 0.0012 | 0.145 | 0.24 | yes |
| RP3-354N19.3 | 0.457 | 0.358 | 0.0012 | 0.0989 | 0.22 | yes |
| ENSA | 0.136 | 0.128 | 0.0013 | 0.00849 | 0.06 | yes |
| F5 | 0.411 | 0.262 | 0.0014 | 0.149 | 0.36 | yes |
| SLC9A8 | 0.306 | 0.195 | 0.0015 | 0.111 | 0.36 | yes |
| MBOAT1 | 0.193 | 0.12 | 0.0016 | 0.0729 | 0.38 | yes |
| FUZ | 0.327 | 0.237 | 0.0020 | 0.0894 | 0.27 | yes |
| MCL1 | 0.242 | 0.115 | 0.0022 | 0.127 | 0.53 | yes |
| RP11-179B15.2 | 0.396 | 0.252 | 0.0023 | 0.144 | 0.36 | yes |
| SLC9A7P1 | -0.247 | -0.242 | 0.0024 | -0.00488 | 0.02 | yes |
| CAMK1 | 0.224 | 0.246 | 0.0025 | -0.0216 | -0.10 | yes |
| CELF1 | 0.154 | 0.12 | 0.0025 | 0.0337 | 0.22 | yes |
| RN7SL473P | 0.436 | 0.297 | 0.0027 | 0.138 | 0.32 | yes |
| MEF2C | -0.194 | -0.18 | 0.0028 | -0.0142 | 0.07 | yes |
| CRIP1 | 0.162 | 0.225 | 0.0028 | -0.0639 | -0.40 | yes |
| MKL1 | 0.229 | 0.19 | 0.0029 | 0.0386 | 0.17 | yes |
| SORD | 0.306 | 0.285 | 0.0030 | 0.021 | 0.07 | yes |
| ADAMTSL4-AS1 | 0.476 | 0.325 | 0.0030 | 0.151 | 0.32 | yes |
| KLHL15 | -0.137 | -0.133 | 0.0031 | -0.00381 | 0.03 | yes |
| BCL2 | -0.228 | -0.171 | 0.0032 | -0.0571 | 0.25 | yes |
| CYTH1 | 0.132 | 0.104 | 0.0032 | 0.0282 | 0.21 | yes |
| LDLR | 0.272 |  0.3 | 0.0033 | -0.0288 | -0.11 | yes |
| NSG1 | -0.555 | -0.518 | 0.0035 | -0.0369 | 0.07 | yes |
| CPSF7 | 0.258 | 0.197 | 0.0037 | 0.0614 | 0.24 | yes |
| EFHC2 | 0.374 | 0.271 | 0.0042 | 0.103 | 0.28 | yes |
| CDK2AP1 | 0.277 | 0.237 | 0.0043 | 0.0406 | 0.15 | yes |
| WSB2 | 0.215 | 0.17 | 0.0046 | 0.0444 | 0.21 | yes |
| CASP8 | 0.199 | 0.113 | 0.0046 | 0.0863 | 0.43 | yes |
| PLAGL1 | 0.297 | 0.243 | 0.0049 | 0.054 | 0.18 | yes |
| LEPROTL1 | -0.107 | -0.1 | 0.0054 | -0.00699 | 0.07 | no |
| OSGIN2 | 0.319 | 0.186 | 0.0054 | 0.133 | 0.42 | yes |
| MPZL3 | 0.343 | 0.164 | 0.0055 | 0.18 | 0.52 | no |
| PLEKHO1 | 0.162 | 0.115 | 0.0057 | 0.047 | 0.29 | no |
| SLC29A1 | 0.489 | 0.491 | 0.0057 | -0.0021 | -0.00 | yes |
| NEK11 | 0.297 | 0.27 | 0.0057 | 0.0264 | 0.09 | yes |
| FLCN | 0.22 | 0.166 | 0.0062 | 0.0545 | 0.25 | yes |
| DYRK1B | 0.251 |  0.2 | 0.0062 | 0.0511 | 0.20 | yes |
| GPR114 | 0.242 | 0.284 | 0.0064 | -0.0413 | -0.17 | yes |
| NACC2 | 0.268 | 0.168 | 0.0065 |  0.1 | 0.37 | no |
| RARA | 0.332 | 0.173 | 0.0068 | 0.159 | 0.48 | no |
| BACE2 | 0.303 | 0.309 | 0.0071 | -0.00603 | -0.02 | yes |
| ERMN | -0.195 | -0.196 | 0.0080 | 0.000541 | -0.00 | yes |
| PLK3 | 0.259 | 0.193 | 0.0080 | 0.0663 | 0.26 | yes |
| FGFR1OP | 0.199 | 0.147 | 0.0080 | 0.0517 | 0.26 | yes |
| MAPK1 | 0.196 | 0.092 | 0.0083 | 0.104 | 0.53 | no |
| GIMAP1 | -0.156 | -0.0996 | 0.0088 | -0.0563 | 0.36 | no |
| RTN2 | 0.443 | 0.302 | 0.0091 | 0.141 | 0.32 | no |
| ROPN1L | 0.446 | 0.222 | 0.0096 | 0.224 | 0.50 | no |
| TUBA1A | 0.343 | 0.176 | 0.0097 | 0.167 | 0.49 | no |
| KIF3C | 0.375 | 0.226 | 0.0098 | 0.149 | 0.40 | no |
| TBC1D15 | -0.15 | -0.146 | 0.0100 | -0.00434 | 0.03 | yes |
| FAM134B | -0.214 | -0.173 | 0.0100 | -0.041 | 0.19 | no |
| SEC24D | 0.168 | 0.113 | 0.0101 | 0.0555 | 0.33 | no |
| SSH3 | 0.296 | 0.201 | 0.0101 | 0.0953 | 0.32 | yes |
| UBE2D1 | 0.385 | 0.249 | 0.0105 | 0.136 | 0.35 | yes |
| NUDT5 | 0.184 | 0.104 | 0.0108 | 0.08 | 0.43 | yes |
| PNPLA6 | 0.206 | 0.193 | 0.0111 | 0.0126 | 0.06 | yes |
| GIT2 | 0.166 | 0.0887 | 0.0111 | 0.0775 | 0.47 | no |
| CATSPER1 | -0.268 | -0.272 | 0.0117 | 0.00387 | -0.01 | yes |
| TMEM237 | -0.277 | -0.235 | 0.0119 | -0.0426 | 0.15 | no |
| SEMA7A | 0.294 | 0.264 | 0.0127 | 0.0301 | 0.10 | yes |
| AC093495.4 | 0.389 | 0.214 | 0.0127 | 0.174 | 0.45 | yes |
| C20orf197 | -0.156 | -0.196 | 0.0134 | 0.0403 | -0.26 | no |
| KIAA1841 | 0.149 | 0.168 | 0.0140 | -0.0188 | -0.13 | yes |
| JMJD1C | 0.244 | 0.127 | 0.0140 | 0.117 | 0.48 | no |
| ZNF853 | -0.26 | -0.218 | 0.0144 | -0.042 | 0.16 | no |
| RILPL1 | 0.428 |  0.3 | 0.0151 | 0.128 | 0.30 | yes |
| FAM214B | 0.295 | 0.212 | 0.0152 | 0.083 | 0.28 | no |
| PLVAP | 0.569 | 0.545 | 0.0152 | 0.0239 | 0.04 | yes |
| RP11-215G15.5 | 0.309 | 0.244 | 0.0152 | 0.0646 | 0.21 | yes |
| MARS | 0.147 |  0.1 | 0.0158 | 0.0466 | 0.32 | yes |
| MAML2 | -0.123 | -0.126 | 0.0161 | 0.00225 | -0.02 | no |
| PRKD2 | 0.245 | 0.17 | 0.0161 | 0.0743 | 0.30 | yes |
| ARHGAP26 | 0.273 | 0.133 | 0.0162 | 0.139 | 0.51 | no |
| TBL2 | -0.133 | -0.113 | 0.0175 | -0.02 | 0.15 | yes |
| MPZL1 | 0.41 | 0.229 | 0.0175 | 0.18 | 0.44 | no |
| EPHB1 | 0.433 | 0.213 | 0.0179 | 0.219 | 0.51 | no |
| HIVEP1 | 0.178 | 0.112 | 0.0186 | 0.0667 | 0.37 | yes |
| SUPV3L1 | -0.14 | -0.1 | 0.0192 | -0.04 | 0.29 | no |
| CTA-217C2.2 | 0.372 | 0.297 | 0.0192 | 0.0746 | 0.20 | yes |
| ZCCHC6 | 0.229 | 0.107 | 0.0194 | 0.122 | 0.53 | no |
| RBM8A | -0.0989 | -0.0935 | 0.0196 | -0.00539 | 0.05 | yes |
| SPEN | 0.165 | 0.105 | 0.0201 | 0.0603 | 0.36 | no |
| OBFC1 | 0.121 | 0.131 | 0.0209 | -0.0103 | -0.09 | yes |
| ZNF66 | 0.316 | 0.246 | 0.0209 | 0.0694 | 0.22 | yes |
| IRF1 | 0.232 | 0.127 | 0.0210 | 0.105 | 0.45 | no |
| BEX4 | -0.204 | -0.155 | 0.0211 | -0.0488 | 0.24 | no |
| PPTC7 | 0.176 | 0.0954 | 0.0211 | 0.0809 | 0.46 | no |
| PLCB2 | 0.239 | 0.163 | 0.0216 | 0.0761 | 0.32 | no |
| FRAT1 | 0.303 | 0.138 | 0.0227 | 0.165 | 0.55 | no |
| KLHL24 | -0.0812 | -0.0891 | 0.0229 | 0.00796 | -0.10 | no |
| ARID3A | 0.295 | 0.146 | 0.0229 | 0.149 | 0.50 | no |
| UBE2R2 | 0.217 | 0.106 | 0.0236 | 0.111 | 0.51 | no |
| NR1D2 | -0.159 | -0.101 | 0.0247 | -0.0575 | 0.36 | no |
| EMR4P | 0.358 | 0.379 | 0.0247 | -0.0207 | -0.06 | yes |
| PDLIM1 | 0.288 | 0.265 | 0.0261 | 0.0224 | 0.08 | yes |
| CCDC115 | -0.129 | -0.108 | 0.0274 | -0.0211 | 0.16 | yes |
| H1F0 | -0.311 | -0.35 | 0.0281 | 0.0384 | -0.12 | no |
| RP11-802E16.3 | 0.263 | 0.161 | 0.0281 | 0.102 | 0.39 | yes |
| ARHGEF11 | 0.25 | 0.146 | 0.0293 | 0.104 | 0.42 | no |
| RP11-54A4.2 | 0.479 | 0.32 | 0.0299 | 0.159 | 0.33 | yes |
| PER2 |  0.2 | 0.193 | 0.0305 | 0.00624 | 0.03 | yes |
| CTC-228N24.3 | -0.163 | -0.118 | 0.0305 | -0.0451 | 0.28 | no |
| PDCL | -0.11 | -0.0971 | 0.0315 | -0.0129 | 0.12 | yes |
| HRH4 | 0.441 | 0.397 | 0.0326 | 0.0437 | 0.10 | yes |
| RP11-50E11.3 | -0.14 | -0.147 | 0.0326 | 0.00703 | -0.05 | no |
| CTB-52I2.4 | 0.185 | 0.192 | 0.0331 | -0.00714 | -0.04 | yes |
| ELL | 0.296 | 0.151 | 0.0336 | 0.144 | 0.49 | no |
| DYNC1LI1 | 0.213 | 0.136 | 0.0336 | 0.0771 | 0.36 | no |
| STX3 | 0.284 | 0.101 | 0.0336 | 0.183 | 0.64 | no |
| ERLIN2 | -0.0914 | -0.0801 | 0.0343 | -0.0114 | 0.12 | no |
| LRP10 | 0.313 | 0.171 | 0.0343 | 0.141 | 0.45 | no |
| INPP4B | -0.18 | -0.131 | 0.0344 | -0.0489 | 0.27 | no |
| DHX34 | 0.301 | 0.137 | 0.0347 | 0.164 | 0.55 | no |
| KRI1 | -0.153 | -0.122 | 0.0352 | -0.0306 | 0.20 | no |
| C7orf43 | 0.257 | 0.143 | 0.0352 | 0.114 | 0.44 | no |
| DHRS9 | 0.424 | 0.261 | 0.0359 | 0.164 | 0.39 | no |
| CD58 | 0.213 | 0.135 | 0.0359 | 0.0782 | 0.37 | no |
| C18orf25 | 0.137 | 0.0831 | 0.0368 | 0.0542 | 0.39 | no |
| KRAS | 0.143 | 0.0974 | 0.0371 | 0.0457 | 0.32 | no |
| CARS | 0.135 | 0.074 | 0.0375 | 0.0609 | 0.45 | no |
| NACC1 | 0.239 | 0.189 | 0.0376 | 0.0495 | 0.21 | no |
| MICU1 | 0.152 | 0.0902 | 0.0380 | 0.0621 | 0.41 | no |
| PTPN23 | 0.178 |  0.2 | 0.0386 | -0.0222 | -0.12 | yes |
| MDFIC | -0.163 | -0.103 | 0.0386 | -0.0596 | 0.37 | no |
| ITGA5 | 0.254 | 0.134 | 0.0386 | 0.121 | 0.47 | no |
| RP11-61F12.1 | 0.422 | 0.253 | 0.0386 | 0.169 | 0.40 | no |
| ADORA3 | 0.349 | 0.339 | 0.0403 | 0.00999 | 0.03 | yes |
| TSTD2 | -0.112 | -0.0788 | 0.0419 | -0.033 | 0.30 | no |
| TMCC1 | 0.263 | 0.12 | 0.0419 | 0.143 | 0.54 | no |
| LETM2 | 0.362 | 0.202 | 0.0422 | 0.16 | 0.44 | no |
| PPIP5K1 | -0.235 | -0.189 | 0.0424 | -0.0457 | 0.19 | no |
| AK2 | -0.0836 | -0.0628 | 0.0434 | -0.0208 | 0.25 | no |
| SERTAD2 | 0.184 | 0.144 | 0.0434 | 0.0407 | 0.22 | no |
| ARAP3 | 0.389 | 0.253 | 0.0439 | 0.136 | 0.35 | no |
| CSK | 0.192 | 0.149 | 0.0447 | 0.0428 | 0.22 | no |
| OAZ2 | 0.234 | 0.11 | 0.0450 | 0.124 | 0.53 | no |
| BAIAP3 | 0.477 | 0.351 | 0.0453 | 0.125 | 0.26 | no |
| CLEC16A | 0.174 | 0.17 | 0.0453 | 0.00365 | 0.02 | no |
| AC009961.3 | 0.256 | 0.171 | 0.0453 | 0.0852 | 0.33 | no |
| TMEM69 | -0.11 | -0.101 | 0.0467 | -0.00899 | 0.08 | yes |
| SNORD12 | 0.301 | 0.182 | 0.0467 | 0.119 | 0.40 | no |
| RALGPS1 | 0.182 | 0.154 | 0.0468 | 0.0282 | 0.15 | yes |
| NFKBIE | 0.137 | 0.12 | 0.0468 | 0.0168 | 0.12 | no |
| ZNF770 | -0.1 | -0.0903 | 0.0468 | -0.0101 | 0.10 | yes |
| SEC61A2 | 0.228 | 0.153 | 0.0468 | 0.0746 | 0.33 | no |
| ADSS | -0.102 | -0.0955 | 0.0478 | -0.00599 | 0.06 | no |
| BTN2A1 | 0.205 | 0.129 | 0.0480 | 0.0753 | 0.37 | no |
| MYC | -0.186 | -0.155 | 0.0480 | -0.0306 | 0.16 | no |
| RNPEPL1 | 0.185 | 0.15 | 0.0480 | 0.0354 | 0.19 | no |
| TMEM55A | 0.261 | 0.107 | 0.0480 | 0.154 | 0.59 | no |
| LDLRAD4 | -0.123 | -0.0983 | 0.0480 | -0.025 | 0.20 | no |
| RP5-1041C10.3 | 0.278 | 0.199 | 0.0485 | 0.0786 | 0.28 | no |
| HMGCR | 0.228 | 0.129 | 0.0492 | 0.0996 | 0.44 | no |

### 16.2 Genes with the largest mediated (indirect) effects — lineage tree

*The 50 genes with a significant total effect whose response is most driven by the cell mix.*

| Rank | Gene | NIE | NIE 95% CI | TE | PM | Category |
|:--|:--|:--|:--|:--|:--|:--|
| 1 | TGM3 | 0.331 | 0.09 to 0.567 | 0.535 | 0.62 | VIA_MIX |
| 2 | CNTNAP3 | 0.278 | 0.0848 to 0.515 | 0.906 | 0.31 | DIRECT |
| 3 | KREMEN1 | 0.278 | 0.0757 to 0.467 | 0.562 | 0.49 | VIA_MIX |
| 4 | ALPL | 0.261 | 0.0694 to 0.448 | 0.419 | 0.62 | VIA_MIX |
| 5 | HECW2 | 0.257 | 0.0773 to 0.44 | 0.312 | 0.82 | VIA_MIX |
| 6 | BASP1 | 0.25 | 0.0721 to 0.428 | 0.349 | 0.72 | VIA_MIX |
| 7 | RP11-925D8.2 | 0.243 | 0.0648 to 0.42 | 0.353 | 0.69 | VIA_MIX |
| 8 | RP11-925D8.6 | 0.233 | 0.0611 to 0.397 | 0.346 | 0.67 | VIA_MIX |
| 9 | DYSF | 0.229 | 0.0641 to 0.393 | 0.391 | 0.59 | VIA_MIX |
| 10 | SIGLEC5 | 0.228 | 0.0697 to 0.396 | 0.354 | 0.64 | VIA_MIX |
| 11 | ADM | 0.225 | 0.0627 to 0.389 | 0.38 | 0.59 | VIA_MIX |
| 12 | KIAA0319 | 0.225 | 0.0652 to 0.382 | 0.413 | 0.54 | VIA_MIX |
| 13 | ROPN1L | 0.224 | 0.065 to 0.383 | 0.446 | 0.50 | DIRECT |
| 14 | FBXL13 | 0.222 | 0.0588 to 0.387 | 0.418 | 0.53 | VIA_MIX |
| 15 | EPHB1 | 0.219 | 0.0638 to 0.366 | 0.433 | 0.51 | DIRECT |
| 16 | RP11-925D8.5 | 0.211 | 0.0578 to 0.359 | 0.345 | 0.61 | VIA_MIX |
| 17 | HCAR2 | 0.211 | 0.0568 to 0.353 | 0.357 | 0.59 | VIA_MIX |
| 18 | PRRG4 | 0.211 | 0.0599 to 0.373 | 0.677 | 0.31 | DIRECT |
| 19 | PANX2 | 0.21 | 0.0612 to 0.36 | 0.322 | 0.65 | VIA_MIX |
| 20 | RP11-455F5.3 | 0.21 | 0.0581 to 0.357 | 0.395 | 0.53 | VIA_MIX |
| 21 | LUCAT1 | 0.208 | 0.0632 to 0.345 | 0.312 | 0.66 | VIA_MIX |
| 22 | TGFA | 0.208 | 0.0577 to 0.349 | 0.311 | 0.67 | VIA_MIX |
| 23 | TMEM252 | 0.207 | 0.065 to 0.354 | 0.492 | 0.42 | VIA_MIX |
| 24 | WLS | 0.207 | 0.0397 to 0.381 | 0.469 | 0.44 | VIA_MIX |
| 25 | ZDHHC18 | 0.205 | 0.0589 to 0.34 | 0.242 | 0.85 | VIA_MIX |
| 26 | ASPRV1 | 0.202 | 0.0587 to 0.345 | 0.287 | 0.70 | VIA_MIX |
| 27 | EPHB4 | 0.202 | 0.0576 to 0.371 | 0.413 | 0.49 | VIA_MIX |
| 28 | ABHD5 | 0.202 | 0.0619 to 0.353 | 0.242 | 0.83 | VIA_MIX |
| 29 | MSRB1 | 0.201 | 0.0543 to 0.341 | 0.251 | 0.80 | VIA_MIX |
| 30 | PYGL | 0.199 | 0.0535 to 0.356 | 0.26 | 0.76 | VIA_MIX |
| 31 | B4GALT5 | 0.195 | 0.0548 to 0.336 | 0.47 | 0.41 | DIRECT |
| 32 | RN7SKP16 | 0.193 | 0.0587 to 0.342 | 0.362 | 0.53 | VIA_MIX |
| 33 | RP1-229K20.5 | 0.192 | 0.0567 to 0.324 | 0.352 | 0.55 | VIA_MIX |
| 34 | KIAA1257 | 0.192 | 0.0536 to 0.326 | 0.366 | 0.52 | VIA_MIX |
| 35 | LIMK2 | 0.191 | 0.0477 to 0.321 | 0.315 | 0.61 | VIA_MIX |
| 36 | AC116366.5 | 0.19 | 0.0548 to 0.322 | 0.398 | 0.48 | DIRECT |
| 37 | CTBP2P8 | 0.188 | 0.0437 to 0.334 | 0.44 | 0.43 | VIA_MIX |
| 38 | ANXA3 | 0.188 | 0.0468 to 0.354 | 0.383 | 0.49 | VIA_MIX |
| 39 | RNF24 | 0.188 | 0.05 to 0.316 | 0.309 | 0.61 | VIA_MIX |
| 40 | PDLIM7 | 0.187 | 0.0498 to 0.322 | 0.348 | 0.54 | VIA_MIX |
| 41 | CCDC19 | 0.187 | 0.0516 to 0.329 | 0.763 | 0.24 | DIRECT |
| 42 | DUSP1 | 0.187 | 0.0538 to 0.316 | 0.323 | 0.58 | VIA_MIX |
| 43 | GCNT1P3 | 0.186 | 0.0514 to 0.33 | 0.274 | 0.68 | VIA_MIX |
| 44 | MBOAT7 | 0.186 | 0.0486 to 0.319 | 0.285 | 0.65 | VIA_MIX |
| 45 | IL1RAP | 0.185 | 0.05 to 0.315 | 0.278 | 0.66 | VIA_MIX |
| 46 | UBN1 | 0.184 | 0.0557 to 0.313 | 0.322 | 0.57 | VIA_MIX |
| 47 | FOS | 0.184 | 0.05 to 0.313 | 0.456 | 0.40 | DIRECT |
| 48 | LSMEM1 | 0.184 | 0.0548 to 0.317 | 0.264 | 0.70 | VIA_MIX |
| 49 | CCDC153 | 0.183 | 0.0466 to 0.334 | 0.404 | 0.45 | VIA_MIX |
| 50 | IGF2R | 0.183 | 0.0452 to 0.315 | 0.254 | 0.72 | VIA_MIX |

### 16.3 Every gene with a significant total effect — lineage tree, by category

**DIRECT — up-regulated by lithium** (204)

*ordered by significance of the total effect*

TSPAN2, PIGB, ARMC2, RFX2, MYADM, LINC00877, FAR2, MIR24-2, ADAM17, SLC31A2, FLVCR1, CCDC19, KBTBD7, FAM126B, ZNF487, VIM-AS1, C12orf55, DMXL2, MPP7, EFHC1, PRRG4, TIPARP, ECT2, CNTNAP3, TAGLN2, FAM63A, GPR157, IDS, TAGLN2P1, ARHGEF40, MED25, LRRFIP1, CYSLTR2, CAMK1D, RBM47, TATDN3, KLF6, SPAG1, NHSL2, ARHGAP19, HIPK3, TAGAP, TMEM216, BAZ1A, RAB11FIP1, VPS37C, USP4, DNAJC3, SLC31A1, KATNBL1, C16orf54, GAS7, KIAA0556, CKAP4, CPPED1, B4GALT5, WDFY2, CD97, DNAJB5, RP11-73E17.2, MTMR12, RPGR, PPP1R15A, RC3H1, EFHC2, ANKRD13A, NFKB2, C3orf62, RN7SL600P, FUZ, ARHGAP9, STX11, RAF1, SLCO3A1, ANKRD28, NUDT5, DNAJC3-AS1, POLD3, FKBP15, LINC00963, HIVEP1, CELF1, RELL1, RP3-354N19.3, FOS, CAMKK2, SLC9A8, CPSF7, MTX1, CYTH1, TSC22D3, F5, RP11-179B15.2, CASP8, SPEN, ADAMTSL4-AS1, RN7SL473P, FLCN, WSB2, C18orf25, PLK3, OSGIN2, FGFR1OP, KIF3C, RILPL1, AC116366.5, CFLAR, ZBTB34, JMJD1C, MPZL1, UBE2D1, ENSA, MBOAT1, AC093495.4, NACC2, SEC24D, MKL1, DYRK1B, CDK2AP1, MARS, PIK3R6, ZNF66, LETM2, TUBA1A, KRAS, SSH3, CASZ1, RP11-54A4.2, RTN2, PPTC7, MSL1, CTA-217C2.2, C7orf43, RP11-802E16.3, PLCB2, PLAGL1, SORD, RP11-215G15.5, FAM214B, ZCCHC6, PRKD2, SNORD12, RP11-61F12.1, MPZL3, SEC61A2, AC009961.3, GIT2, HMGCR, MICU1, EPHB1, TMCC1, MCL1, DHRS9, RARA, MAPK1, NEK11, TPPP3, BAIAP3, DHX34, DYNC1LI1, LRP10, BTN2A1, CARS, CD58, ARHGAP26, ITGA5, SERTAD2, PLEKHO1, UBE2R2, ELL, ARHGEF11, FRAT1, RP5-1041C10.3, ARID3A, ARAP3, ROPN1L, SEMA7A, SLC29A1, TMEM55A, NACC1, IRF1, BACE2, HRH4, OAZ2, CSK, RNPEPL1, PDLIM1, PLVAP, RALGPS1, PNPLA6, LDLR, STX3, PER2, ADORA3, CAMK1, NFKBIE, OBFC1, CLEC16A, EMR4P, KIAA1841, CTB-52I2.4, GPR114, PTPN23, CRIP1

**DIRECT — down-regulated by lithium** (41)

*ordered by significance of the total effect*

IL8, KLHL11, GIMAP1, ASB13, BCL2, SUPV3L1, NSG1, KLHL15, BEX4, TMEM237, TBL2, ZNF853, NR1D2, MEF2C, CCDC115, PPIP5K1, TSTD2, KRI1, LDLRAD4, AK2, ERMN, FAM134B, TBC1D15, PDCL, RBM8A, LEPROTL1, ERLIN2, MDFIC, CTC-228N24.3, SLC9A7P1, ADSS, ZNF770, CATSPER1, TMEM69, MYC, MAML2, RP11-50E11.3, KLHL24, INPP4B, H1F0, C20orf197

**VIA_MIX — up-regulated by lithium** (585)

*ordered by significance of the total effect*

CREB1, CEP97, MTCH1, RP11-290F5.1, NUP50, TMEM252, THBS3, KCMF1, RRAGC, FAM107B, SEMA4B, FNIP1, HEBP2, ZNF592, CTC-479C5.12, HSH2D, CREBBP, DSC2, RNF169, RAC1, SIGLEC10, ABCA7, KLF7, CYP4F12, CLN3, CCDC153, RHOB, RNF111, NFKBIA, XKR8, EP300, SNORD12C, IPMK, RNU7-45P, ABHD3, UBN1, DNAJB6, STRN4, ATXN1, PIK3CB, SIRPB1, SGK1, IKBIP, WAS, KDM6B, KATNBL1P6, LPAR2, KREMEN1, TIGD3, MARK2, RHOQ, ZFP36, LRRFIP1P1, RASGRP4, RP11-177C12.1, TCP11L2, PLEKHO2, SH2D3C, WIPF2, TBC1D1, KDM2A, HIPK1, PDLIM7, KIAA1257, RGS14, RP11-274B18.2, C14orf159, IQSEC1, ARHGDIB, ANTXR2, PSENEN, KIAA0319, AE000661.37, SLC8A1, ZBTB7B, PRR14, MED13L, TET3, STARD10, STXBP2, CRK, PAK2, RAB31, RERE, ITGAX, CBL, ATP6V0E1, MAP3K5, IQCE, LPGAT1, FMNL1, EPHB4, GMIP, PTOV1, TM6SF1, RP11-134L10.1, ELF2, PXN, N4BP1, SLC36A1, VSIG10, TMEM43, ZNF467, WDR1, BANP, LASP1, IL6R, SPAG9, PKN2, GALNS, RELB, SRPK2, C10orf12, FAM53B, MAST3, SIPA1L1, CTB-152G17.6, PHOSPHO1, NDE1, APOL2, SP110, CDC123, SPI1, NEK6, DUSP1, ARAP1, UIMC1, RP11-205M20.7, TMEM234, APC, E2F3, UBE4B, RP11-455F5.3, RABIF, ZYX, HCAR2, NUFIP2, TBC1D10B, MOB3A, ZNF33A, GRK6, CTBP2P8, TGFA, DGCR2, PCBP1-AS1, CMTM6, CDK19, SLK, IFNGR2, WLS, SLC35E3, MID1IP1, EHBP1L1, DENND5A, LILRA6, AP003733.1, TYK2, RP11-454L1.2, ASAP1, ATG16L2, UBASH3B, RP2, PILRA, HRH2, CWC25, B3GNTL1, LIMK2, CERS2, LUCAT1, DYSF, ANKRD13D, ADPGK, AP5B1, AKT1, SVIL, RBM5, BCL3, PPP4R1, TGFBR2, LXN, PRDM10, CASC3, UFD1L, MKNK2, GSK3A, TGM3, MGEA5, RP5-981O7.2, RP1-224A6.9, STK40, SLC25A44, SNORD12B, CNN2, PHTF1, TMEM154, SIRPD, SEMA4A, IRF2, PHF21A, CARS2, USP3, CCDC147, ZNFX1, PTK2B, NAP1L4P1, KCNE1, RNF24, SLC6A6, SLC22A1, TFEB, DAPP1, KDM4B, RABEPK, VEGFA, RP11-175O19.4, THOC5, MIR4435-1HG, FBXL13, MAP3K11, EGLN2, CD46, RAB1A, CLIP1, AL512428.1, GLT1D1, UBR2, NR6A1, ARHGAP26-IT1, PIK3R5, PSTPIP1, VMP1, RELT, RP11-448G15.3, RFWD2, VCPKMT, LRRFIP2, BIN3, RP11-326I11.5, NFIC, RP5-1085F17.3, PTBP3, RP1-229K20.5, TIAM2, RP11-68I3.11, MBP, RC3H1-IT1, MARCH7, MPPE1, MAP7D1, FLOT2, NPEPPS, RP11-70C1.1, DEF8, AGO2, DEDD2, RP5-1086K13.1, IP6K1, LAMTOR4, CYB5R4, TRIOBP, PAK1, EZH1, TMX4, MSRA, DCP2, MKNK1, SPOPL, C1orf132, APAF1, ITPK1, PRR14L, NUP214, BRI3, RPS16P5, RSPH3, WHSC1L1, RP11-318K15.2, STXBP5, NRD1, MIR29B1, C1orf63, SAMD4B, ZBTB18, AC116366.6, GBA2, BSDC1, NOTCH1, SBNO2, LRFN1, PLCG2, RP11-288H12.4, MYO9B, LGALSL, FBRS, CCDC17, KLHL21, RP11-212I21.4, RP11-326I11.3, PPP3CA, IL1RAP, MTMR3, USP32, H3F3AP4, STK35, TXN, KMT2D, RER1, CANT1, RP11-288H12.3, SAP130, ALPL, APBB1IP, UNC13D, IFNAR2, HLX, ERGIC1, 7SK, AMD1, ZNF646, CDC42SE1, ADM, RP11-47G11.2, NCOA2, RP11-925D8.2, EGLN1, MBOAT2, SLC19A1, RP11-399O19.9, DDIT3, NSFL1C, RTN3, GNAQ, NDEL1, MIDN, PTTG1IP, HSD17B11, FLVCR1-AS1, FAM45B, CAMK2G, DOK3, RP11-378A13.1, ADAM19, RASSF2, KLF2, DYRK1A, AFF1, USB1, SERPINA1, SLC12A6, RHOG, PTAFR, RP11-446E9.2, MBOAT7, RN7SKP16, VPS8, ST3GAL2, ATG7, SRF, U73166.2, VPS37B, RP11-695J4.2, FOXN3, TFE3, SIGLEC5, GNAQP1, DGCR11, NT5C2, AC083799.1, JMJD6, TNIP1, BRWD3, SNAP23, GPR107, RXRA, VASP, IL10RB, SNX10, RTEL1P1, MEGF9, FAM45A, RNF130, RP11-396K3.1, TMEM50A, ATXN7, SHISA5, ZNF341, IL17RA, WIPF1, TET2, LILRB3, NR2C2, AGO4, AC092620.2, ADAM8, KMT2E, PANX2, TNFRSF1A, SNX11, STX7, TRA2A, RSBN1L, RP11-344B5.2, CTD-2020K17.1, NPTN, NFKBID, BTBD10, CTBP2, NUAK2, DHRS13, ZFP36L1, GLTSCR1L, GPSM2, IMPDH1, TOPORS, MARCKSL1, RAB3D, CD300LB, RP4-647C14.3, SETD1B, KIAA0247, DCAF15, LILRA1, THEMIS2, AC009506.1, BASP1, RP11-473I1.10, GTPBP1, RP11-16E12.1, ATP6V0A1, NFAM1, RP11-925D8.5, SELPLG, MCU, SIK3, SSH2, ST20, AGTPBP1, AREL1, FAM53C, GNB2, ANXA3, SERTAD3, SEC14L1, AC079210.1, ITPRIP, MSRB1, RAPGEF2, C15orf39, TCIRG1, ALAS1, TRIM25, DOCK2, SLC12A9, NIT1, RAC1P2, CTA-217C2.1, RAB7A, RP1-90J20.12, WAC, RN7SL239P, FBXL20, TOP3A, H3F3A, RALB, HGS, HS1BP3, ZDHHC5, FTH1, RP11-925D8.6, SORL1, KIAA0513, TMEM81, BEST1, SIRPB2, LY75, NQO2, CEP63, SCARNA3, LRTOMT, RHOQP2, TRIM56, CSNK1D, MYO1F, AC008984.2, PRKACA, KCNE3, FAM160A2, ATXN7L3, CELF2, GCNT1P3, LSMEM1, ZNF746, RPS6KA5, IGF2R, ANKRD44, DDX17, EXTL3, GOLGB1, WWC3, ARHGEF2, RP11-44F14.4, RP11-373D23.2, DIP2B, GNB1, ARPC1B, MAML1, RP11-349N19.2, RCN3, C6orf62, MAPKAPK2, PYGL, TADA3, RNF122, UBTD1, TCEB1P19, KLF7-IT1, BAZ2B, SHOC2, ABI1, FAM157A, PCYT1A, MAPK13, DUSP18, ZNF281, C16orf72, ELOVL5, TMEM127, TMUB2, HDAC7, MARCKS, CTD-2210P15.2, ZBED6, CHST11, F11R, MAP3K3, ACTR3, ACTN1, DOCK5, GDPD3, GNAI2, NOTCH2, ITSN2, SH3BP2, KIF13A, MPZL2, TIMP2, VRK3, CHRNA10, ASB7, TMLHE, TSC22D4, ARHGAP1, PSEN1, TBL1X, PTPN12, EFHD2, SLC22A4, ZSWIM8, ARPC5, SIRPA, CTBS, ZDHHC18, PHKA2, RHOA, B3GNT5, FAS, ZNF396, TLE4, ASPRV1, HECW2, OPRL1, TANK, FGD3, CYTH4, EPS15L1, GNA13, PATL1, SLC18A2, ABHD5, ICAM3, C16orf93, ANXA11, DKFZP761J1410, STAT3, MGRN1

**VIA_MIX — down-regulated by lithium** (199)

*ordered by significance of the total effect*

PPIH, FER, SLC25A12, C2orf44, XPOT, TRUB1, DPY19L4, AKR1B1, TRIM4, TXLNG, HMG20A, COMMD7, SUCLG2, CLIP4, ACTL6A, TRMT5, DDX28, ARMCX5, MRPL21, CDC23, KLF9, GLOD4, AP3M2, ZNF7, NHP2L1, CLASP2, TARBP1, MAK16, EIF3G, DPH6, WDR3, THUMPD3, CRY1, TPP2, MAT2A, JMY, AC074286.1, HDDC3, NDNL2, SOCS2, ACVR2A, METAP1, QARS, PCCA, DDX18, ARL4C, HDGFRP3, IMPAD1, CAAP1, AC104820.2, STAT4, EXD2, HIST1H1D, UFL1, FASTKD3, NANP, LRPPRC, PEX1, WDPCP, IGIP, PLAC8, RFC1, ARHGEF9, ACADSB, GIMAP6, NCOA7, PTPLB, RNMT, NHLRC3, TWISTNB, ALKBH8, WDR43, TRMT61B, KIAA1324L, CUL1, POLR1B, USPL1, SPIN1, GTPBP4, MKI67IP, RBBP7, HLCS, PDE6B, CENPC, C8orf33, AKAP1, EIF5A2, U2SURP, TMEM97, ADO, EID2, SUMF2, MRPL24, FKTN, AP1S1, KIF3A, SERBP1, ECI2, PSMD5, ALKBH3, TMEM192, GTF3C4, FADS1, FTO, NUCKS1, KDM1A, IMP4, SEPHS1, C1orf109, THOC1, BMS1, TMPO, MUTYH, DLAT, OGFOD3, KDSR, TCERG1, LANCL1, RP11-700P18.1, FAM135A, RPS2, NUDT19, CTDSPL2, AP1AR, DDHD2, PFKM, DIEXF, MTRR, SDAD1, SHPRH, PTDSS1, SLC35B4, POLR3C, CYP2U1, MTIF2, FAM60A, TOMM20, TMEM186, PCNXL4, NT5C3B, PRPS1, C21orf33, PPIA, THUMPD2, HIBADH, UMPS, API5, CIAPIN1, DARS, DDX1, CECR5, DSTN, LRRC58, RP13-977J11.5, DHX32, SLC2A13, NCBP1, DDX10, TBC1D19, AKT3, IMPDH2, SLC16A1, POLR1E, FASTKD2, ATG14, B3GALTL, PEX19, MYBL1, HIST1H1E, NRIP1, HIST2H2AC, POGLUT1, DHX29, LSM11, AHCY, GPATCH2, TIMM50, TOMM70A, MKS1, TTLL1, PRDX1, PTER, EXOC6B, MTX2, EIF3FP3, TMEM161A, C11orf83, UBFD1, SSBP1, SLC20A2, ARMC1, ORMDL1, PCBP4, NARG2, EID1, PRMT3, TRIM44, SHQ1, IRF2BP2

**UNRESOLVED — up-regulated by lithium** (99)

*ordered by significance of the total effect*

RTN4, WDR13, ZMIZ1, TRERF1, TMEM86B, ARRDC4, MICAL1, TTC7B, FCGR1B, FAM193A, SPPL2A, ZFAS1, GNAS, BMP6, RP11-295H24.3, LGALS12, DGKZ, SLC24A3, TTLL3, AF001548.5, RUNDC3A, LYL1, DUS2, MAP2K2, ELF4, CNOT3, INTS3, CDK11A, PTGDR2, VIM, FLJ27365, FAM219A, RP11-815J21.1, VOPP1, REC8, GATAD2B, CLDN15, RIPK2, RIN3, BRD4, TFPI, KSR1, POLR2A, RBCK1, RFX1, RP11-262D11.2, JPX, C5orf56, EPN1, RABGAP1L, EPAS1, PABPC1P3, KIAA0226, KIF1C, LINC01011, NCOA6, SRGAP2B, RP11-73M18.6, RAP1GAP2, MFSD11, NPHP4, EXOC3, LINC00152, IER5, ERP27, CD274, PPP2R3C, PI4KB, INPP4A, YWHAH, DOPEY2, NDUFB3, SLC37A1, MED12, SMOX, ATP2B4, NRBP1, YBX1P6, RAB44, ATP2A3, ARID1A, AKNA, TBC1D25, TREML1, FAM110A, CDC42-IT1, MFAP3L, CLN6, CAMTA2, TLN1, ZGPAT, RPL32P3, ILK, ATF7IP, RP11-320L11.2, CMIP, FAXDC2, DIAPH1, MED15

**UNRESOLVED — down-regulated by lithium** (144)

*ordered by significance of the total effect*

TRAPPC6A, LDHA, QDPR, TNFSF8, GLT8D1, GALNT12, MOAP1, EIF3K, CLDND1, TGIF1, ADAM22, RNF25, TMEM25, EPM2A, H2AFJ, KBTBD3, AC093323.3, ZNF252P, ADNP2, KIAA1279, RP11-553L6.5, C2CD5, CNOT6L, TUBE1, PAIP2B, SURF2, LRRC8C, MPHOSPH10, GLUD1, ALG1, VWA5A, PSMD12, RP11-705C15.2, ACN9, SLC46A3, FNIP2, ZNF512, APEX1, HAVCR2, DST, CTD-2228K2.5, KIAA1107, TMTC2, ZNHIT2, SLC25A36, AP3M1, ATP13A3, FBXO32, RP5-874C20.3, ATP8A1, NIPSNAP1, COPS8P2, HTATSF1, MAPK8, HCCS, MS4A4A, MIR181A1HG, LONP2, ZNF251, DHRS3, POFUT1, BTF3L4, RNF217, ZNF684, ALDH1B1, C3orf18, BBS9, ADK, ZNF225, ZNF736, SRSF8, ACKR3, ALG6, ZSCAN30, PHYH, PLAA, PLD6, CNPY4, C16orf52, CTSF, RNF113A, MB21D2, KCTD3, CEACAM21, RPGRIP1L, APEH, AMOTL1, MINA, PARD6B, RRP1, AKR1C3, SLC16A7, TRPC1, MDH2, PLCD1, PIGK, HS2ST1, SEC62, ARHGAP5, PAQR3, TUT1, DHX57, TRIM68, CLYBL, GIMAP5, NUDT21, C1orf213, SFT2D3, HIST3H2BB, MTMR11, EML5, AUH, ATAD1, CDK5, RPS4XP11, DYNC2H1, ZMAT1, DGKE, BTBD8, ZBTB20, PQBP1, PDCD4, AFMID, GZF1, STAU2, PSIP1, DNAJA2, PDP2, SRSF11, SIRT1, CNOT2, TMEM161B, RP11-705C15.3, CRYL1, MAL, TRBV7-3, CTD-2576F9.1, RP5-935K16.1, KLF12, ZNF3, TNFRSF18, PIK3R1, FAM120C, AC019097.7

### 16.4 Opposing and interaction genes

*OPPOSING (lineage): significant direct effect, cancelled by an indirect effect of opposite sign.*

| Gene | NDE | NIE | TE | TE q |
|:--|:--|:--|:--|:--|
| DGAT2 | -0.178 | 0.231 | 0.0529 | 0.747 |
| RP11-76E17.3 | -0.504 | 0.217 | -0.287 | 0.143 |
| TMEM91 | -0.203 | 0.134 | -0.0692 | 0.538 |
| SLC16A6 | -0.164 | 0.0538 | -0.111 | 0.098 |
| RAD18 | 0.14 | -0.0418 | 0.0986 | 0.082 |
| ST8SIA4 | -0.139 | 0.119 | -0.0202 | 0.836 |
| CD14 | -0.253 | 0.0723 | -0.18 | 0.113 |
| CCDC28A | -0.138 | 0.0533 | -0.0852 | 0.144 |
| RP11-59D5__B.2 | -0.269 | 0.183 | -0.0865 | 0.545 |
| CXCR4 | -0.157 | 0.0812 | -0.0759 | 0.289 |
| CA4 | -0.269 | 0.275 | 0.00514 | 0.984 |
| PELI2 | -0.111 | 0.128 | 0.0169 | 0.871 |
| HSPBAP1 | -0.133 | 0.0653 | -0.0679 | 0.336 |
| ETFA | 0.0994 | -0.0322 | 0.0673 | 0.110 |
| KBTBD11 | -0.165 | 0.0371 | -0.128 | 0.072 |
| SEPT11 | 0.134 | -0.038 | 0.0958 | 0.094 |
| PLB1 | -0.322 | 0.15 | -0.171 | 0.279 |
| SULT1B1 | -0.171 | 0.161 | -0.00986 | 0.948 |
| CACNB4 | -0.356 | 0.0869 | -0.269 | 0.073 |
| FAM217B | -0.0992 | 0.0266 | -0.0726 | 0.092 |
| OSBPL9 | -0.0888 | 0.0396 | -0.0492 | 0.266 |
| PSMD6 | -0.0689 | 0.0133 | -0.0556 | 0.054 |
| GALNT14 | -0.436 | 0.235 | -0.201 | 0.361 |
| SNX18 | -0.124 | 0.0815 | -0.0425 | 0.591 |
| FRRS1 | 0.28 | -0.0527 | 0.227 | 0.054 |
| MIR223 | -0.143 | 0.132 | -0.0117 | 0.928 |
| PHF17 | -0.0956 | 0.0287 | -0.0669 | 0.128 |
| CXCL1 | -0.253 | 0.168 | -0.0845 | 0.590 |
| HDAC5 | -0.127 | 0.0693 | -0.0579 | 0.418 |
| CR1 | -0.167 | 0.178 | 0.0106 | 0.951 |
| CD93 | -0.213 | 0.129 | -0.0837 | 0.542 |
| WEE1 | 0.185 | -0.0421 | 0.143 | 0.091 |
| HIST1H4I | 0.201 | -0.061 | 0.14 | 0.120 |
| TMCC3 | -0.132 | 0.203 | 0.0715 | 0.599 |
| DTX4 | -0.203 | 0.0739 | -0.129 | 0.168 |
| SMAD7 | 0.165 | -0.0368 | 0.128 | 0.138 |
| MCM8 | 0.125 | -0.0597 | 0.065 | 0.308 |

**INTERACTION (lineage)** (1)

FAM3A

### 16.5 Signature (LM22) genes with a significant total effect

*From the circularity analysis (section 11), lineage tree, BH within the 350 expressed signature genes.*

| Gene | LM22 cell type it marks | TE | TE q | NDE | PM | Cell-mix partial R² |
|:--|:--|:--|:--|:--|:--|:--|
| MAN1A1 | Plasma cells | 0.459 | 2.4e-10 | 0.354 | 0.23 | 0.350 |
| RRP12 | Eosinophils | 0.468 | 4.7e-07 | 0.323 | 0.31 | 0.410 |
| IL5RA | Eosinophils | 0.748 | 5.5e-07 | 0.742 | 0.01 | 0.024 |
| MAK | Neutrophils | 0.528 | 1.1e-05 | 0.31 | 0.41 | 0.582 |
| PBXIP1 | T cells CD4 memory resting | 0.268 | 3.6e-05 | 0.189 | 0.30 | 0.230 |
| CDK6 | NK cells activated | -0.203 | 0.0001 | -0.15 | 0.26 | 0.377 |
| FAM212B | Neutrophils | 0.414 | 0.0002 | 0.207 | 0.50 | 0.774 |
| PTGER2 | NK cells activated | -0.183 | 0.0003 | -0.112 | 0.39 | 0.354 |
| DAPK2 | Eosinophils | 0.348 | 0.0008 | 0.154 | 0.56 | 0.711 |
| CLC | Eosinophils | 0.614 | 0.0008 | 0.647 | -0.05 | 0.017 |
| TREM1 | Neutrophils | 0.407 | 0.0008 | 0.191 | 0.53 | 0.730 |
| IDO1 | Macrophages M1 | 0.595 | 0.0009 | 0.613 | -0.03 | 0.063 |
| CEACAM3 | Neutrophils | 0.36 | 0.0016 | 0.162 | 0.55 | 0.699 |
| SMPD3 | Eosinophils | 0.507 | 0.0020 | 0.515 | -0.02 | 0.031 |
| ALOX15 | Dendritic cells resting | 0.799 | 0.0020 | 0.805 | -0.01 | 0.013 |
| C5AR1 | Neutrophils | 0.354 | 0.0020 | 0.155 | 0.56 | 0.729 |
| PPBP | Macrophages M0 | 0.445 | 0.0021 | 0.456 | -0.03 | 0.003 |
| P2RY14 | Eosinophils | 0.376 | 0.0023 | 0.36 | 0.04 | 0.045 |
| EPN2 | Eosinophils | 0.241 | 0.0053 | 0.21 | 0.13 | 0.029 |
| DPP4 | T cells CD4 memory resting | -0.209 | 0.0053 | -0.157 | 0.25 | 0.345 |
| EMR1 | Eosinophils | 0.303 | 0.0063 | 0.304 | -0.00 | 0.200 |
| CSF1 | Macrophages M0 | 0.303 | 0.0096 | 0.241 | 0.20 | 0.100 |
| MXD1 | Neutrophils |  0.3 | 0.0105 | 0.0813 | 0.73 | 0.908 |
| DACH1 | Eosinophils | 0.259 | 0.0109 | 0.217 | 0.16 | 0.067 |
| RENBP | Macrophages M2 | 0.208 | 0.0117 | 0.14 | 0.33 | 0.363 |
| NLRP3 | Monocytes | -0.191 | 0.0147 | -0.223 | -0.17 | 0.398 |
| LRMP | Eosinophils | 0.171 | 0.0149 | 0.0616 | 0.64 | 0.663 |
| QPCT | Neutrophils | 0.317 | 0.0162 | 0.106 | 0.67 | 0.760 |
| PLEKHG3 | Neutrophils | 0.227 | 0.0164 | 0.0825 | 0.64 | 0.543 |
| NFE2 | Neutrophils | 0.207 | 0.0173 | 0.0908 | 0.56 | 0.532 |
| FRMD8 | T cells regulatory (Tregs) | 0.163 | 0.0173 | 0.0912 | 0.44 | 0.393 |
| RCAN3 | T cells CD4 memory resting | -0.165 | 0.0213 | -0.103 | 0.38 | 0.530 |
| CREB5 | Neutrophils | 0.297 | 0.0213 | 0.0987 | 0.67 | 0.621 |
| SEPT8 | Mast cells resting | -0.183 | 0.0213 | -0.116 | 0.37 | 0.149 |
| FCGR2B | B cells memory | 0.292 | 0.0232 | 0.246 | 0.16 | 0.045 |
| IL1RL1 | Mast cells activated | 0.414 | 0.0280 | 0.383 | 0.08 | 0.020 |
| HPSE | Monocytes | 0.211 | 0.0280 | 0.138 | 0.35 | 0.329 |
| HSPA6 | Neutrophils | 0.264 | 0.0300 | 0.063 | 0.76 | 0.749 |
| ZBTB10 | T cells follicular helper | -0.142 | 0.0300 | -0.0826 | 0.42 | 0.449 |
| TRAC | T cells follicular helper | -0.172 | 0.0308 | -0.135 | 0.22 | 0.223 |
| IGSF6 | Macrophages M0 | 0.205 | 0.0324 | 0.0464 | 0.77 | 0.743 |
| HDC | Mast cells activated | 0.705 | 0.0325 | 0.481 | 0.32 | 0.161 |
| SLC15A3 | Macrophages M2 | 0.19 | 0.0332 | 0.0703 | 0.63 | 0.719 |
| GNG7 | Plasma cells | 0.182 | 0.0365 | 0.118 | 0.35 | 0.151 |
| TREML2 | Neutrophils | 0.229 | 0.0414 | 0.0493 | 0.79 | 0.726 |
| IL2RA | T cells CD4 memory activated | -0.169 | 0.0421 | -0.143 | 0.16 | 0.132 |
| ACAP1 | T cells regulatory (Tregs) | 0.0954 | 0.0484 | 0.0801 | 0.16 | 0.051 |
| CD247 | NK cells resting | -0.128 | 0.0484 | -0.0675 | 0.47 | 0.402 |

### 16.6 Per-gene indirect-effect significance — proposal tree (flagged)

> **Not discoveries.** Listed because they were requested and for completeness. The lineage tree gives none of these, and the permutation null reaches thousands. Use only alongside the warning in section 6.

**Proposal tree: NIE significant by all three methods (delta, bootstrap, permutation)** (3385)

*ordered by delta-method significance; each method's full list is in results/genes/genes_proposal_all.csv*

FGR, NFYA, KRIT1, LASP1, SNX11, CFLAR, ARF5, SKAP2, LAMP2, GDE1, JHDM1D, IFRD1, LUC7L, MYLIP, CDKL5, NADK, SPAG9, TRAF3IP3, ZNF200, ABHD5, TTC19, PLAUR, TYROBP, MBTPS2, EXTL3, ELOVL5, MVP, TACC3, ZMYND11, WAS, BID, RALBP1, SLC11A1, MNAT1, SERPINB1, RTFDC1, ALAS1, AKAP11, UBR2, RRAGD, FAS, TNFRSF1B, VEZT, FUT8, GAB2, PIAS1, SLC4A7, TIMP2, ZZZ3, BOD1L1, ZFYVE16, PHF23, TDP1, YTHDC2, KIAA0556, WWC3, RRM2B, CELF2, FAM120A, R3HDM1, VAMP3, TNFRSF9, RFC2, LIMA1, LETMD1, HERPUD1, CYBA, MPHOSPH9, TTC17, ANAPC4, KIF1B, RELT, NOP58, ZNF280C, TMCC3, MTA3, CDK14, YIPF1, NDC1, TBXAS1, TARBP1, SLC2A3, CTDP1, CCAR1, OGFR, ACAA1, DGAT2, VMP1, EIF4B, CASP8, CNN2, SBNO2, IPO5, MAP2K4, CYB5R4, PRKCQ, TBC1D1, DIP2B, NFYC, SPI1, TRMT11, MYO9A, KLF6, TNFRSF1A, RHOA, PDK3, HDAC4, TFE3, ACSL4, ATP11A, KIF2A, NUCKS1, BCL3, RAB27A, NUCB2, DGCR2, JMJD6, NDST1, AP3M2, ST6GALNAC2, EIF2B3, MAP4K4, SEL1L, PRKACA, ACTN1, UBE2D1, MPP5, CHFR, NDE1, MRVI1, PANX2, ALPK1, PDE8A, CLCN4, PTGS2, FRY, PICALM, TSG101, BCS1L, SCARF1, VPS9D1, FNDC3B, FOSL2, DLG1, RAB7A, KIFAP3, KLHL20, ANKRD13A, ICAM3, NT5C2, STXBP2, NFKB2, PPP1R12B, DYNC1I2, APBB1IP, CAPZB, MCCC1, N4BP2, GNB1, MLLT10, NRD1, ITCH, TP53INP2, THOC1, MKNK1, DNM2, KIAA0020, PSEN1, RBL1, PHLPP1, WDR70, RNF13, BCKDHB, KAT6A, ZCCHC6, OXCT1, SMAP2, SSH1, NCOA1, KIF3C, CD59, CD82, HACE1, MAP3K4, PILRA, DDHD2, SNX10, RBM22, TXLNG, LAT2, PPP1R15A, FTL, ATXN7L3, PGS1, UIMC1, LPCAT2, CASS4, EDEM2, FKBP1A, XRN2, MAPKAPK5, TASP1, PXN, GRAMD1A, CMTM1, GMIP, MLF2, RCOR1, YPEL3, MAEA, ICAM1, IRAK3, PDCD7, SPG21, CMTM6, HAUS4, RPGRIP1, TOX4, TINF2, TRPM7, WDR76, SNAP23, NUP50, UPRT, WAC, SRPK1, SH3GLB1, MYO9B, MTAP, MKNK2, MAPK1, PPM1F, CYTH4, TRIOBP, GTPBP1, PACSIN2, TOM1, THOC5, TSPO, MTMR3, NCF4, CSF2RB, ABHD4, RBM23, PRMT5, SOS2, NIN, PYGL, ITPK1, DHRS7, HIF1A, KIAA0247, DICER1, ZC3H14, RPS6KA5, ARHGAP5, CINP, NFKBIA, GMPR2, TRPC4AP, UQCC1, C20orf24, STK4, TPD52L2, DNAJC5, RNF24, RASSF2, CDS2, SIRPB1, DNTTIP1, MAP1LC3A, PIGU, APMAP, ZNF516, SMCHD1, TBL1X, ALG13, AMMECR1, BMX, SLC35A2, RP2, TNFSF13B, KLF5, DNAJC3, VWA8, DHRS12, MAPK3, N4BP1, POLR2C, USB1, SETD6, TANGO6, MON1B, COTL1, USP10, CRISPLD2, CPPED1, PYCARD, STX4, BCKDK, ZNF106, CPQ, RELB, FCGRT, STX10, LILRA1, ASF1B, VRK3, URI1, DENND3, CEACAM4, PLIN3, TYK2, SIGLEC5, RAB3D, JAK3, ELL, ERF, AVL9, NAMPT, PIK3CG, PON2, TTC26, PLEKHA8, IMPDH1, PHF14, YKT6, FKTN, SPIN1, TMEM245, MEGF9, TLE4, KDM4C, DOCK8, UBE2R2, NSMCE4A, PLEKHA1, C10orf54, ANKRD26, MAP3K8, FAM208B, CCNY, PPIF, KRT23, PSMB3, CWC25, CASC3, KPNB1, CBX1, INTS2, PFN1, RNF167, BLMH, CPD, AKAP10, RAB5C, MLX, PRKAR1A, DHRS7B, UNC119, PHF12, TMEM33, UBE2D3, GAB1, KLHL2, TBC1D19, RAPGEF2, KLF3, HTATIP2, C11orf63, EHD1, ST3GAL4, RNF141, HIPK3, CARS, NUP98, TCIRG1, SELPLG, CORO1C, CAMKK2, BIN2, LIN7A, MANSC1, CREBL2, NAA25, LTBR, ARHGDIB, GTF2H3, DDX55, VDR, TBC1D30, CAND1, NUP107, UHRF1BP1L, LPCAT3, LDHB, RAB35, TDP2, DSE, NEDD9, SLC26A8, MAPK14, SOD2, RNF8, ICK, E2F3, QKI, TFEB, DUSP22, EXOC2, PAPD7, HMGCR, COL4A3BP, RNF130, DROSHA, LMNB1, ARRDC3, RAD1, RAD50, ST8SIA4, NUP155, NR3C1, H2AFY, TCERG1, ERGIC1, CPEB4, HRH2, BCL6, NIT2, FAM162A, KPNA1, UBE3A, CEP70, PFKFB4, USP4, ACAP2, GNAI2, GBE1, UMPS, ATP6V1A, ABTB1, MAPKAPK3, ABCC5, ZBTB47, ACTR3, ITGA4, GCA, DOK1, GALNT3, MRPL19, LANCL1, WDR75, EHBP1, CHMP3, IL1R2, NOL10, STRN, DCAF17, SLC25A12, PLEK, ARID3A, NFE2L2, FARSB, RPL22, RNF19B, DOCK7, KIAA2013, RGS2, SRSF11, PHTF1, EXOC8, RRAGC, MTR, SIPA1L2, AKT3, PADI2, CTBS, ECE1, CD46, ATP6V0B, DIEXF, RPS6KA1, ESYT2, ATF6, FASTKD2, KLF7, HMGN3, UBE3D, SGK1, RNF146, FBXL5, TMEM5, UBN1, KLF12, C2orf43, CYP20A1, GTF3C3, UBE2B, HEATR1, PTBP3, RAD23B, FKBP15, CNTRL, FBXW2, NEK6, HSDL2, KDSR, VPS4B, IRF2BPL, AREL1, MLH3, COQ6, RHOQ, OGFRL1, IFIT2, GPAM, FAM45A, CCDC147, GNA13, DUSP1, NUP43, CYSTM1, PLXDC2, TAF12, KBTBD7, FAM53C, APAF1, PPP3CC, TRIM25, SPOP, TSHZ3, PLBD1, PDS5A, CLCC1, GPSM2, ANXA11, BBS9, GLIPR2, RECK, MED13L, WWP1, CD97, ARHGAP9, TTC21B, BAZ2B, LPGAT1, KCNJ2, B9D2, PFKFB2, PREX1, ZNFX1, PPP4R1L, ZNF45, HIST1H1D, NQO2, ATXN1, LRRFIP1, DOK4, IRF1, ARMC7, MBOAT7, GPR108, RTN2, VASP, TGM3, CD93, TMX4, ZNF133, SIRPD, S1PR4, HNRNPR, PCID2, BECN1, ATG14, L3HYPDH, NUP214, HNRNPH2, KDM4B, TUBA4A, PEX1, LRFN1, ZFP36, TST, APOL2, RAC2, DOCK4, LRRC4, MTX2, OSGEPL1, TTBK2, NARG2, ALDH1A2, MBD4, SLC44A2, CDKN2D, SIGLEC9, SEC14L1, TMEM8A, ERMARD, PHF10, MOSPD2, PGLS, TULP4, MLLT1, NDUFA10, ARPC1B, ZSWIM6, LSP1, TNNI2, OSBPL2, CHMP2A, SESN2, MED18, THEMIS2, CLIP1, NOL11, HABP4, CAP1, PPT1, RAB11FIP4, HACL1, SLC6A6, PSME3, NINJ1, IL13RA1, RARA, PRKAB2, CHSY1, TRIM21, RAF1, TIMM10B, PTCD3, MRPL35, PTPRE, TBC1D14, H3F3B, KDM6B, FLOT2, ARHGEF11, SERINC3, PPP1R3D, CASP9, DCTN4, XPO4, ALOX5AP, SCO1, CHIT1, GPALPP1, MORC4, CSNK1G2, RTN3, SBF2, C14orf159, NUMB, CATSPERB, THUMPD3, VAV3, PPHLN1, KIDINS220, LPIN1, ECHDC3, PHC2, AGO4, ELP2, DHX34, TPP2, CARS2, APC, OSTF1, NAA35, TMEM2, AGTPBP1, ADAM19, CCNJL, TRAFD1, SRPK2, MTO1, AKIRIN2, PHF21A, PRRG4, CD63, OS9, REPS1, STX11, DYSF, MDM2, KIAA0513, PCNXL2, EGLN1, COG2, GLUL, STX6, RNASEL, FAM129A, ACBD6, RC3H1, SP110, DOCK10, USP37, CAB39, TMEM127, ANKRD36, CKAP4, PLXNC1, PWP1, ITM2B, RCBTB2, LCP1, SETDB2, RAC1, DBNL, CCM2, MTHFS, TM6SF1, RSAD1, ACTL6A, TRA2B, HLX, RASSF5, IL1RN, ABI1, TOR1B, SLC31A2, TLR4, TEX10, XPA, NCBP1, POLR1E, DENND4C, CNPY3, PPIL1, KIF13A, YIPF3, RNF144B, MTCH1, FAM8A1, PRCP, SDCBP, SORL1, SQRDL, CTDSPL2, EPT1, RAB1A, PREPL, LRPPRC, STAMBPL1, ENTPD1, STAT4, CARF, HECW2, CIR1, SSFA2, WDR12, ABI2, SLC40A1, DIRC2, PARP9, VWA9, PPCDC, AGPAT9, RAP1GDS1, PRDM5, NUP54, ANXA3, PPP3CA, RNF185, ETV6, GABARAPL1, YARS2, FGD4, FAM60A, C1RL, SLC15A4, GIT2, NUPL1, CSAD, PELI2, JDP2, SLC12A6, TCF12, TLE3, ANP32A, HMG20A, NCOA2, ABHD2, MCTP2, IQGAP1, ITGAX, CPNE2, CMTM2, MBTPS1, ZCCHC14, DEF8, PCTP, NPEPPS, SSH2, ARSG, SS18, TAF4B, IMPA2, ARRB2, SLC16A3, CSNK1D, NARF, SECTM1, VAV1, TMEM91, IFNAR1, AKT1, MYO1F, NLRP12, ZNF787, EFHD2, PGD, SH3BGRL3, BCL10, PRKACB, CD53, RFWD2, FCGR2A, TOR1AIP1, TUFT1, SNX27, MCL1, FAM63A, CERS2, EIF2D, TPM3, SLC39A1, RIT1, LYST, NVL, DEGS1, ARF1, LBR, RHOB, MRPS5, EXOC6B, SFXN5, RALB, SPOPL, NBEAL1, RAB5A, CTDSP1, CSRNP1, IQSEC1, SHQ1, NFKBIZ, SRPRB, OSBPL11, CENPC, PLAC8, ENOPH1, TBCK, SPATA5, ROPN1L, BDP1, ARHGAP26, RNF145, FBXO38, TNIP1, G3BP1, FARS2, CDKAL1, DOK3, PRIM2, RPL7L1, SHPRH, SLC12A9, TMEM140, GPR174, RBMX, MCPH1, ZNF185, RPL10, ATP6V1B2, BIN3, DOCK5, TACC1, UHRF2, HIATL1, NR6A1, NOTCH1, NACC2, NRBF2, CAMK2G, RPP30, DNAJB12, ITPRIP, ADM, PTPRJ, PAK1, LAMTOR1, MPZL2, COMMD7, CDC42EP2, PPP4C, TMCO3, LATS2, SEC24D, RILPL2, SLC7A11, IPMK, FER, CDC123, THYN1, FAM160B1, ADAM8, FLI1, WWC2, ACSL1, CACUL1, GLT1D1, ARL11, TADA1, JMY, SUV39H2, NMT2, GPD1L, PELO, GPR180, HNRNPDL, ANAPC1, RASSF3, FAM49B, UBP1, RMND5A, TUBGCP5, TOMM70A, MRPL39, PPP4R1, MPPE1, EPHB1, TMEM55A, CDK19, MARCKS, TTC39B, NUP205, PIK3AP1, RBM45, C10orf12, FAM126B, SLA, MCU, GNAQ, DPY19L4, BACH1, ZFAND3, NPTN, RAB11FIP1, HIAT1, GHRL, EXOG, ST3GAL2, CACNA1D, TSC22D3, DYRK1A, KCNJ15, ETS2, SLC37A3, MRPL17, TPRG1L, EIF5B, B4GALT5, NCF1, ZFAND2B, PPP1R15B, TAGLN2, ELK4, F11R, USF1, ALG8, IFNAR2, IFNGR2, BTG2, HK2, THEM4, RGL4, ATP6V0D1, ZYX, SSU72, SHKBP1, DEDD2, MPZL3, AMICA1, ZBTB7B, NLRX1, ADAR, IL6R, SLC25A44, IER2, ZNF333, ORAI2, SH2B2, FBXL13, PCYT1A, ACOX1, ITGA5, CXCL16, STX5, NOL9, KLHL21, LAPTM5, ALPL, FUBP1, GBP2, AGL, ZNF281, ARPC5, BROX, MAPKAPK2, NUP35, C2orf44, H3F3A, SMARCAD1, TNFAIP8L2, RNF149, CDC42EP3, S100A11, ARHGAP25, TGFA, SLC22A15, EIF4E3, PROK2, ARPC2, TGFBR2, CLASP2, NUAK2, IFI16, RYBP, C3orf17, WDFY3, PPM1K, GMPS, SLMAP, RBM47, MTMR14, CXCL1, SNRK, PLB1, WDR43, ZDHHC3, TKT, PRKCD, METAP1, SHISA5, PPM1M, NAA15, PRRC1, F2RL1, UTP15, CREBRF, SFXN1, TAGAP, BRI3, LMTK2, SUN1, SLC25A32, TP53INP1, UBAP1, SYK, ABCA1, NFIL3, MID1IP1, NCF1C, DDX26B, WRN, SUGT1, RPUSD4, METTL17, ARHGEF40, METTL3, TRUB1, FRAT1, UBTD1, TC2N, TCP11L2, ALKBH3, ZNF202, WBP1L, TAF10, RNF169, PRKCB, CLEC4E, NDEL1, AP1G1, ANPEP, STAT6, PATL1, STX3, MTMR10, DIS3L, CCNDBP1, CNPY4, PHB, CRK, POP5, ZNF646, CA4, RAB8A, MIDN, ZNF641, DHRS13, TUBA1A, SERTAD3, PLIN4, SLC43A2, BEST1, FTH1, ATG16L2, TRANK1, FADD, CCDC88B, PAFAH1B2, RBPJ, ZCCHC4, TTC39C, PXK, IRF2, MLKL, RAB31, TMUB2, STAT3, TET2, HNRNPH1, XPO6, RGS14, LMAN2, RAB24, CA5B, PGM2, PTAFR, AGPAT2, ITGAM, TOR1AIP2, MAP3K2, MYO7B, TMEM154, NIPA1, USP38, SLC16A5, USP47, GABARAP, PRDM10, FOS, PFKFB3, USP32, PPM1D, GPR27, KBTBD2, KIAA0232, TMEM43, TSEN34, NAIF1, TRIM8, JUNB, LRG1, CANT1, CHST11, POLR1C, LRRC8C, PIK3CD, RGS19, RASGRP4, PTEN, ZNF217, JMJD1C, RAB33B, MOB3A, MALT1, ISG20, MBOAT1, CEBPB, SUCLG2, GNB2, PPP1CA, ZMAT3, CFL1, TMCC1, NAA16, FAM192A, DCP2, MYD88, MIR4435-1HG, TADA2B, ADRBK1, KDM2A, ESRRA, MTX1, PPP1R3B, MAP3K11, TRIB1, EHBP1L1, NABP1, SULT1B1, NMNAT1, SLC19A1, STAT5B, SLCO4C1, CEP19, GNG5, TLR1, TLR6, FAM174A, MGA, PITPNA, RALGAPA1, MSL2, KIAA1551, NR1D2, CTBP2, ATR, SH3BP5L, TP53I11, ZNF519, ZNF25, CCDC14, LRRC25, KCNE3, ERCC4, B3GNTL1, AURKAIP1, ORAI3, IP6K1, RTTN, EID2, KCMF1, SLCO3A1, KIAA2018, B3GNT5, BASP1, METRNL, PNMA1, RHOG, ZBTB34, TALDO1, NAP1L4P1, B3GNT8, TOP3A, IL17RA, AGTRAP, ANXA2R, GRB2, RP1-232L22__B.1, GALNT11, THBD, OPLAH, DPY19L3, C17orf62, ZNF552, SNX18, RP11-349N19.2, ZNF664, NSUN7, CLK3, ZBTB18, GCC1, PCBP1-AS1, MYADM, OAZ2, PNPLA1, HCLS1, ZNF609, CCDC66, KCNE1, C9orf139, FUT7, HIST1H2AC, CUEDC1, KCTD2, ST20, LSMEM1, DHTKD1, ZNF746, FRAT2, DDX60L, AATK, ZNF467, RNF135, TMEM252, RELL1, RNF41, IST1, NCF1B, LIMK2, TTC3, HCAR2, C16orf72, CEP63, SLC8A1, ABAT, RABIF, PTTG1IP, DAZAP2, CECR6, CCDC125, ASB7, MX2, ZNF438, FAM101B, TMEM50A, KREMEN1, DENND5A, CMSS1, WDR27, PDE4B, SNN, APOBR, UBE2G2, ZFP90, SEMA4B, FAF1, IFITM2, TNFAIP2, PRMT3, KIAA0825, UBALD2, TARSL2, SP1, EVI2B, C16orf54, SETD4, IRS2, BICD2, ZNF566, MARC1, RXRA, FCAR, ARHGAP30, CYP4F3, LILRA5, SIRT7, TLR5, MAGEH1, TET3, SOWAHD, ZFP69, TMPPE, ZC3H6, LAMTOR4, C19orf35, PLSCR1, C3orf62, H3F3C, RALGAPA2, MTF1, NHLRC3, MSL1, LRRK2, KCTD21, LITAF, TMEM120A, HN1, KAZN, TMEM194B, IL1RAP, ACADSB, STK40, SEMA4A, SIRPB2, PPIA, IARS, GIMAP5, CD55, NTNG2, EPHB4, MAN2A2, MME, SULF2, TECPR2, HSH2D, MAML3, PPTC7, KIAA1211L, PDLIM7, CASP4, IGF2R, SLC22A4, SERPINA1, FITM2, LRP10, PELI1, MAP3K5, MYO5A, SIPA1L1, PIGN, TOPORS, C5orf42, CDC42SE1, RPE, KLHDC1, SLC9A8, FAM49A, ADH5, PLCG2, S100A6, MPZL1, VPS13A, MBP, FCGR1B, MAK16, SIRPA, MAN1A2, HSD17B11, CACNA1E, TMEM116, SDAD1, ZKSCAN8, PIM3, ASPH, WWP2, TXNRD1, LRBA, BAZ1A, TTC37, FAM179B, ECI2, F5, MSRB1, FAM169A, TMEM184B, LRIG2, GK, RUNDC1, MAP3K3, TBKBP1, PJA2, PDCD4-AS1, CR1, FCGR3A, RP11-290F20.1, LINC00963, NHSL2, ZDHHC18, SPIN3, LILRB3, SLC35B4, PSENEN, IPO7, NAP1L4, HMGN1, WDR52, MIR223, DKFZP667F0711, RP11-288H12.3, GIMAP1, MXD3, ARHGAP19, SIPA1, SRA1, LTB4R, AC078899.1, GNAQP1, PLEKHM1P, C19orf38, PLIN5, STARD10, RP5-1053E7.3, TTC3P1, UBXN2B, DDX3X, FNIP1, RP1-229K20.5, FTH1P8, CEBPD, FAM45B, RN7SKP16, RP4-717I23.3, SNHG14, AC009506.1, RP11-290F20.3, UBE2D3P1, NPM1P39, AC079630.4, KB-1027C11.4, LINC00694, TMEM185B, USF1P1, FTLP3, RALGAPA1P, NAMPTL, RP11-456P18.2, GK3P, AC084290.2, NFE4, GCNT1P3, RP4-673D20.3, SDCBPP3, SNHG15, FAM157B, RP11-295G20.2, RP11-454L1.2, AC104820.2, AC087499.5, NFAM1, H3F3AP4, KIAA0040, TTC28-AS1, RP4-814D15.1, RP11-59D5__B.2, SUCLG2P2, FAM157A, PTENP1, RP11-443B7.1, AC116366.5, RP11-550F7.1, LINC00877, ARPC4, ARPC1A, PLEKHO2, PISD, RP11-242C19.2, ARFGAP3, RP4-576H24.2, RPS2P44, IL10RB, PBX2P1, LILRA6, AC096772.6, ASPRV1, LINC00861, LINC00920, OIP5-AS1, SDCBPP2, TWF2, AC008984.2, LUCAT1, TOX4P1, CCDC153, RP11-701P16.1, RP11-321E2.4, PVT1, RAC1P2, RP11-455F5.3, ZNF345, RTN3P1, RTEL1P1, DDX3YP3, RP11-875O11.3, TAGLN2P1, ALG1L13P, LYN, AP5B1, RP11-676M6.1, CHMP1B, RP11-832A4.5, RP11-750H9.5, HCAR3, SNHG1, RP11-561P12.5, RP11-695J4.2, CAPNS2, RP11-1220K2.2, RP11-256L6.3, CUX1, RP11-443B7.3, RP11-299A16.1, RP11-73E17.2, CEP95, RP4-647C14.3, RANBP20P, RP11-925D8.5, RP11-925D8.3, RP11-925D8.2, CTD-2330J20.1, RP11-925D8.6, AC083843.1, RP11-4F5.2, TPRKBP2, FAM157C, 7SK, CCPG1, RP11-38G5.4, RP5-1057J7.6, RP11-473M20.9, RP11-473I1.9, RP11-473I1.10, AL512428.1, CTB-167B5.2, PGDP1, CTD-2267D19.2, CTC-510F12.4, RP11-47L3.1, KCNJ2-AS1, CTC-503J8.6, CTB-61M7.2, LINC00528, CTD-3023J11.2, RP11-65L3.2, RP11-326I11.5, RP11-373D23.2, RP11-288H12.4, RP4-769N13.6, RP4-673D20.6, RP11-446E9.2, RP11-134L10.1, CTA-384D8.34, NUP160, SMCR8, RP11-792D21.2, SH2D3C, RQCD1, DKFZP761J1410, UTP20, POR, VPS13C, TANK, DUSP19, NOC3L, NLRP6, TRMT2B, CEP290, RP11-290F5.1, CTD-2530H12.2, RN7SL600P, TMEM38B, DUS4L, ELF2, MAST3, PTPDC1, MARK2, CHIC2, CMAHP, AC093673.5, PIK3R5, OTUD4, PPIH, ZDHHC12, AMPD2, SAT1, ZNF585B, ASNS, ANTXR2, RBM5, GCFC2, AGK, JARID2, BRCA1, TTC27, TYMP, IFNGR1, FAM160A2, SP100, PPAP2A, HYAL2, DHRS9, CLNS1A, PTPLAD1, FMO4, DDX1, STX7, IPO11, DDX18, NSFL1C, STRN4, MSH2, BTAF1, ACTR5, DHX35, LPIN2, NDFIP2, MGRN1, KAT8, INTS10, RABAC1, COPE, HBP1, TADA2A, RPL19, GALNT7, FNBP4, CBL, SLC38A1, HCFC2, VEGFA, CLIP4, DNAJC16, KIAA1324, ARHGEF2, DPH5, IFIT3, CCDC91, ATP7B, NAGK, ZNF331, AOC2, ILKAP, SDF2, OSER1, DPH6, SYNCRIP, HNRNPA1, APTX, NARS2, ACTR1A, HERC3, PARVG, VPS37B, SLAIN1, RAB20, UBE2Q2, PSTPIP1, FTO, PDPK1, ZNF287, ZFP14, POGK, ZNF687, WDPCP, TRPC1, TRIM4, ADAM17, CWF19L2, SPEF2, LONRF1, TSEN2, LSM11, ZDHHC5, MAPK13, ZFYVE9, FCER1G, PAXBP1, NAE1, G6PD, NAPEPLD, SPICE1, TIPARP, CCR1, ZNF589, BRWD3, FBXO33, NUDT5, API5, RIMKLB, SAAL1, GFM1, FAM153A, ACYP2, PAQR8, TRMT61B, BANP, RAB37, PDE12, VCPIP1, PCCA, MCTP1, GBA, PAK2, CAPS2, TIGIT, RAD51B, EXT1, IGIP, NDNL2, CXorf38, IFITM1, KPNA4, MRPS21, PTMA, BLOC1S5, ZNF785, SVIL, AC083799.1, RUFY2, STK38L, SNORD12, SUPT4H1, LYRM4, SH3D21, CKLF, RP5-1086K13.1, AHCTF1P1, RP11-297L17.6, ADSL, ZNF487, RP11-582J16.5, ZNF486, RP11-25K21.6, POLA1, NGDN, CARHSP1, CRELD1, SFMBT1, RASGRP1, SYNJ2BP, RPL19P12, CTD-2037K23.2, CTC-523E23.11, SEPSECS, RP3-402G11.26, NUBPL, CYHR1, RRP15, SUPT7L, SEC31A, ABCG1, RAB4B, SPEN, PGM1, ZNF213, TMEM14A, DECR1, MICU1, SH2B3, ASCC3, CREB1, CASD1, UBE4B, KIF3A, CHD1L, DSC2, DGKE, IKBIP, S1PR1, PGAM1, RAB1B, NKRF, B3GALTL, ZNF548, RP11-380G5.3, RAB43P1, RP11-61F12.1, CTB-31N19.3, LPP-AS2, GSK3B, RP11-809N8.4, HLTF, LRCH4, TOLLIP, HPCAL1, PCNXL4, FBXO25, KAT6B, GLMN, YTHDF3, LIN52, RPL18AP3, LINC00211, KLHL42, ZMAT5, IL1R1, AMIGO2, RP11-563J2.2, TNFSF10, ZMYM4, MAP4K5, MDH1, MAPK9, CUL1, CD84, PIGB, DAPP1, WIPI1, TUBE1, CEACAM1, FYB, RRN3, NFX1, PTPN4, RBM41, RFFL, PCNX, USP14, PRPS2, ELP4, HSPA8, C12orf57, ASF1A, HARS2, MRPS27, THG1L, EIF1B, CYTIP, SRSF7, GPX7, WLS, CD58, ANKRD13C, SOCS2, ZRANB3, USP45, SLC36A1, PANK2, RBCK1, EMC1, STARD3, USPL1, TMEM243, DYNC1LI2, LIMD2, TOR1A, TSTD2, KIAA0319, EXOC6, SLC38A6, IGF1R, KATNB1, MAP1LC3B, MRPL24, SDE2, HAUS6, BUB3, SUPV3L1, SAMD8, PSMD4, C1orf123, ARFIP1, KIF27, ZNF319, CEP120, TADA3, FAM98B, PLRG1, TCEAL1, SPRYD4, DDX10, RP11-179B15.2, NPM1, ADO, C8orf33, DCAF4L1, LINC00085, HDDC3, XPOT, DGAT1, NAP1L1, ZNF573, ZNF682, ZNF248, IPO9, CHIC1, ZBTB48, TRBC2, CCDC7, RP3-425P12.4, RPL10AP6, RPS14P3, DIRC3, CAP1P2, RP11-47G11.2, NKAPP1, RP11-229P13.20, TTN-AS1, AC092620.2, RP11-44F14.4, RP11-296I10.6, RP11-212I21.4, AC079210.1, RP11-1109F11.3, POLR2J, FAM214B, CREBBP, CROT, TNFRSF12A, COX10, MATR3, SLC45A4, PLEKHO1, HSF2, RNASET2, TFB1M, FAM136A, RFC1, AGA, TRIT1, CTPS2, THOC3, PPP1R12A, POLR3E, GNAI3, TESK2, MGAT4A, WDR1, IKBKG, KDM5A, LMAN1, PALB2, TDRD3, SH3BP2, DOCK9, SNX5, RBBP9, OSBPL8, WDR7, UNC13D, GPATCH2, HIVEP1, JAG1, RNF125, RGCC, CORO1A, MTFMT, ARMC1, FBL, RPS19, HIBADH, SETX, CCSER2, SHOC2, GALK1, C1QBP, WSB1, UFSP2, ATG2A, MRPS30, ARL6, SLC25A36, NCBP2, KIAA1257, ALMS1, WDR77, SLC16A7, MSANTD2, EXOSC8, KDM3B, TMPO, PDLIM2, SEC22A, ZNF639, MTERFD2, FAM126A, AMD1, EXOSC9, OPRL1, KTN1, RPL36, CAMSAP1, HELZ2, ZNF337, POMT1, ZNF236, TMEM204, MCCC2, COG6, PRAM1, LARS, DPF2, IAH1, DOCK2, MYBPC3, CLOCK, MDFIC, CREBZF, ALKBH8, ATRAID, TMOD3, CCNG2, REM2, WDR61, GLYR1, RPS2, SMYD2, CAMKMT, DIS3L2, EAF1, TBC1D32, AC004453.8, C7orf43, GSN, TMEM218, TRAPPC8, ZDHHC7, SLC26A2, PHF6, FBXO32, GALNT14, DEDD, ADPGK, ZNF235, C21orf33, SAP30BP, GPR155, PIGX, SLC9B2, ARSK, KIAA1430, RICTOR, STXBP5, OSGIN2, YWHAZ, METTL2B, C9orf89, EML5, FAAH2, BEND7, PDZD8, BTRC, ATF7IP2, TSC22D4, C15orf39, RPL13, SIGLEC7, SLC33A1, ZNF562, AFF1, RARG, ZNF483, CCDC41, ANKRD36C, CTDSP2, GOLGA8A, EID2B, MYO1D, CD300LB, SRSF8, SNRNP35, SOCS3, CYP2R1, MIR22HG, DYNC2H1, NHS, ZNF559, NBR1, FAM120AOS, LIN54, SUPT3H, NIF3L1, ARID5A, EIF4BP6, MCMBP, NOL8, PRIM1, ITSN2, INPP5F, ARMCX6, HIST2H2BF, GNB2L1, CNEP1R1, TRAJ19, ERV3-1, DNAJC9, IRF9, CTD-2006C1.2, SNORD12B, AC009948.5, RP11-34E5.4, AC093495.4, YWHAZP3, RP11-37E23.5, RP11-175O19.4, PFN1P1, FTH1P2, PRKCQ-AS1, RPS2P5, RN7SL288P, RP11-259P15.4, KB-431C1.4, NPM1P27, NAIP, RP11-318K15.2, AF131215.2, MIR4519, RP11-44F14.8, RP1-178F10.3, RP11-635N19.2, SNORD64, RP11-65L3.4, UPF1, ZNHIT6, NCOA4, RP5-1041C10.3, CAAP1, MPHOSPH10, RBM17, RABGGTB, GPR107, RER1, RP11-380G5.2, RP11-3P17.3, HEATR5A, INADL, KLHL18, ATP6V0A1, U2SURP, ATXN7L3B, TMBIM1, LRRC58, CEP68, TG, DNM1L, DCUN1D4, SDHB, NR2C1, STYXL1, EFTUD1, ABCE1, CTD-2192J16.15, ZSWIM8, AC074289.1, PPFIA1, SNRPA1, ARMCX5, RBBP4, MRFAP1, TMBIM6, DNPH1, SLC16A6, CPSF6, SLC10A7, ARF3, ACTR2, LYSMD4, MARCH5, TSEN15, RP11-864N7.2, ACAD11, CTC-228N24.3, KB-431C1.5, DDHD1, OGFOD1, ZNF259, ZNF592, ZNF91, SLC22A1, NBN, ZNF586, KLF7-IT1, TAS2R30, CCDC109B, DIMT1, KLF9, CLEC1A, RPGR, POLR2H, ATG7, ZNF43, PRR13, RP11-75L1.2, RP11-747H7.3, RP11-67L2.2, SEC22C, ACVR2B, HSPE1, UBA2, RLIM, KRAS, PTPRN2, SLFN5, PCBP1, ANKRD13D, IRF7, ZMYM1, LDHBP2, DDX10P1, ILF3-AS1, CTB-152G17.6, SRGN, ATPAF1, SIGLEC10, TMEM71, RSL1D1, TMA16, MIEN1, SLC25A30, PLCB2, DLAT, FBRS, GTPBP4, ATP6V0E1, EIF4BP7, KRR1, ZNF599, HNRNPA3, QSOX1, ZNF18, GRK6, MAX, DYRK2, PTPN12, ANKAR, FKBP1C, SEC22B, HK2P1, OTUD5, POLD3, NLRC4, PTK2B, CCDC65, CCDC19, RP11-889L3.1, ZMPSTE24, ABCB7, CTD-2228K2.5, CDC42P6, LSG1, LPAR2, THUMPD1, MNT, PLXNA2, RNMT, RPL18A, INVS, MRS2, CHMP1A, OASL, HNRNPA1L2, LPP, AHCTF1, CALM3, ANO10, NBEAL2, RRM1, BMI1, C2orf68, WIPF2, MRPL14, PLAG1, C19orf59, PRPF39, BMPR2, MIR29B1, RP11-742N3.1, RP11-732M18.3, EIF2S1, SRSF2, THOC2, POGLUT1, RFX7, MAP4K3, GOPC, DGKG, MRPL28, MTG2, MDN1, DBN1, BPHL, PRPS1, SLC38A9, HIST1H4E, METTL10, RP11-44F14.2, ISG20L2, ABHD3, TMEM161B, TRUB2, LETM2, TIMM21, SUV420H1, RABL3, EIF3EP1, IFT74, CAMK2D, RP11-196G18.22, ELP6, MAU2, ACSL3, ZXDB, SEH1L, GLS, ACVR2A, PI3, ZNF710, CYP4V2, MAP1LC3B2, OXSR1, KLHL11, RPS19P1, CAMKK1, TTLL4, LINC00173, AC005795.1, WDR60, VPS45, PTPLB, AD000092.3, KLHDC2, HSPBAP1, TNPO3, TAB3, PEX3, INPP5A, TIMP1, VPS8, SENP2, AKAP13, C11orf31, RN7SL239P, PPP1R9B, TPM4, SMYD3, EPB41L4A-AS1, PPIE, CEP41, PDE4DIP, ZNF563, EXOSC7, RTF1, C12orf23, PEF1, DUSP18, NNT-AS1, MEF2BNB, ZNF175, SOD1, ZNF34, RP11-580I16.2, RAPH1, SEPT6, PTPN6, CSE1L, MARCH7, USP53, RSBN1L-AS1, MRE11A, BPNT1, HECA, EMC7, MSI2, RPL23AP79, ARNTL, ALG2, DYNC1LI1, COG5, DMXL1, COA4, MGEA5, DHRS4-AS1, hsa-mir-6080, RP11-378A13.1, RP11-396K3.1, RBM26, PHKA2, CDC42, ZBTB25, AFTPH, FGD3, WDFY4, NOM1, NPAT, CSRP2BP, MOCS2, NLRC3, TP53RK, C9orf66, RP5-874C20.3, KATNBL1P6, LRRC37A15P, YAE1D1, UBA6-AS1, VPS28, PPP2R5A, MICALCL, PPP1R21, RP11-70C1.1, TWISTNB, ATXN7, RPL10A, RIPK1, PIGK, GS1-358P8.4, TTC5, CCDC6, MED25, EPS15L1, FIGNL1, ATIC, CIRH1A, PRKCI, TRIM56, CRLF3, PIPSL, OR52K3P, ANKIB1, RPL37, UBE2J2, MRVI1-AS1, HEBP2, ZCCHC2, RASGRF2, MLH1, RPL38, BIN2P1, PCED1B-AS1, RPL3, TSPAN2, AC116366.6, SIRT5, ASUN, ANKS1A, DOPEY1, ATMIN, RBMS1P1, RPS14, DHX36, SPATA2L, SRF, RP11-418J17.1, MID2, PDCD5, FBXL20, MARCH9, CD226, OXSM, BCL2L11, SYAP1, MRPL48, FAM73A, RP11-421L21.3, RP11-100N21.1, ZNF674, RP11-76E17.3, RP4-714D9.5, CCT4, KIAA1328, TP53TG1, HDDC2, LCP2, ELOVL1, BEX4, DDX58, INPP4B, BRD8, EEPD1, MIER3, SGCB, TMEM67, ARF4, LCMT2, GCSAM, TRAJ12, ZNF512, TAS2R18, PWAR6, TRAF1, UHMK1, MKL2, TYW1, OTUB1, VPS41, C17orf75, EXOSC2, PNRC1, IMMP1L, PYGM, PTPRC, SH3KBP1, MZT2B, CENPH, SNX20, NOP10, MAP7D3, USP40, CYB5B, DNAJB6, FKBP14, PANK3, DNAJB5, NAA30, EIF2A, ASAP1, APOOL, XKR8, U2AF1, KLC2, FANCF, ZNF121, ZNF69, NR3C2, CTD-2576F9.1, RP11-290F24.6, CDC23, CYP2U1, MKS1, IQCE, FAM122B, CD300LF, RSBN1L, UROS, KIAA1279, FAM120C, SNRNP40, CASP1, ZSWIM5, ALG10B, ZNF329, UPF2, DENND2D, NPM1P6, RP11-809O17.1, MTPN, SLC2A3P2, SLC35A3, CTSD, VIPAS39, PRIMPOL, GSR, CDC14A, SLBP, DNAJC2, ITM2A, RPL18, LINC00630, IMMP2L, PPM1H, AL162151.3, RP11-7F17.3, SEPHS1, NOTCH2, OSCAR, SSBP1, ZNF7, COMMD1, EVL, ABCF2, LONP2, CCDC64, RPL32, AC006129.2, STX17, CHCHD3, RPS20, PGM3, MAST4, METTL8, MMS22L, MSRB2, N6AMT1, TPCN2, MEI1, RPSA, TIGD3, ATG13, RPL14, SEPT7P2, KATNBL1, SMC2, UBALD1, THAP7, ESD, RP11-452F19.3, CCDC132, NUDT16, CRLS1, GNPNAT1, MRPS9, ALDH9A1, AASDH, TMEM42, TMX2, RPL13AP5, RNASEH2B, POLR3C, RP11-177C12.1, KIF5C, PURA, C15orf41, WDR59, PPP3CB, ARAP3, MGST3, CCT2, PRELID1, CLP1, ASB8, NAP1L2, ZNF569, ZXDA, RP11-68I3.11, GNG5P2, DTD2, CTC-479C5.12, CASP6, MLLT3, TUBB4B, ZNF568, RUFY3, SUDS3, USP15, KRBOX4, MIOS, NDUFB6, CTD-3088G3.8, CTBP2P8, LRRFIP1P1, RBBP4P1, RALA, KBTBD4, RALY, PHF20L1, ZCCHC11, SRFBP1, ZNF83, BOLA1, RPL23A, LL22NC03-86G7.1, ARID4A, PLD1, LRRK1, COQ4, IDS, CCNB1IP1, ORMDL1, RP5-1115A15.1, ANKRD44, ZNF827, LRRN1, RP11-173M1.8, KNOP1, EMG1, ZC3H11B, C18orf25, UBR1, PSMD6-AS2, AC093323.3, GBA2, MEAF6, RPL37P2, ZMAT1, NAF1, ZNF558, MTRR, RP11-641D5.1, ZNF33B, MANBA, KIAA0319L, RP11-553L6.5, SH3YL1, CEP128, RPL37A, ADH5P4, UFL1, CSPP1, PCED1A, RYK, SP3, EWSR1, ZFP1, SBDSP1, PRR14L, STK3, CHP1, MTX3, ZNF550, REXO2, NBAS, SASH3, NIT1, ZNF33B, AC007969.5, SMURF1, PSMA5, VAPA, ANKEF1, PAN2, NCSTN, PLEKHA2, RCSD1, CHUK, RP11-196G18.24, DGCR14, DDX50, EEA1, ULK1, KCNQ1, STK10, UBA5, MTMR2, CRNKL1, ELF4, STX2, TMEM30A, MOB1A, RPS15, CRYZ, PROSER1, ZNF706, HVCN1, SRSF6, RIOK1, ATXN10, PCBD2, C12orf29, LMO2, TGS1, SCAPER, TMEM50B, TMEM87B, EIF4A2, WIPI2, SAFB, ZC3H12A, MFSD8, DCBLD1, FBXO22, FAM21C, SCAI, PTP4A3, RPL12, ZNF607, RP11-344B5.2, RP11-17A4.1, RP11-902B17.1, PTGES2, TNFRSF10A, SCAP, RPL21, CHD6, FUNDC2, KNTC1, CRYZL1, RP13-977J11.5, DGCR11, MAN2B2, ALG9, HBS1L, UBE2L6, C6orf89, UNKL, ODF2L, COMMD10, GPN1, NOV, FAHD2A, HUS1, ATP7A, ZBTB20, AK3, STK35, NEDD1, STMN3, TRBV29-1, MPP6, HSPH1, ZSCAN20, SLC30A7, TRAV17, CRKL, WIPF1, TIGD6, PSMG1, NUP88, ZNF546, ITGB2, TMEM41A, THNSL1, CTD-2270N23.1, EXTL2, ZNF571, ACVR1, ARAP1, C17orf103, ASTN2, ATP2C1, ERP29, FAM178A, IL11RA, DTX3L, RP11-74E24.2, RC3H1-IT1, POLR1B, RRAGB, FLCN, ZNF146, LSM10, RP11-137L10.6, TAOK1, RP11-51F16.5, GNL3, ZFYVE1, WBP11, C12orf65, RPL7, ZDHHC21, RP11-1319K7.1, C7orf41, RPS4X, RPL37AP1, RPL7AP11, RPS27, DYRK4, SLC39A8, OSBPL1A, GUF1, ZNF846, GIMAP7, TBC1D22A, B4GALT1, XIAP, VPS18, NDRG1, MAP2K6, BNIP2, KDELC2, PSMB10, ZNF880, ZNF529, PDAP1, ACVR1C, ACTN4, ADPRHL2, PM20D2, CUL5, WBP2, CYB5A, OFD1, CARD8, COPS8P2, GTDC1, RP11-16E12.1, ZNF570, RAN, CTSZ, PPAT, MAP1S, DBT, IBTK, TMEM206, NDUFAF5, RPS12, ZNF22, ZNF792, FMNL1, CTD-2287O16.1, ZNF451, ZFP62, RNF111, HDAC5, NUP153, ZBTB5, RP11-711M9.1, RP11-3K24.1, FXN, MGAT1, CDK8, SUCLA2, SIAH1, RP11-159C21.4, SPAG16, TNFRSF25, CTD-2528L19.6, PLOD1, TM9SF2, EI24, DCK, NUDT9, CTC-575D19.1, AC073046.25, DNAL1, CCDC25, CCDC28A, NELL2, BHLHB9, HAUS2, GRAMD3, ANGEL1, PBDC1, C4orf21, RPL14P1, ARL8A, EBAG9, RPRD2, PPAPDC2, MRPL11, PCIF1, SEC61A2, COPS8, RPL3P4, POMGNT1, C10orf88, FAR2, CDR2, BBC3, PSD4, RPL13A, ZNF326, RPN1, PRMT10, LGALS9, ST7L, EZH1, KLF16, MTERFD1, RPLP0P6, SH3RF3, WDR92, SEC63, GDPD1, TAF9B, RAB21, CPOX, PLP2, HTRA2, ZNF549, TESPA1, TMLHE, UBE2M, HNRNPK, GRSF1, EIF3L, GRAMD1C, MIR3671, ARID4B, GLTSCR2, ADAP1, STAT5A, MBD6, C15orf61, PPIAP22, RP11-326I11.3, ZNF793, TRMT6, RN7SL172P, SNX1, SLC2A9, STK38, CHRNA10, SLC39A10, METTL6, SCMH1, RIC8B, PPP2R1B, RNMTL1, CD2AP, IL23A, DNAJC14, MRPL42, BCCIP, SPATA13, RP3-354N19.3, KPNA5, ZC3H4, ZNF280D, KANSL2, VPS39, NACA, NUFIP2, DARS, AC007038.7, RP11-351I24.3, GOLGB1, AF131215.9, ATP6V1E1, RPLP0, DPP8, METTL2A, ACIN1, RPGRIP1L, MKLN1, CCDC136, USP3, OCIAD2, MSN, TMEM192, ZBTB14, SMIM11, RPL15P3, CD302, RP11-284N8.3, SPTLC1, RP11-314A20.1, MRPL43, ACTB, SENP1, SAYSD1, CCNG1, DYNC2LI1, MTMR6, THEMIS, RP11-171I2.2, PRKD3, MSANTD3, MTIF2, RASGEF1B, ZNF461, BRAP, ANKRD40, DCAF16, AC000089.3, CPEB2, NDUFA9, MEX3C, EGLN2, EXOSC6, TP73-AS1, NIPA2, PVRL3, SKP2, SMIM19, CRY1, UBC, ORC3, IL16, PRKRA, COL18A1, DCAF13, DTWD2, TDRKH, EIF5A2, WWOX, FAM217B, PRPF4B, SLK, COX11, ZNF257, RPS20P14, MSRA, ACAT2, SYTL3, TRIM35, HIST1H4C, TMEM19, CCDC167, NME3, ZNF714, ZNHIT3, APH1B, SNHG6, NEK1, HNRNPKP2, PRR14, PDP2, FYTTD1, SLC30A6, FAM208A, ORC5, RPL13P12, SKIV2L2, STAM2, NR1H2, TMEM18, CSTF3, RPSAP58, TMEM242, TIA1, HSPB11, PIGM, FEM1C, CALCOCO1, WDR83OS, PSMD14, WDR26, ADAM22, MRPS35, ROCK1, SELO, CLUAP1, CCDC9, PPWD1, VAMP4, TSFM, AVIL, METTL5, PPA2, SCAF11, C21orf59, ZNF382, HELQ, ZNF143, ADAM9, ZNF540, FAM91A1, HIST2H2BE, ADAT2, ZNF33A, GAS5, RP11-557J10.3, RP11-658F2.8, ZNF195, SPDL1, USP13, MAPRE1, MKI67IP, GAA, SMC5, RPL4P4, RP11-384C4.6, SP4, RNF170, FAM135A, LIG1, CTD-2302E22.5, SNAPC3, LGALSL, RPS3, PIM2, ZNF528, TXNDC12, UXT, VPS13B, IRAK1BP1, CAMK1D, NOL12, RP3-375P9.2, GPA33, TATDN1, ACSL6, AC009961.3, RP11-58H15.1, SRSF4, ARL1, RPAIN, PARP8, ZNF567, CTB-31O20.2, CEP76, RPS25, RSPH3, PINK1, BRIX1, ZNF302, ZC3HAV1, UBE2J1, LRRC40, RWDD1, DHX8, MANBAL, MRRF, EIF4BP3, RPL7AP6, PDCD2, MTBP, TRAJ17, PUS3, POT1, SLC35E3, MRPL20, SAP30L, PSMG3, DZIP3, STYX, LRPAP1, CCDC174, PRR4, ZNF736, RAB5B, TMEM161B-AS1, RHOBTB3, IDE, IKZF2, CDH26, RPS4XP11, RGS3, HSCB, RP11-640M9.2, SC5D, RPF2, GLO1, RP11-16F15.1, STK24, SLPI, CASC7, OCIAD1, MIR4426, ALG10, PCGF6, INTS6P1, TBCCD1, ZNF227, MUTYH, PPP2R2D, ZNF584, TRAPPC2L, TPST1, ABCD2, RPS23P8, ABCD3, EIF3M, GEMIN6, GCHFR, LARP1B, GSTCD, NT5C3B, NUDCD2, TRNAU1AP, ZNF544, RPL17, ZNF17, RP11-802E16.3, ATP6V1C1, TYW3, C17orf89, DBF4, FAM84B, MTMR8, ADCK4, UST, DANCR, RP11-10J21.4, PARP16, RHOQP2, YWHAZP2, RPL4, UPF3AP3, SNHG5, PKNOX1, CTD-2373J19.1, UBE2MP1, RABEP1, MAGI3, KIAA1377, RBBP5, UTP18, MINPP1, CENPV, MARS, CTD-2341M24.1, ZYG11B, RPS29, ARRDC1, ZNF624, TMEM260, DNAJC10, LACE1, POLE3, ZNF131, UBE2N, KDM5B, TXNIP, EBNA1BP2, NBPF10, ZSWIM7, COX10-AS1, KLHL8, AGO2, USP9X, NUPL2, FCHO2, CCNL1, MSANTD4, RPL23AP2, WDR86-AS1, STXBP3, ACSS2, PLK3, RP1-197B17.3, LYPLA2, MPHOSPH6, HAX1, STXBP4, GAPVD1, SS18L1, PMS2P1, MRPS25, CHCHD6, GRN, NEK3, RPL10AP2, FAM53B, ASAH1, FASTKD1, AIMP1, PGK1, PRKD2, CBX4, PLEKHM3, ZNF485, AATF, RNLS, CDK5RAP1, HEXIM1, PGAP1, EPHA4, MCM8, RP5-1112D6.4, CORO7, NANP, ZNF10, TAB2, EARS2, TSGA10, TMEM261, DNAJC21, MTRF1, TRIM32, TGOLN2, CCDC104, FLII, NSUN6, CLDND1, TRAF6, ZNF510, ACADM, STX18, RP5-1085F17.3, RPL41, TMCO4, SERTAD1, COX16, ELAC1, SNIP1, PHB2, TXNDC15, TRMT13, APPL1, RPS2P46

### 16.7 Proposal tree — direct-effect genes not already listed

**DIRECT on the proposal tree only** (3)

FAM107B, LDHA, SLC46A3

---

## 17. Files

| Path | Holds |
|:--|:--|
| PRESPEC.md | the protocol, fixed before each analysis ran (incl. addenda A and B) |
| METHODS.md | methods and the incident log |
| DIVERSITY_RESEARCH.md | the diversity research, with sources |
| MEDIATION_MASTER.md / .html | this document |
| scripts/med.R | the estimator — one implementation used everywhere |
| scripts/00_validate_med.R | 27 correctness checks |
| scripts/02_effects.R | the primary result, both trees |
| scripts/03_bootstrap.R, 04_permutation.R, 04b_perm_pm.R | bootstrap and permutation inference |
| scripts/05–07, 09 | sensitivity, calibration, E-value, power, exploratory |
| scripts/10_circularity.R, 10b_… | circularity sensitivity analysis |
| scripts/11–13 | diversity: prep, alpha, beta |
| scripts/14_gene_lists.R | gene categories |
| results/genes/ | every gene with its category and statistics (CSV) |
| results/circularity/ | circularity results |
| results/diversity/ | diversity results |

## 18. Key references

- Pearl J (2001). Direct and indirect effects. *UAI*.
- VanderWeele TJ, Vansteelandt S (2014). Mediation analysis with multiple mediators. *Epidemiologic Methods* 2:95–115.
- VanderWeele TJ (2014). A unification of mediation and interaction: a 4-way decomposition. *Epidemiology* 25:749–761.
- Valeri L, Lin X, VanderWeele TJ (2014). Mediation analysis when a continuous mediator is measured with error. *Statistics in Medicine* 33:4875–4890.
- VanderWeele TJ, Ding P (2017). Sensitivity analysis in observational research: introducing the E-value. *Annals of Internal Medicine* 167:268–274.
- Egozcue JJ et al. (2003). Isometric logratio transformations for compositional data analysis. *Mathematical Geology* 35:279–300.
- Newman AM et al. (2015, 2019). CIBERSORT; CIBERSORTx. *Nature Methods* 12:453–457; *Nature Biotechnology* 37:773–782.
- Krebs CE et al. (2020). Whole blood transcriptome analysis in bipolar disorder reveals strong lithium effect. *Psychological Medicine* 50:2575–2586.
- Martínez O, Reyes-Valdés MH (2008). Defining diversity, specialization, and gene specificity in transcriptomes through information theory. *PNAS* 105:9709–9714.
- García-Nieto PE, Wang B, Fraser HB (2022). Transcriptome diversity is a systematic source of variation in RNA-sequencing data. *PLoS Computational Biology* 18:e1009939.
- Jost L (2007). Partitioning diversity into independent alpha and beta components. *Ecology* 88:2427–2439.
- Anderson MJ (2001, 2006); Anderson MJ, Walsh DCI (2013). PERMANOVA; PERMDISP; heterogeneous dispersions. *Austral Ecology*; *Biometrics*; *Ecological Monographs*.
- Zapala MA, Schork NJ (2006). Multivariate regression analysis of distance matrices for testing associations between gene expression patterns and related variables. *PNAS* 103:19430–19435.
- Liberzon A et al. (2015). The MSigDB hallmark gene set collection. *Cell Systems* 1:417–425.
- Zhao X, Lynch JG, Chen Q (2010). Reconsidering Baron and Kenny: myths and truths about mediation analysis. *Journal of Consumer Research* 37:197–206.
- Full lists with verification marks: `DIVERSITY_RESEARCH.md` §10 and `METHODS.md`.

