# Mediation arm — analysis protocol

**Status: fixed before any mediation estimate was computed.**
Written 2026-10-08, immediately after `01_prep.R` and immediately before `02_effects.R`.

**The governing document is `Implementation/Proposal Feedback assignment/AliAva_Proposal.docx`**
(the full proposal, 2026-10-07). Its "Mediation Arm" section specifies the primary model,
and this protocol implements it. Where the proposal and this document differ, the
difference is named and justified in section 6.

This document also exists because of a direct instruction in the pre-proposal feedback:

> "Swapping exposures, mediators and outcomes until something is significant is a
> multiple-testing problem. **Pre-specify the primary mediation model now.**"

**On the source of the quotations below.** They are the 39 comment annotations embedded in
`Syllabus/AliAva_PreProposal_Feedback.pdf`. The PDF's annotation author field records
`ahmetay` for all of them, but Ali has said they may be the team's own review notes rather
than the instructor's. They are therefore cited here neutrally as **"review notes"**, with
no claim about authorship. Nothing in the analysis depends on who wrote them — they are
used as a checklist of methodological obligations, and each one is either met or explicitly
declined below.

Everything below is declared in advance. Anything not in this document, if it appears in
the results, is labelled post hoc there.

---

## 0. What was already known when this was written

Full disclosure, because this is not a true pre-registration and should not be called one:

* The differential-expression arm had already run. Lithium has a large whole-blood
  signature (1,426 genes at BH 5 %; a 1,092-gene core robust across 31 filter settings).
  Bipolar-off-lithium vs control has no signature surviving any correction.
* Lithium's DE signature was already shown to be *largely compositional*: 92.8 % direction
  agreement with the granulocyte fraction, 40.8 % median shrinkage on adjustment.
* Reliability of the mediator was already measured (`00_reliability.R`):
  b1 0.720, b2 0.585, b3 0.642, b4 0.516.
* **No mediation quantity — no NDE, NIE, proportion mediated, or interaction test — had
  been computed in any form.**

So: a large total effect and a compositional flavour were expected. The *size of the
mediated fraction*, the *interaction*, and every number in the results were not.

---

## 1. Primary analysis — one model, declared now

| | |
|:--|:--|
| **Sample** | the 226 bipolar I subjects (152 lithium users, 74 non-users) |
| **Exposure** `X` | lithium use, 0/1 |
| **Mediator** `M` | the **four ILR balances, jointly, as one vector-valued mediator** |
| **Outcome** `Y` | voom log-CPM, **one model per gene**, on the **12,018 non-LM22 genes** |
| **Covariates** `C` | age, sex, tobacco, RIN, plate, seqPC1-3 (assessment group is constant here and is dropped) |
| **Estimator** | regression-based natural effects **with exposure-mediator interaction** |
| **Primary estimand** | the **joint natural indirect effect (NIE)** through the whole mediator vector |
| **Primary inference** | delta-method SE -> z -> **Benjamini-Hochberg at 5 % across the 12,018 genes** |
| **Headline summary** | **median proportion mediated among genes with a significant total effect** |

### 1.1 The models

Mediator models, k = 1..4, sharing one design matrix:

```
M_k  =  b0_k  +  b1_k * X  +  b2_k' * C  +  e_k
```

Outcome model, per gene g:

```
Y  =  t0  +  t1 * X  +  SUM_k t2_k * M_k  +  SUM_k t3_k * (X * M_k)  +  t4' * C  +  e
```

`M` and `C` are centred at their sample means, so the model is evaluated at the mean
covariate profile and the intercepts are directly interpretable. Centring provably does
not change NDE, NIE or TE; it only improves conditioning. `00_validate_med.R` checks this
numerically rather than asserting it.

### 1.2 The estimands

```
NDE  =  t1  +  SUM_k t3_k * b0_k                  natural direct effect
NIE  =  SUM_k (t2_k + t3_k) * b1_k                natural indirect effect  <- PRIMARY
TE   =  NDE + NIE                                 total effect
PM   =  NIE / TE                                  proportion mediated
```

These are Pearl's (2001) natural effects, after Robins & Greenland (1992), in the
regression form of VanderWeele & Vansteelandt (2009, 2010), extended to a **vector**
mediator by VanderWeele & Vansteelandt (2014) — the reference already in the proposal.

Because the mediator is a vector, the quantity above is the **joint** indirect effect
through the whole cell mixture. Splitting it into path-specific effects through individual
balances would require an assumed causal ordering among the balances, which we do not have.
Per-balance terms are therefore reported as a **descriptive decomposition of the joint
effect, not as four separate causal effects**, and the report says so.

Also reported, the **four-way decomposition** (VanderWeele 2014), which splits TE into the
part needing neither mediation nor interaction, the part needing only interaction, the part
needing both, and the part needing only mediation:

```
TE = CDE + INT_ref + INT_med + PIE
CDE     = t1                        neither
INT_ref = SUM_k t3_k * b0_k         interaction only
INT_med = SUM_k t3_k * b1_k         mediated interaction
PIE     = SUM_k t2_k * b1_k         pure indirect
```

### 1.3 Why not product-of-coefficients, which is what the proposal's flowchart says

Because the feedback told us not to:

> "Mediation flowchart: product of coefficients assumes no lithium-by-cell-mix interaction.
> **Test that interaction before reporting a proportion mediated.**"

The estimator above **contains** product-of-coefficients as the special case `t3 = 0`, so
nothing from the proposal is abandoned. We additionally:

* run a **joint test of H0: t3_1 = t3_2 = t3_3 = t3_4 = 0** per gene (moderated F, limma),
  BH-corrected, and report how many genes show interaction; and
* report **both** estimators side by side, so the cost of the proposal's original
  assumption is a number rather than an argument.

### 1.4 Thresholds, declared now

| Quantity | Threshold | Source |
|:--|:--|:--|
| "substantial" proportion mediated | **PM > 0.50** | the feedback asked for a number and suggested this one |
| gene-level significance | **BH q < 0.05** | proposal flowchart |
| interaction present | **BH q < 0.05 on the joint 4-df test** | feedback |
| bootstrap replicates | **B = 1000**, subject-level resampling | |
| permutation replicates | **B = 1000**, lithium labels shuffled | proposal flowchart |
| primary CI | **percentile bootstrap, 95 %** | |
| seed | **481**, threaded everywhere | CLAUDE.md 4.5 |

PM is computed **only for genes whose total effect is significant**. A ratio whose
denominator is indistinguishable from zero is not interpretable, and reporting it would
manufacture spurious values of plus or minus infinity.

### 1.5 Direction of the primary result, declared now

The proposal predicts lithium raises **myeloid relative to lymphoid only** — that is,
`b1_1 > 0` and `b1_2 = b1_3 = b1_4 = 0`. The feedback challenged the word "only":

> "Why 'only'? Any other lineage shift would then falsify the hypothesis; is that intended?"

We therefore test **all four** balances, report all four, and treat a shift in b2, b3 or b4
as a finding about the hypothesis rather than a nuisance.

---

## 2. Pre-declared sensitivity analyses

Every one of these is declared now; none is contingent on the primary result.

| # | Analysis | Why |
|:--|:--|:--|
| S1 | **LM22 marker genes (n = 350) analysed separately** | CLAUDE.md 4.2. The mediator is estimated *from* the outcome matrix; for signature genes this is near-tautological. The gap between the 350 and the 12,018 *measures* the circularity. |
| S2 | **Second, independent deconvolution** — Krebs et al.'s deposited CIBERSORT fractions put through the identical lineage/ILR construction | feedback: "Rerun with different ways to find the cell fractions." |
| S3 | **Regression calibration** for mediator measurement error, using the per-balance ICCs, following Valeri, Lin & VanderWeele (2014) | the proposal names this method and this reference, and says to report the uncorrected figure **as a lower bound**. |
| S4 | **E-value** for the total and indirect effects, benchmarked against the obesity-to-neutrophil association | the proposal names this benchmark explicitly (Furuncuoglu et al. 2016), because bodyweight is not recorded in GSE124326. |
| S5 | **All 20 tobacco imputations, pooled by Rubin's rules** | review notes asked which imputation technique. |
| S6 | **Zero-replacement delta sensitivity**, delta in {0.1, 0.3, 0.5, 0.65, 0.8} | feedback: "does that choice change the ratios?" |
| S11 | **The proposal's own 6-part / 5-balance ILR tree**, run end to end as a parallel primary | the proposal text specifies this tree; the implementation uses a 5-lineage / 4-balance tree. See section 6. |
| S7 | **Alternative multiple-testing corrections**: BH, BY, Bonferroni, Holm, Storey q | consistency with the DE arm. |
| S8 | **Product-of-coefficients (no-interaction) estimator**, reported alongside | feedback. |
| S9 | **Unadjusted and age-only-adjusted fits** | age is the one covariate imbalanced with lithium (p = 0.0007). |
| S10 | **Power / minimum detectable NIE at n = 226** | feedback: "Add a power calculation for the mediation arm with 226 patients." |

## 3. Pre-declared secondary analyses — labelled secondary, not primary

The proposal says: *"Will do this mediation analysis a couple of different ways, changing
the exposure to maybe age only or trying out different combinations."* The feedback
warned that doing this freely is a multiple-testing problem. So the list is closed here:

| # | Analysis | Source |
|:--|:--|:--|
| T1 | **Age as exposure**, cell mix as mediator, same outcome set | proposal flowchart note |
| T2 | **Moderated mediation**: does the age-mediated effect differ between lithium users and non-users? | proposal hypothesis + Salarda et al. (2021) |

T2 is reported **with its power**, because the feedback asked directly:

> "This is a separate hypothesis about moderated mediation. With 74 non-users, do you have
> power to detect a difference in mediated effects?"

No further exposure/mediator/outcome combination will be run. If one is wanted later it is
a new, separately declared analysis.

## 4. Assumptions, stated plainly

Identification of natural effects needs (Pearl 2001; Imai, Keele & Yamamoto 2010):

1. no unmeasured exposure-outcome confounding, given C;
2. no unmeasured mediator-outcome confounding, given C and X;
3. no unmeasured exposure-mediator confounding, given C;
4. no mediator-outcome confounder that is itself affected by the exposure.

**(2) and (4) are not guaranteed here and cannot be made so.** Illness severity, mood state
and other psychotropics plausibly affect lithium prescription, cell mix and expression
together, and **none of them is in this dataset**. BMI is not recorded either — checked, it
does not exist in GSE124326 — so it is bounded rather than adjusted, exactly as the
feedback anticipated. This is why S4 (E-value) is not optional.

Every estimate here is **associational under a causal interpretation**, and the write-up
says so.

## 5. What would falsify the hypothesis

| Result | Reading |
|:--|:--|
| median PM > 0.5 among genes with significant TE, and b1 the balance that moves | hypothesis supported |
| median PM <= 0.5, direct effect dominant | hypothesis **not** supported |
| balances other than b1 shift | the "only" in the proposal is wrong; report it |
| PM in the 350 marker genes >> PM in the 12,018 | circularity is inflating the estimate; the 12,018 figure is the defensible one |

---

## 6. Where this protocol departs from the proposal text, and why

Three departures. All are recorded here rather than quietly absorbed.

### 6.1 The mediator is 5 lineages / 4 balances, not 6 cell types / 5 balances

The proposal says:

> "five balances follow the hematopoietic lineage: myeloid {neutrophils, monocytes} vs
> lymphoid {CD4 naive, CD4 memory, NK, B}; then neutrophils vs monocytes; T {CD4 naive,
> CD4 memory} vs {NK, B}; CD4 naive vs CD4 memory; and NK vs B cells"

The implemented mediator instead aggregates all 22 LM22 types into **5 lineages**
(granulocyte, monocyte, T, NK, B), giving **4 balances**. The reason is zero-inflation,
measured rather than assumed:

| | 22 types | 6-part (proposal) | 5 lineages (implemented) |
|:--|--:|--:|--:|
| cells that are exactly zero | 4,811 / 10,428 = **46.1 %** | 72 / 2,844 = **2.5 %** | 1 / 2,370 = **0.042 %** |
| samples needing zero replacement | nearly all | 71 / 474 = **15.0 %** | 0 |
| share of total cell mass used | 100 % | **92.4 %** | 100 % |
| median reliability (ICC vs a second deconvolution) | 0.389 | not measured | **0.693** |

Aggregating to 5 lineages nearly doubled the reliability of the mediator (0.389 to 0.693).
A log-ratio built on a part that is zero in 100 % of samples (NK activated) or 98 % of
samples (B memory) is not a measurement of anything.

**This does not make the proposal's tree unusable, and we do not discard it.** It is
computable, so S11 runs the *entire primary analysis* a second time on the proposal's exact
6-part tree, and the report shows both side by side. If they agree, the choice did not
matter. If they disagree, the reader sees the disagreement rather than our preference.

**Action item for the written proposal**: the sentence above should be corrected to describe
the tree actually used, or the 6-part result promoted, once S11 has run.

### 6.2 A factual error in the proposal that should be corrected

The proposal states:

> "we found that even though the cell type that were estimated using different
> deconvolution techniques were different from each other and from the ground truth, the
> ILR ratios that were computed were more accurate/consistent across different
> deconvolution techniques"

**This is not what the data show.** `00_reliability.R` measured it directly against a second,
independent deconvolution (Krebs et al.'s deposited CIBERSORT fractions):

| | reliability (ICC, absolute agreement) |
|:--|--:|
| 22 raw LM22 proportions | 0.389 |
| 5 aggregated lineage proportions | **0.693** |
| the 4 ILR balances | 0.720, 0.585, 0.642, 0.516 (median **0.613**) |

The ILR transform did **not** improve consistency across deconvolution methods; it is
slightly *worse* than the lineage proportions it is built from. **Aggregation** from 22
types to 5 lineages is what improved consistency. The correct justification for ILR is the
one the proposal already gives elsewhere and which stands on its own: the proportions sum
to one, so they cannot all enter a regression, and a rise in one forces a fall in the
others. That argument needs no reliability claim.

### 6.3 "Age as exposure" is a contingency in the proposal, and stays a contingency here

The proposal lists it under the work plan as a fallback:

> "if we run into any issues such as completely uninformative/uninterpretable results, we
> will add a variable such as age along with the Lithium so the input is continuous instead
> of binary"

T1 and T2 in section 3 are therefore run and reported as **exploratory**, exactly as the
proposal's own rule requires ("We specify a single primary model here and designate every
other configuration exploratory"). They are run regardless of whether the primary result is
interpretable, because running them only on a bad result would make them a contingency on
the outcome -- which is the multiple-testing problem the feedback warned about.


---
---

## 7. ADDENDUM A — circularity sensitivity analysis

**Declared 2026-10-18, before any of it was run.**
Requested by Ali: *"add another mediation that is done on only the 547 signature genes, and
then on the entire filtered gene set including the 547 signature genes, so that we can do a
sort of circularity sensitivity analysis."*

### 7.1 Which of the 547 signature genes exist in this data — measured, not assumed

| Step | Genes | Note |
|:--|--:|:--|
| LM22 signature matrix | **547** | Newman et al. 2015 |
| present in the GENCODE v19 count matrix | **541** | absent: FLJ13197, LINC00597, LINC00921, LOC100130100, LOC126987, TRBC1 |
| pass the expression filter (>10 counts in ≥90 % of samples) | **350** | the genes that carry signal in blood |
| present but fail the filter | **191** | median gene has >10 counts in only **3.8 %** of samples; they mostly mark cell types that are rare or absent in blood (Tregs 25, plasma cells 18, macrophages M1 17, mast cells activated 15, dendritic cells resting 13, …) |
| all-zero in all 474 samples | 2 | cannot be modelled at all |

So "the 547 signature genes" can mean, in blood, the **350 that are expressed**. The 191 others
contributed almost nothing to the deconvolution either: CIBERSORT fits on the linear scale,
where genes near zero carry essentially no weight.

### 7.2 Design — one common fit, so only the gene set differs

| Element | Choice |
|:--|:--|
| subjects | the 226 bipolar I subjects, as the primary |
| normalisation | TMM factors and library sizes computed **once**, on all 12,368 filtered genes |
| precision weights | voom **once**, on all 12,368 |
| variance moderation | empirical Bayes **once**, across all 12,368 |
| mediator, estimator, covariates | identical to the primary |
| trees | both (lineage 4-balance and proposal 5-balance) |

| Code | Outcome gene set | Genes |
|:--|:--|--:|
| **C1** | non-marker genes | 12,018 |
| **C2** | expressed LM22 signature genes | 350 |
| **C3** | entire filtered set, markers included | 12,368 |
| C4 | every LM22 gene with any counts (adds the 189 unexpressed ones via a second voom on 12,557 genes) | 539 |

Multiple testing: BH **within each set's own family** (that is what "a mediation on only these
genes" means); C2's genes are also reported within the full 12,368 family.

Inference: bootstrap B = 1000 over subjects for the median proportion mediated of C1–C3,
taken from the **same** replicates so the difference C2 − C1 gets its own paired interval;
permutation null B = 300 per set for the median PM on a fixed gene set (as in `04b`); point
estimates and delta-method only for C4.

### 7.3 Why one common fit — and why this supersedes S1

The earlier marker check (S1, `05_sensitivity.R`) normalised the 350 markers **on their own
counts**: library size = the sum of the marker counts. That sum itself rises and falls with
cell composition, so each marker was measured relative to the other markers, and composition
entered the denominator. That is a confound, not a sensitivity analysis. It is superseded.

Normalising on the whole expressed transcriptome is the standard and is what the DE arm did.
With normalisation, weights and the variance prior held fixed, any difference between C1, C2
and C3 comes **only from which genes are summarised** — the cleanest form of the question.

The primary analysis normalised on its own 12,018 genes. C1 repeats it on the 12,368 basis,
so the effect of that choice is measured, not assumed.

### 7.4 What would indicate circularity — fixed now

1. **The direct signature.** Per gene, the share of expression variance explained by the cell
   mix beyond lithium and covariates (partial R² of the balance block). Circularity predicts it
   is much larger for markers than non-markers. Tested by Wilcoxon rank-sum.
2. C2's median proportion mediated above C1's, with the paired bootstrap 95 % interval of the
   difference excluding zero.
3. C3's median PM pulled toward C2's.

If (1) holds but (2) does not, circularity inflates the mediator–outcome link without
inflating the *share* mediated — the pattern S1 hinted at, now tested properly.

---

## 8. ADDENDUM B — alpha and beta diversity, mediated through cell mix

**Declared 2026-10-18, before any of it was run.** Research: `DIVERSITY_RESEARCH.md`.
Ali delegated every open choice ("make these decisions on your own … make sure it is
justifiable"). Each is recorded with its reason.

### 8.1 Decisions

| # | Question | Decision | Why — briefly |
|:--|:--|:--|:--|
| 1 | Which meaning of alpha / beta? | **All three, ranked.** (A) per-person / between-person = **primary**. (B) within- vs between-cell-type decomposition = secondary. (C) diversity of the cell mix itself = secondary | (A) is the standard omics usage and answers the question with one test. (B) is the only exact split of bulk diversity into "cells changed" vs "mix changed" (Jost 2007), but rests on noisy bMIND estimates, so it cannot be primary. (C) is cheap and explains *why* (B) behaves as it does: the between-cell-type term can never exceed the diversity of the cell mix |
| 2 | Beta: location, dispersion, or both? | **Both.** Mediation of the location effect is primary; dispersion is mediated only if lithium changes dispersion at α = 0.05 | PERMANOVA is unreliable under unequal dispersion in **unbalanced** designs (Anderson & Walsh 2013) and ours is 152 : 74, so dispersion must be checked anyway; dispersion is also its own question (differential variability, Ho et al. 2008) |
| 3 | Exposure | **Lithium = primary.** Bipolar-off-lithium vs control = **exploratory** | Lithium is the arm's pre-specified primary. The pre-proposal feedback warned that swapping exposures is a multiple-testing problem, so bipolar is labelled exploratory and reported as such. A global measure is the natural follow-up to a per-gene DE that found almost nothing |
| 4 | Pathway database | **MSigDB Hallmark** (50 sets) primary; **KEGG legacy** sensitivity | Hallmark was built specifically to remove redundancy (Liberzon et al. 2015), and 50 programmes are interpretable. KEGG legacy is flat (not nested), so fractional counting stays meaningful; Reactome and GO are deeply nested and would make it meaningless |
| 5 | Pathway representation | **Aggregated abundance** (a gene in *k* sets gives 1/*k* of its TPM to each), not GSVA | Diversity needs non-negative abundances; GSVA scores can be negative and no published work takes the entropy of them. Aggregation is the direct analogue of collapsing species into higher taxa, which is the framework the request comes from. It also avoids a large Bioconductor dependency |
| 6 | Alpha measure | **Shannon H** (Hill q = 1) primary; Hill q = 0 and q = 2 secondary | Shannon is the measure that partitions exactly (Jost 2007) and is insensitive to depth, unlike richness |
| 7 | Gene universe | **protein-coding ∩ filtered ∩ non-LM22** primary | Structural non-coding RNAs dominate the counts (RN7SL1 6.86 %, RN7SL2 4.79 %, RN7SK 3.45 %); markers excluded for circularity, matching the primary rule. Sensitivities: markers included; all genes; mitochondrial genes excluded |
| 8 | Scale | **TPM** primary; read counts sensitivity | TPM approximates transcript proportions; read counts match García-Nieto et al. 2022 |
| 9 | Covariates | the arm's set **+ log library size** | Transcriptome diversity is a known technical variable (García-Nieto et al. 2022) |
| 10 | Beta distance | **Aitchison** primary; Jensen–Shannon and Bray–Curtis secondary | RNA-seq is compositional; Aitchison is the same geometry as the ILR mediator |
| 11 | Mediator tree | lineage 4-balance primary; proposal 5-balance sensitivity | as the arm |
| 12 | Evidence for mediation | **the indirect-effect test**, never the proportion mediated | the arm showed PM sits near 0.4 under pure noise (permutation p = 0.49) |
| 13 | Resampling | alpha: B = 10,000 bootstrap and 10,000 permutations; beta: 9,999 permutations, B = 1,000 bootstrap | alpha is one outcome per person, so this is cheap |
| 14 | Multiplicity | **primary family = 2 tests** (alpha NIE, beta mediation; lithium) → Holm. All else descriptive or exploratory | keeps the confirmatory claim small |
| 15 | Packages | install `vegan` and `msigdbr` only | `GSVA` is not needed (decision 5); `energy` (MODIMA) is not needed because MODIMA cannot adjust for confounders and age is imbalanced here |

### 8.2 The tests, exactly

| # | Question | Statistic | Null | Inference |
|:--|:--|:--|:--|:--|
| **A1** | Does lithium change alpha diversity? | total effect of lithium on Shannon H | no effect | t-test in `H ~ lithium + C`; bootstrap CI |
| **A2** | **Is that mediated by cell mix?** | natural indirect effect, `med.R` | NIE = 0 | percentile bootstrap 95 % CI (B = 10,000) + permutation (B = 10,000) |
| A3 | Does cell mix modify lithium's effect? | joint F, 4 lithium × balance terms | no interaction | F-test |
| **B1** | Does lithium shift whole-transcriptome location? | PERMANOVA pseudo-F, lithium given C | no shift | 9,999 permutations, by margin |
| **B2** | Does lithium change heterogeneity? | PERMDISP F | equal dispersion | 9,999 permutations |
| **B3** | **Is the location shift mediated?** | (i) partial dbRDA, lithium given C + balances; (ii) share of lithium's adjusted R² shared with cell mix; (iii) joint significance of both paths | no remaining direct effect / no mediation | (i) permutation; (ii) bootstrap CI; (iii) both paths at α = 0.05 |
| B4 | Is a dispersion change mediated? *(only if B2 significant)* | NIE, outcome = distance to the non-user centroid | NIE = 0 | as A2 |
| D-B | Within- vs between-cell-type split | lithium's effect on each Shannon component; residual reported | — | regression + bootstrap |
| D-C | Cell-mix diversity | lithium's effect on Shannon of the 5 lineage fractions, and on Aitchison dispersion of the cell mix | — | regression; PERMDISP |
| P | Pathway level | A1–B3 on Hallmark aggregates | as above | as above |
| BD | Bipolar vs control *(exploratory)* | A1, A2, B1, B2, B3 | as above | as above |

### 8.3 Reported whichever way it comes out

Every row above, both trees, every sensitivity, and the residual of the D-B decomposition.
Nothing is dropped for being null.


---
---

## 9. ADDENDUM C — alpha and beta diversity of the cell-type-specific estimates

**Declared 2026-10-27, before any of it was run.** Runs on another machine from
`mediation/CELLTYPE_DIVERSITY_HANDOFF.md`; code `scripts/16_celltype_diversity.R`.

### 9.1 What changes when the outcome is a cell type

In whole blood the question was "does lithium change diversity, and does that run through
the cell mix?". The bMIND estimates are, by construction, expression **within** one lineage,
with the mix already taken out. So:

| | Whole blood (addendum B) | Cell-type estimates (this addendum) |
|:--|:--|:--|
| primary question | is the change mediated by the cell mix? | **does lithium change diversity *inside* a cell type at all?** |
| a lithium effect means | the blood changed, by any route | the cells changed — or the **sub-types inside the lineage** shifted (e.g. eosinophils within granulocytes, neutrophil maturity) |
| mediation through the cell mix | the scientific question | a **leakage check**: it should be ~0. A non-zero indirect effect means the estimates still carry composition |

### 9.2 Decisions

| # | Question | Decision | Why |
|:--|:--|:--|:--|
| 1 | Input | the bMIND posterior means already used by the DE arm (`data/analysis_matrices/bmind/`, exported from `de_analysis/results/bmind_profiles.rds`, verified to reproduce stored values to 2e-5) | same estimates as every other cell-type result; re-running bMIND is ~4 h and not bit-reproducible |
| 2 | Lineages | gran, mono, T **interpretable**; NK and B run but **flagged uninformative** | DE-arm detection floors: 0.73 / 0.74 / 0.77 log2FC vs **4.07 (NK) and 26.11 (B)** |
| 3 | Exposures | lithium within bipolar I (226) **primary**; bipolar-off-lithium vs control (308) **exploratory** | as addendum B |
| 4 | Gene universe | protein-coding ∩ filtered ∩ non-LM22 (10,288), as whole blood | comparability; markers included as a sensitivity |
| 5 | Scale | un-log `2^A − 1`, clipped at 0, renormalised to proportions within each sample and lineage | bMIND worked on log2(TPM+1); Shannon needs linear proportions |
| 6 | Alpha | Shannon H (primary), log Hill q = 2; Hallmark pathway Shannon (fractional counting) | as whole blood. Richness is not computed: bMIND has no detection concept |
| 7 | Beta distance | Euclidean on each sample's centred log2(TPM+1) profile (a log-ratio / Aitchison distance with pseudocount 1 TPM); Jensen–Shannon secondary | as whole blood; the estimates arrive on that log scale |
| 8 | Covariates | the whole-blood set: age, sex, tobacco, RIN, plate, seqPC1–3, log library size; + assessment group for bipolar vs control | identical to the whole-blood analysis |
| 9 | **Shrinkage confound** | every total-effect test is repeated **adjusting for the lineage's own log fraction and its mean posterior SE** | bMIND shrinks a sample's estimate toward the cohort profile more when the lineage is scarce in that sample. Lithium changes fractions, so a dispersion or diversity difference could be manufactured by differential shrinkage. Lithium raises granulocytes: their estimates get *less* shrunk in users |
| 10 | Primary tests | **A1** (lithium → H per lineage) and **B1** (PERMANOVA per lineage), **Holm across the 3 interpretable lineages** within each test | small confirmatory family; NK and B reported but outside it |
| 11 | Secondary | B2 PERMDISP (also on residuals after the shrinkage covariates), A2/B3 mediation through the cell mix as the **leakage check**, pathway layer | |
| 12 | Resampling | alpha: B = 10,000 bootstrap, 10,000 permutations; beta: 9,999 permutations, B = 1,000 bootstrap for the shared fraction | as whole blood |
| 13 | Bridge to whole blood | the per-lineage H already enters the whole-blood decomposition (section 8, D-B); report how A1 for gran/mono/T compares with that decomposition's "within-lineage" component | the two must tell one story |

### 9.3 What would mean what — fixed now

| Result | Reading |
|:--|:--|
| no A1/B1 effect in gran, mono, T | lithium's diversity effect in blood is compositional, consistent with the whole-blood decomposition |
| an effect in a lineage that **survives** the shrinkage covariates | lithium changes that lineage's transcriptome (cell-intrinsic or sub-type shift) |
| an effect that **disappears** with the shrinkage covariates | an artefact of bMIND precision, not biology |
| a non-zero indirect effect through the cell mix | the estimates still carry composition; their per-lineage results must be read with that caveat |
| anything in NK or B | uninformative (precision), whatever its p-value |
