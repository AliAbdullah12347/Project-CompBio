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
