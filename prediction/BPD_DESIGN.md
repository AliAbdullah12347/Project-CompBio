# Predicting bipolar status from whole blood — design and decisions

Fixed before any AUC was computed. The question is narrow: **does whole-blood
gene expression carry enough signal to separate bipolar I from healthy
controls, beyond what demographics alone already give?**

Every decision below has a reason, and several of them lower the AUC we will
report. That is the point — an AUC that is easy to get here is almost
certainly measuring something other than illness.

---

## 1. Who is a case

**74 bipolar I subjects who are not on lithium.**

Not all 226 BP1. The project's whole premise is that a "bipolar signature" in a
medicated cohort may be a lithium signature, so training on 152 lithium users
would answer a question we have spent the semester arguing is unanswerable from
that contrast. The 74 unmedicated patients are the entire reason this dataset
can address illness at all.

**BP2 excluded.** All 14 are lithium-free and would add cases, but bipolar II
is a different diagnosis and 14 subjects cannot support a separate arm. Adding
them would blur the label for a 19% increase in cases.

## 2. Who is a control — the decision that matters most

**125 controls from assessment group A. Not all 234.**

The deposit splits controls across two assessment groups:

| | BP1 | BP2 | Control |
|:--|--:|--:|--:|
| group A | 226 | 14 | 125 |
| group B | 0 | 0 | **109** |

Group B contains 109 controls and **zero cases**. So in a 74-vs-234 contrast,
47% of controls carry a label that is deducible from their assessment group
rather than their biology.

That would only matter if assessment group were visible in the expression. It
is, overwhelmingly:

```
predicting assessment group from expression, controls only
  PC1 (23.1% var)   AUC 0.787    p = 7e-16
  PC3 ( 6.7% var)   AUC 0.781    p = 4.3e-13
  logistic on top 10 PCs, 5-fold CV:  AUC 0.939
```

Assessment group is one of the strongest signals in the matrix — stronger than
anything we have found for lithium or illness. A classifier trained on the full
control set would find it immediately.

The arithmetic is worth stating plainly. A model that perfectly identifies the
109 group-B controls and is at chance on everyone else scores roughly

```
AUC ≈ (109/234)·1.0 + (125/234)·0.5 ≈ 0.73
```

**An AUC near 0.73 is available here with no biological signal whatsoever.**
Reporting that as evidence that blood predicts bipolar disorder would be wrong,
and it is exactly the result a naive run would produce.

So the primary analysis is **group A only: 74 cases vs 125 controls, n = 199,
prevalence 37.2%**. The full 308-subject version is run as a *secondary*
analysis and reported beside it, because the gap between them is itself the
finding.

## 3. What remains confounded after that

Restricting to group A does not make the contrast clean:

| | cases (74) | controls (125) | p |
|:--|--:|--:|--:|
| age | 53.4 | 46.9 | 0.00076 |
| tobacco | 35.1% | 19.2% | 0.00047 |
| sex (F) | 47.3% | 59.2% | 0.138 |
| RIN, plate, library size, seqPC1–2 | — | — | balanced |

Age and tobacco both have large, replicable effects on the blood
transcriptome. A classifier can score well on either without touching illness.

This is why the covariate-only baseline below is not a formality. It is the
number the gene model has to beat.

**Also note** (CLAUDE.md §2): the deposited QC flag drops 6 controls and 0
BP1 subjects from this contrast, so the control group is filtered in a way the
case group is not. Nothing to do about it; it is stated so it is not
discovered later.

## 4. Feature sets

All four run through the identical CV structure, so their AUCs are comparable.

| name | contents | what it answers |
|:--|:--|:--|
| `covariates` | age, sex, tobacco, RIN, plate, seqPC1–3 | the floor. What demographics alone give. |
| `genes` | filtered genes, TMM + log2 CPM | the question as asked |
| `genes_cov` | both | does expression add to demographics? |
| `genes_resid` | genes with covariates regressed out, fitted on train | is there signal left once confounders are removed from the features? |

`genes` beating `covariates` is the minimum bar. `genes_cov` beating
`covariates` is the honest test of added value. `genes_resid` is the strictest
reading and the one least able to flatter us.

## 5. Cross-validation

Same shape as Figure 5, for the same reason: 25 repeats of stratified 3-fold
CV, every subject out-of-fold exactly once per repeat, 75 train/test pairs.
A single holdout of ~20 cases cannot resolve anything at this sample size.

**Inner loop:** 5-fold CV inside each training fold, for hyperparameters only.
The outer test fold is touched once, to score.

**Everything cross-sample is fitted inside the training fold** — gene filter,
TMM factors and reference, log-CPM scaling, residualisation coefficients,
feature standardisation, hyperparameters. The test fold is transformed with
train-derived values and never contributes to any of them.

**Stratification** is on case/control, which keeps roughly 25 cases in every
test fold. Unstratified 3-fold would occasionally hand a fold 15.

## 6. Models

**Elastic net** (`glmnet`) and **random forest** (`randomForest`), as the
proposal specifies, with the same two applied to every feature set so the
comparison is about the input rather than the algorithm.

- elastic net: α ∈ {0.1, 0.5, 0.9, 1.0}, λ by inner CV
- random forest: 1000 trees, `mtry` ∈ {√p, p/10, p/3} by inner CV,
  class weights inverse to prevalence

## 7. Null model

Labels permuted **within the training fold only**, then the entire pipeline
re-run. 50 permutations. This gives the AUC distribution obtainable from this
cohort, this CV structure and this feature count when no signal exists, which
is a sharper reference than 0.5 — small-n CV with 12,000 features does not
centre at 0.5 in practice.

A result is reported as signal only if it clears the permutation distribution,
not merely 0.5.

## 8. What gets reported

- Pooled out-of-fold AUC per repeat, then the mean across 25 repeats
- 95% CI by bootstrapping **over subjects** (2,000 resamples), not over
  predictions — subjects are the independent unit
- The permutation null distribution beside it
- Every feature set × model combination, including the ones that lose

Pre-declared: we treat whole blood as carrying usable illness signal only if
`genes_cov` beats `covariates` by a margin whose 95% CI excludes zero. If the
transcriptome adds nothing over age, sex and tobacco, that is the result and it
gets reported as such.

## 9. What this cannot answer

With 74 cases, the CI on any AUC will be roughly ±0.07 at best. We can
distinguish "clearly better than demographics" from "no better", but not fine
gradations. A null result here means *this cohort at this size cannot
demonstrate signal*, not that no signal exists.

---

# Covariates — which ones earn their place

Added 14 November 2026, after the instruction to use the covariates that make
sense for lithium rather than one block reused everywhere.

The original code built one covariate matrix — age, sex, tobacco, RIN, plate,
seqPC1–3 — and handed it to both contrasts. For the case-control contrast that
is defensible. For lithium it is not, and it costs us AUC.

## The two criteria

Both are fixed before any AUC is looked at, and neither involves the outcome,
so this is covariate *selection*, not covariate *fishing*.

**(a) A confounder must be causally prior to the exposure.** Lithium was
prescribed years before the blood was drawn. RIN, sequencing plate, library
size and the sequencing-metric PCs are all created by the assay, after the
exposure and after the outcome was biologically determined. They cannot
confound the lithium–expression relationship. Whatever they are, they are not
confounders.

**(b) Among the variables that pass (a), adjustment only helps if the variable
is actually imbalanced.** Threshold |SMD| > 0.1, the usual convention.

## The audit

`18_covariate_audit.R`, output in `results/covariate_audit.txt`:

| variable | non-user | user | SMD | p | prior to exposure? |
|:--|--:|--:|--:|--:|:--|
| **age** | 53.4 | 48.0 | **−0.462** | 0.00072 | yes |
| **sex (M)** | 52.7% | 40.8% | **−0.240** | 0.12 | yes |
| tobacco | 35.1% | 31.6% | −0.075 | 0.70 | yes |
| RIN | 7.49 | 7.48 | −0.017 | 0.91 | no — assay |
| library size | 8.43M | 7.88M | −0.218 | 0.13 | no — assay |
| seqPC1 | 0.007 | −0.000 | −0.156 | 0.27 | no — assay |
| seqPC2 | 0.002 | 0.006 | 0.065 | 0.64 | no — assay |
| seqPC3 | 0.008 | 0.006 | −0.031 | 0.83 | no — assay |
| plate | — | — | — | 0.89 | no — assay |
| assessment group | A | A | — | — | constant, dropped |

**In: age and sex.** Both prior to prescribing, both imbalanced.

**Out: tobacco** — prior, but balanced at |SMD| 0.08. Nothing to correct. It is
a confounder in the case-control contrast (p = 4.7e-04) and that is where the
intuition comes from, but it does not transfer.

**Out: RIN, library size, plate, seqPC1–3** — post-exposure assay variables.

## What the kitchen sink costs

Same CV, same subjects, same model; only the covariate block changes.

```
age alone                     1 var    AUC 0.613
age + sex  [adopted]          2 vars   AUC 0.611
+ tobacco                     3 vars   AUC 0.602
full block [original code]   11 vars   AUC 0.571
```

Eight balanced nuisance variables are worth **−0.042 AUC** against age alone.
They cannot reduce confounding, because there is none to reduce, so all they do
is spend degrees of freedom and dilute the one variable that carries signal.
The reported lithium baseline of 0.568 was understating what demographics
alone give.

This matters beyond the baseline. `genes_resid` residualises expression on this
block, and residualising on eight variables unrelated to the exposure strips
variance without removing bias. One pilot split: `genes_resid` scores 0.887
under age + sex against 0.835 (all splits, full block).

## Covariate sets in the code

`11_run_bpd_split.R --covset`:

| set | contents | used for |
|:--|:--|:--|
| `confounder` | age, sex | **lithium, primary** |
| `confounder_tech` | + RIN, plate, seqPC1–3 | lithium, sensitivity — noise reduction is a separate claim from confounder control, so it is reported separately |
| `full` | + tobacco | **case-control, primary**; also kept so earlier results reproduce |

The case-control contrast keeps `full` on purpose. There the assay variables
genuinely are imbalanced — plate is near-deterministic (chi-square p = 1.6e-09)
because group B was sequenced apart — so criterion (b) catches them even though
criterion (a) says they cannot be causal confounders. A variable that predicts
the label through the study's logistics still has to be controlled.

## One sensitivity run, declared in advance

The full-block lithium results are already reported, so we are not in a clean
pre-specification position and should not pretend to be. The honest handling:
the restricted set is adopted on exposure-side reasoning alone, **both** sets
are reported side by side, and the restricted one is primary because of the
argument above and not because of its score.

---

# What we could add, and what we cannot

## Added: cell composition — as a mediator, not a covariate

The deposited classic CIBERSORT fractions, summed onto a haematopoietic lineage
tree (gran / mono / T / NK / B) and expressed as four ILR balances. Classic
fractions, not B-mode CIBERSORTx: B-mode batch correction is fitted across the
whole cohort at once, which makes it transductive and bars it from Arm 1
(CLAUDE.md 4.1). Zero replacement is cross-sample, so its detection limit is
taken on the training fold only.

**It is not a covariate.** Lithium raises neutrophils — 0.378 to 0.441 here,
p = 2.3e-04, with lymphocytes falling 0.361 to 0.319, p = 0.0038. Composition is
therefore on the causal path from drug to expression. Adjusting away a
confounder removes bias; adjusting away a mediator removes the effect. So it
enters as its own feature set, where the question it answers is asked out loud:

| feature set | question |
|:--|:--|
| `cells` | can a blood count alone identify lithium users? |
| `cells_cov` | cells plus age and sex |
| `genes_cells` | does expression add to the count? |
| `genes_decell` | what survives projecting the cell mixture out of expression? |

Pilot numbers, 25 repeats x 3 folds: `cells` alone **0.615**, `cells_cov`
**0.631** — about what age alone gives, and far below `genes` at 0.906. On one
split `genes_decell` scored 0.893 against `genes` 0.917, which would mean the
lithium signature is mostly *not* cell mixture. One split is not a result; the
full run settles it.

One honest caveat: seqPC1 correlates with neutrophil fraction (Spearman
rho = −0.169, p = 0.011), so the original full block was partly adjusting away
the mediator as well. Small, but it is another reason the block was wrong here.

## Not added, because the deposit does not contain it

Checked directly against the series matrix: all 35 characteristics fields are
already in `metadata_474.csv`. Nothing was left behind. These belong in the
limitations section:

| missing | why it matters for lithium |
|:--|:--|
| **serum lithium level, dose, duration** | no dose-response check is possible; "user" lumps someone started last month with someone at 20 years |
| **concomitant medications** | the largest gap. Lithium users are likelier to be on antipsychotics, valproate or antidepressants, all of which move the blood transcriptome. Any "lithium signature" may be a polypharmacy signature |
| **thyroid function (TSH)** | lithium causes hypothyroidism in 20–30% of users; thyroid hormone has direct, broad transcriptional effects |
| **renal function (creatinine, eGFR)** | lithium is renally cleared and nephrotoxic with long use |
| **measured CBC / white cell differential** | would let us validate the deconvolution against ground truth instead of comparing to published reference intervals |
| **time of blood draw** | neutrophil proportion swings diurnally by more than lithium's effect |
| **BMI, alcohol, fasting state** | each has replicated transcriptome effects |
| **mood state / episode at draw** | euthymic and manic states differ; lithium users are likelier euthymic, so this confounds with treatment |
| **genetic ancestry** | not derivable: we have an HTSeq count matrix, not reads, so no variants can be called |

The first two are the ones to state plainly. Without dose or concomitant
medication, the contrast is "prescribed lithium or not", not "exposed to
lithium at dose X", and it cannot separate lithium from the treatment regimen
that accompanies it.

---

# Results: all three covariate sets, 15 November 2026

`results/summary_cells.csv`. Lithium ran under two covariate sets, the
case-control contrast under `full`.

## The restricted set moves both affected numbers up

| lithium, elastic net | full (11 vars) | confounder (age+sex) | confounder_tech |
|:--|--:|--:|--:|
| `covariates` | 0.568 | **0.589** | 0.571 |
| `genes_resid` | 0.835 | **0.854** | 0.841 |

## The sensitivity arm settles the noise-reduction question: no

Paired by repeat, `confounder_tech` minus `confounder`, elastic net:

```
covariates    -0.018  [-0.028, -0.007]  p = 0.0017
cells_cov     -0.046  [-0.060, -0.032]  p = 0.0000
genes_resid   -0.013  [-0.025, -0.000]  p = 0.0486
genes         -0.005  [-0.011, +0.001]  p = 0.0965
genes_cov     +0.002  [-0.003, +0.007]  p = 0.4552
genes_cells   +0.000  [-0.005, +0.005]  p = 0.9826
genes_decell  +0.001  [-0.008, +0.009]  p = 0.8793
```

Adding RIN, plate and seqPC1-3 never helps. Where it changes anything it
makes things worse, and it is worst exactly where those variables enter the
model or the residualisation. So the assay variables fail on both counts:
they cannot confound a prescription written years earlier, and they do not
buy precision either. `confounder` is the primary set, and this is the
sensitivity analysis that justifies it rather than merely asserting it.

## Arm 1's primary hypothesis is rejected

Do compositional coordinates predict as well as the full transcriptome?

| contrast | cells | genes | difference | margin |
|:--|--:|--:|:--|:--|
| lithium | 0.614 | 0.906 | **-0.292** [-0.303, -0.282] | +/-0.05 |
| bipolar | 0.573 | 0.819 | -0.245 [-0.261, -0.230] | +/-0.05 |

Not equivalent, and not close. A routine blood count does not carry what the
transcriptome carries. Adding cells to genes changes nothing (0.902 vs 0.906),
so the composition contributes no information the genes do not already hold.

Cell composition is not *nothing* for lithium: 0.614 against a measured null
whose 95th percentile is 0.564. It is a real but small signal. For bipolar
status it is indistinguishable from chance (0.573, CI [0.498, 0.644]).

## How much of the lithium signature is cell mixture? Less than expected

Projecting the cell mixture out of the genes costs **-0.046** [-0.055,
-0.038]: 0.906 to 0.860. The drug's expression signature is overwhelmingly
*not* a composition shift, which is the opposite of what the project's
framing anticipated. The neutrophil shift is real (0.378 to 0.441,
p = 2.3e-04) but it is not what the classifier is using.

Worth stating plainly in the paper: we set out expecting cell mixture to
explain much of lithium's apparent transcriptional effect, and it does not.

## Caveats

- Permutation nulls were spent on the feature sets new to each run, so some
  rows have no null of their own. Rows without one are marked `-` in the
  summary; do not quote a null from a different covariate set for
  `covariates` or `genes_resid`, whose features change with the adjustment.
- The random forest still leads elastic net on bipolar `genes_resid`
  (0.572 vs 0.497). The matched-feature run with 1,000 permutations settled
  this: the margin does not survive multiplicity correction, and the same
  forest loses to elastic net on every lithium feature set. Not evidence of
  non-linearity.

---

# The proposal's prediction arm, run as written — 15 November 2026

Earlier runs used four ILR balances on a five-lineage tree. The proposal
specifies **five** balances on six parts (p.7, quoted exactly):

> "five balances follow the hematopoietic lineage: myeloid {neutrophils,
> monocytes} vs lymphoid {CD4 naive, CD4 memory, NK, B}; then neutrophils vs
> monocytes; T {CD4 naive, CD4 memory} vs {NK, B}; CD4 naive vs CD4 memory;
> and NK vs B cells"

and the hypothesis is literally about *five* ratios. The equivalence test
reported before was therefore against the wrong object. `cells5` is now the
primary; the four-balance version is kept as a sensitivity.

Decisions taken where the proposal names cell classes and the deposit names
LM22 types: CD4 memory = resting + activated; NK = resting + activated;
B = naive + memory + plasma cells. CD8 and gamma-delta T are excluded, because
the proposal's "T" is explicitly {CD4 naive, CD4 memory}; the six parts are
renormalised, which is a legitimate sub-composition.

## The five inputs

| # | name in code | status |
|:--|:--|:--|
| 1 | `genes` | run, 75 splits x 2 contrasts |
| 2 | `cells5` | run |
| 3 | `genes_decell5` | run |
| 4 | `ct_neut` | run separately: bMIND, 4,000 genes |
| 5 | `ct_all` | run separately: bMIND, 4,000 genes |

## The pre-declared equivalence test — REJECTED

Equivalent only if the whole 95% CI for (cells5 - genes) lies inside
+/-0.05. A CI containing zero is not equivalence; a CI below -0.05 rejects.

| contrast | model | cells5 - genes | verdict |
|:--|:--|:--|:--|
| lithium | elastic net | **-0.306** [-0.319, -0.292] | rejected |
| lithium | random forest | -0.240 [-0.251, -0.229] | rejected |
| bipolar | elastic net | -0.256 [-0.272, -0.240] | rejected |
| bipolar | random forest | -0.268 [-0.282, -0.254] | rejected |

Not equivalent, in both contrasts, under both models, by a margin five to six
times the one we declared. **Five compositional coordinates are not a
substitute for the transcriptome.** This is the arm's primary result and it
contradicts our stated expectation.

Absolute numbers, lithium, elastic net: `cells5` 0.601 [0.534, 0.663] against
`genes` 0.906 [0.865, 0.940]. Composition is not nothing -- 0.601 clears its
measured null (95th percentile 0.557) -- but it is a small fraction of what
the transcriptome carries.

## Input 3 is the surprise

Projecting the cell mixture out of the genes costs almost nothing:

| contrast | model | decell5 - genes |
|:--|:--|:--|
| lithium | elastic net | -0.051 [-0.064, -0.037] |
| lithium | random forest | -0.016 [-0.026, -0.007] |
| bipolar | elastic net | -0.056 [-0.070, -0.042] |
| bipolar | random forest | **+0.008** [-0.000, +0.016] |

Lithium's signature survives at 0.856 after its cell mixture is removed. We
expected the opposite: the proposal's hypothesis was that "most of the gene
level changes are due to" the cell mix. The neutrophil shift is real
(0.378 to 0.441, p = 2.3e-04) but it is not what the classifier uses.

Both halves of the arm point the same way. Composition cannot do the
transcriptome's job, and removing composition barely dents the transcriptome.

## Five balances vs four

| contrast | cells5 - cells4 |
|:--|:--|
| lithium | -0.016 [-0.025, -0.008] |
| bipolar | -0.012 [-0.023, -0.000] |

The proposal's six-part tree performs slightly *worse* than the five-lineage
one, though both are far from the genes. Dropping CD8 and gamma-delta T costs
a little, which is worth one line in the write-up and changes no conclusion.

## Caveat on this run's nulls

Permutations covered `covariates`, `genes`, `cells5` and `genes_decell5` --
the sets the equivalence test needs. Rows without a null are marked `-` in
`results/summary_proposal.csv`; do not borrow one from another feature set.

## Inputs 4 and 5: cell-type expression (lithium)

bMIND per split, train and test folds deconvolved separately, labels never
passed. 4,000 highest-variance genes selected in the training fold -- against
~12,210 for inputs 1-3. That gap is the cluster's 24-hour job deadline, not a
design choice, and it must be stated wherever these numbers are quoted.

| input | elastic net | vs genes, +/-0.05 margin |
|:--|--:|:--|
| 1 genes (12,210) | 0.906 | - |
| 2 five ILR balances | 0.601 | -0.306 rejected |
| 3 genes, composition out | 0.856 | -0.051 inconclusive |
| 4 neutrophil expression | 0.856 | -0.051 inconclusive |
| 5 all cell types | 0.814 | -0.092 worse |

Input 4 is the finding. Neutrophil-specific expression alone, one lineage and
a third of the gene universe, matches whole blood with composition removed and
comes within the margin's width of the full transcriptome. Krebs identified
neutrophils as the cell type whose PROPORTIONS lithium changes; this says
neutrophils are also the cell type whose EXPRESSION carries the signal, and
that those are two separate facts.

Input 5 is worse than input 4 (0.814 vs 0.856) while containing it plus four
more lineages at five times the features. The other lineages add noise. That
agrees with the DE arm: 109 lithium genes in granulocytes against 28 in T, 20
in monocytes, 4 in NK and 0 in B.

Inputs 3 and 4 are reported as INCONCLUSIVE rather than equivalent. Their CIs
straddle -0.05, so at n = 226 this design cannot resolve whether they fall
inside the margin. Saying "equivalent" would be claiming a precision we do not
have.

**Not run: inputs 4 and 5 for the case-control contrast.** The bMIND worklist
covered lithium only. The bipolar arm therefore stands at three inputs of five.
