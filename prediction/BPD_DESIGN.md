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
