# method_sweep — vary the estimator, hold everything else fixed

**Status: methods and rationale final; results section filled in on completion.**

## 1. What this experiment varies

One axis only: the *estimator* used to test the lithium coefficient. Nothing
else moves. Every configuration sees the same 226 samples, the same design
matrix, the same gene filter and the same TMM normalisation, so any difference
in the answer is a property of the estimator and not of the pipeline around it.

Five estimators, each run twice (`robust = FALSE` and `robust = TRUE`):

| id | what it is |
|:--|:--|
| `voom` | `limma::voom` → `lmFit` → `eBayes`. The project baseline. |
| `voomWQW` | `limma::voomWithQualityWeights`. voom plus per-sample quality weights. |
| `trend` | limma-trend: `cpm(log = TRUE, prior.count = 3)` → `lmFit` → `eBayes(trend = TRUE)`. |
| `QLF` | edgeR quasi-likelihood: `estimateDisp` → `glmQLFit` → `glmQLFTest`. |
| `LRT` | edgeR likelihood ratio: `estimateDisp` → `glmFit` → `glmLRT`. |

`robust` does not mean the same operation in all five, and it is worth being
explicit because it is easy to read the column as if it did:

- `voom`, `voomWQW`, `trend` — `robust = TRUE` is passed to `eBayes`, which
  makes the variance-shrinkage prior resistant to individual hypervariable
  genes (gene-specific `df.prior`).
- `QLF` — passed to both `estimateDisp` and `glmQLFit`, so it affects the
  dispersion trend *and* the QL dispersion shrinkage.
- `LRT` — `glmFit`/`glmLRT` take no `robust` argument, so the flag reaches only
  `estimateDisp`. `LRT_rob` therefore differs from `LRT` solely through the
  dispersion estimates. This is a real difference, not a no-op, but it is a
  smaller intervention than in the other four.

## 2. Fixed configuration

Taken from `common.R`'s `CONTRASTS$lithium`, unchanged:

- **Samples**: bipolar I only, n = 226 (74 lithium non-users, 152 users).
- **Design**: `~ lithium + age + sex + tob + rin + plate + seqpc1 + seqpc2 + seqpc3`
  → 13 columns, full rank 13. Assessment group is constant within BP1 and is
  dropped by `build_design`; `dropped_constant` and `dropped_aliased` were both
  empty, confirming nothing else was silently removed.
- **Tobacco**: completion column `tobacco_imp_01` (the engine default).
- **Filter**: > 10 counts in ≥ 90% of the 226 → 12,173 genes.
- **Normalisation**: TMM.
- **Coefficient tested**: `lithium`.

The `voom` configuration reproduces the stated project baseline exactly —
12,173 genes, 1,382 DEG at FDR 0.05, π₀ = 0.5161 — which is the check that the
sweep is varying what it claims to vary.

## 3. Metric definitions

Several of these are easy to define two ways, so the exact choice is recorded.

- **DEG**: `adj.P.Val < 0.05`, Benjamini–Hochberg, as produced by `topTable`
  (limma) or `topTags` (edgeR). No log-fold-change threshold unless the column
  name says otherwise.
- **π₀**: Storey's proportion of true nulls at λ = 0.5, from `common.R`. Also
  reported at λ = 0.8, because a p-value histogram with a non-flat right tail
  makes the λ = 0.5 estimate method-dependent in its own right.
- **λ_GC**: `median(qchisq(1 - p, 1)) / qchisq(0.5, 1)`, the genomic-inflation
  factor. Under genuine signal this is legitimately > 1, so it says nothing in
  isolation. It is used here only *across* methods on identical data, where the
  differences are attributable to the estimator.
- **Signed test statistic**: limma methods supply `t` directly. `F` (QLF) and
  `LR` (LRT) are non-negative by construction, so the signed statistic is
  `sign(logFC) · sqrt(F)` and `sign(logFC) · sqrt(LR)`. `sqrt` puts them on a
  t-like scale and is monotone, so every Spearman correlation involving them is
  unaffected by the transform; the sign is what actually matters, because
  without it two methods that disagree about direction would look concordant.
- **Jaccard**: |A ∩ B| / |A ∪ B| over DEG sets at FDR 0.05.
- **Containment** (asymmetric): fraction of the *row* method's DEGs that the
  *column* method also calls. Reported alongside Jaccard because Jaccard
  confounds "disagrees with" and "is a smaller set than".
- **Top-K overlap**: fraction shared between the two methods' K smallest
  p-values, K ∈ {100, 500, 1000, 2000}. This measures agreement about ranking,
  which is independent of where each method happens to put its threshold.
- **Kish effective sample size**: (Σw)² / Σw², with w the voomWQW sample
  weights. limma scales array weights to mean 1, so the spread of w *is* the
  information lost, and n − n_eff expresses it in samples.

## 4. Choices made, and why

**Common gene universe, asserted not assumed.** The filter is a cross-sample
operation, so a different sample set yields a different gene list. Here the
sample set is identical in all ten configurations, so the gene list must be
identical too — the script `stopifnot`s this before computing any Jaccard. If
it ever failed, every concordance number would be silently comparing different
universes.

**CPU seconds, not wall clock, as the runtime metric.** Five to ten sibling
experiments were running on this machine throughout. Wall clock therefore
measures contention, not the estimator. Both are recorded; `cpu_seconds` is the
one to read, and even it is only accurate to within the memory-bandwidth
contention that CPU time cannot subtract. Treat the timings as order-of-
magnitude.

**Permutation null of the exposure, not of the samples.** The `permute_exposure`
path in `de_fit` permutes the lithium column of the design and leaves every
covariate attached to its own sample. That is the null "lithium carries no
signal", which is the hypothesis being tested, rather than "the samples are
exchangeable", which is false here (age differs by lithium status, p = 0.0034).
All methods are given the *same* permutations, so the comparison is paired.

**Why a permutation null at all.** A DEG count cannot by itself say which method
to believe: a method can lead the table because it is more powerful or because
it is anti-conservative, and those look identical from the count alone. Under a
correct null, 5% of p-values should fall below 0.05 and essentially no gene
should survive FDR. The ratio of the observed fraction to 0.05 is reported as
`excess_type1_at_05`.

**Why split-half as well.** The permutation null detects miscalibration but not
over-fitting of the variance model. A method that calls more genes because it
is genuinely more powerful will also agree with itself better across two halves
of the cohort; one that calls more genes because it is unstable will not. The
split is stratified on lithium so both halves carry the same exposure ratio —
an unstratified split can hand one half far fewer non-users, and the
reproducibility drop would then be a sample-size artefact rather than a
property of the estimator. The gene filter re-runs inside each half (it is a
cross-sample operation and freezing it would leak), so the halves are compared
on their intersected gene lists.

**Decomposing voomWQW.** `voomWithQualityWeights` does two things at once: it
estimates sample weights, and it re-runs the voom mean–variance fit with those
weights in place. To separate them, plain voom's observation weights were
multiplied column-wise by voomWQW's sample weights and refitted. The three-way
comparison (voom / weights-only / full voomWQW) says which of the two
operations is responsible for the difference. A small synthetic check first
confirmed that `voomWQW$weights / voom$weights` recovers the sample weights
column-wise, i.e. that the multiplication is the right decomposition.

**Downweighting versus excluding.** The project deliberately retains samples
that its own depth and RIN screen flags. voomWQW is the only estimator in the
sweep that can react to them. So the sweep also asks whether downweighting is
doing the job an exclusion would have done: plain voom on the 217 clean samples
is compared with voomWQW on all 226. A second exclusion — the nine
*lowest-weighted* samples, whatever they are — is run alongside, so that
"exclude what voomWQW distrusts" and "exclude what our thresholds distrust" can
be told apart.

**A permutation null is conservative when real signal exists — read it
comparatively, not absolutely.** Permuting the lithium column moves the
variance explained by a genuine lithium effect into the residual, which
inflates the residual variance of exactly the genes that carry signal. Because
`eBayes` and edgeR's dispersion estimators share information across genes, that
inflation is partly passed on to *every* gene, and the whole null p-value
distribution shifts right. With ~11% of the tested genes called at FDR 0.05 in
the observed fit, this is not a negligible effect. So "0 DEG under the null" is
necessary but not sufficient evidence of calibration, and a λ_GC below 1 under
permutation does not on its own prove a method is conservative. What the null
*can* do validly is rank the five estimators against each other, because all
five see the identical permutations and identical data.

**An unplanned reproducibility check, worth recording.** The `norm_sweep`
experiment in this folder, written independently, also runs a permutation null
on `voom` + TMM with `set.seed(20261218)`. Its first ten permutations and this
experiment's first ten are identical to four decimal places on every reported
quantity — p<0.05 fraction 0.0283, 0.0171, 0.0139, 0.0357, 0.0106, 0.0271,
0.0372, 0.0654, 0.0299, 0.1087, and π₀ 1.000, 1.000, 1.000, 1.000, 1.000,
1.000, 0.979, 0.955, 1.000, 0.775 in both. Two scripts that share only
`common.R` and a seed therefore draw the same permutations and land on the same
numbers. That is the seed-threading rule in section 4.5 of the project brief
working as intended, demonstrated rather than asserted.

**Why `voomWQW` is so much more expensive.** `voomWithQualityWeights` defaults
to `method = "genebygene"` with `maxiter = 50`, an incremental REML update that
loops over all 12,173 genes in R, and it runs the weight estimation twice
(once on the initial voom fit, once after re-voom). The cost is a property of
that default, not of the weighting idea — `method = "reml"` would be far
cheaper. It was left at the default because the experiment's premise is to
measure what a reasonable user gets from the package as shipped.

**Biotype composition.** A method that is picking up low-count noise shows it
as an excess of pseudogene and non-coding calls relative to the tested
background. `data/gene_annot.csv` (written by a sibling experiment in this
folder) supplies the annotation; it is read defensively and the stage is
skipped rather than guessed at if the file is absent or has changed shape.

## 5. What was deliberately *not* done

- **DESeq2 was not added.** It is not installed in this R library, and adding a
  dependency is not this experiment's call to make. The sweep therefore covers
  the limma/edgeR family only. That family shares a normalisation and a
  shrinkage philosophy, so the concordance numbers below are an *optimistic*
  bound on method dependence — a genuinely independent estimator would very
  likely sit further out than anything here.
- **No change to `common.R`.** Every fit goes through the shared `de_fit`. The
  three places that reach past it (extracting voomWQW's sample weights, the
  weights-only decomposition, and the edgeR dispersion parameters) rebuild
  their inputs with `build_design` and `filter_genes` from `common.R`, so they
  are the same objects `de_fit` uses and not a lookalike assembled a second
  way.
- **No tuning.** Every estimator runs at its package defaults. The point is to
  measure the spread a reader would get by making a reasonable-looking choice,
  not to find the configuration that gives the largest number.

## 6. Files

| file | contents |
|:--|:--|
| `summary.csv` | master table, one row per configuration |
| `per_method.csv` | full per-configuration metrics |
| `concordance_*.csv` | Jaccard, containment, Spearman, top-K matrices |
| `oddness.csv` | mean distance of each configuration from the rest |
| `voomWQW_sample_weights.csv` | all 226 weights with depth, RIN, exposure |
| `weight_correlates.csv`, `weight_group_tests.csv`, `weight_joint_model.txt` | what predicts the weight |
| `effective_sample_size.csv` | Kish n_eff overall and by group |
| `voomWQW_decomposition.csv` | weights vs re-fitted mean–variance trend |
| `downweight_vs_exclude.csv` | voomWQW vs hard exclusion of the flagged samples |
| `shrinkage_parameters.csv` | df.prior / s2.prior / dispersion per configuration |
| `deg_by_abundance_decile.csv`, `discordant_gene_profile.csv` | where on the abundance axis the methods disagree |
| `deg_biotype_composition.csv`, `most_disputed_genes.csv` | what kind of gene each method calls |
| `permutation_null_*.csv` | type-I calibration |
| `splithalf_*.csv` | reproducibility |
| `fig_*.png` | weight vs quality, p-value histograms, concordance heatmap |
| `console.log` | full run log |

## 7. Results

*(filled in on completion)*
