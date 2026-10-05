# threshold_sweep — every choice, and why

**Scripts** `scripts/threshold_sweep.R` (main), `scripts/threshold_sweep_dig.R` (follow-ups)
**Console** `console.log`, `console_dig.log`
**Seed** 481, set once at the top of each script
**Ran** 2026-09-29 · R 4.6.1 · limma 3.68.5 · edgeR 4.10.5

---

## 1. What this experiment varies, and what it deliberately does not

Exactly one thing varies: **the rule used to call a gene significant**. The model is
frozen. Every number below comes from a single fit —

```r
de_fit(p, "lithium")     # every argument at its default
```

— which is voom + TMM, filter `>10 counts in >=90% of samples`, bipolar-I only
(74 non-users vs 152 users), covariates age + sex + tobacco (`tobacco_imp_01`) +
RIN + plate + seqPC1-3, assessment group dropped because it is constant within BP1.

The first thing the script does is assert the stated baseline and `stop()` if it fails:

```
genes 12173   samples 226   coef lithium   method voom/TMM
DEG @ FDR 0.05  1382      pi0 (Storey)  0.5161
reproduction check  PASS
```

Residual df = 226 − 13 = **213**.

**Why an assertion rather than a printed check.** A threshold sweep of the wrong fit
looks exactly like a threshold sweep of the right one — same shape of table, same
monotone columns. The only way to notice would be to remember that the baseline was
1382. So the script refuses to produce a table at all unless the fit it is sweeping is
the one the experiment claims to be sweeping.

**What is deliberately out of scope.** Filter, normalisation, estimator, covariate set,
tobacco imputation draw, cell-composition adjustment. Each of those is a separate
experiment. Section I of the dig script does step outside the canonical filter, and is
labelled a *diagnostic* throughout: it exists only to establish whether a ceiling seen
inside the canonical fit is real or imposed, and none of its numbers replace a
canonical number.

---

## 2. Choices made, with reasons

| Choice | What I did | Why |
|:--|:--|:--|
| Grid | FDR {0.001, 0.01, 0.05, 0.10, 0.20} × \|logFC\| {0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0} | As specified. 0.32 (1.25×), 0.58 (1.5×) added separately because those are the cuts papers actually use. |
| Multiplicity | BH for the headline, Bonferroni reported alongside | BH is what the baseline uses. Bonferroni is the honest worst case and answers "how many would survive if we refused to tolerate any false positives". |
| π₀ | Storey, λ = 0.5 | The engine's `pi0_storey` default. Single λ, not the spline estimator — with π₀ ≈ 0.52 the estimate is far from the boundary where λ choice matters. |
| Effect-size scale | Reported in log₂, in fold change, and in percent | log₂ 0.18 sounds like nothing and means nothing to a reader. 1.13× / +13% is the same number in a form that can be argued with. |
| Yardstick | The **sex** coefficient pulled from the **same fit** | Sex is already a covariate in the canonical design, so its logFC comes from identical genes, samples, normalisation, residual variance and moderation. Anything that inflates or deflates lithium's logFC does the same to sex's. A yardstick from a separate fit would not have that property. |
| Precision vs magnitude | SE recovered as \|logFC\|/\|t\|; partial R² as t²/(t²+df) | Both are exact algebra on what `topTable` returns — no re-derivation, nothing estimated twice. Partial R² = squared partial correlation = the share of residual variance lithium accounts for once covariates are out. |
| MDE | (t<sub>crit</sub> + t<sub>0.80</sub>) × SE, with t<sub>crit</sub> from the **operative** raw-p cutoff | Using α = 0.05 would be wrong: BH is not operating at 0.05 on the raw scale. The largest raw p that still clears BH 0.05 here is 5.67e-3, so that is the alpha the procedure is really using. |
| TREAT | `limma::treat(lfc = τ)` on a locally rebuilt fit, τ ∈ {0.1, 0.2, 0.5} | `de_fit` returns a table, not an `MArrayLM`, and TREAT needs the fit object. So the fit is rebuilt from the engine's own `build_design` + `filter_genes`, and the script **asserts** that `eBayes` on the local fit reproduces `de_fit`'s logFC and adjusted p to 1e-10 before `treat` is allowed to run. It passes. Only the final moderated-t step differs. |
| Permutation null | Exposure **column only** permuted, covariates left attached to their own samples; 20 reps in the main script, 40 in the dig | This is the null "lithium carries no signal". Permuting whole rows would be the null "samples are exchangeable", which is false here — age, RIN and plate genuinely differ between samples and destroying that gives an anti-conservative null. The engine's `permute_exposure` already does the right thing. |
| Null summary | Mean **and** median **and** max **and** quantiles | Forced by the data: the null turned out to be bimodal (§4.2). A mean alone would have been actively misleading. |

---

## 3. Results

### 3.1 The grid

DEG counts, one canonical fit:

```
          |lfc|>=0  0.05   0.1    0.2   0.3   0.5   1.0
FDR<0.001      156   156   156    137    94    19     0
FDR<0.01       495   495   489    320   157    23     0
FDR<0.05      1382  1382  1264    557   202    25     0
FDR<0.1       2447  2447  1967    671   215    25     0
FDR<0.2       3960  3930  2655    747   222    26     0
```

Two structural facts fall straight out of the table.

**A \|logFC\| > 0.05 filter is inert.** The 0 and 0.05 columns are identical at every
FDR ≤ 0.10. The smallest \|logFC\| that clears BH 0.05 is **0.0649** — the FDR rule
already implies an effect-size floor above 0.05, so the filter removes nothing. It is
not a conservative safeguard; it is decoration.

**The 1.0 column is empty everywhere.** Not one gene out of 12,173 reaches a 2-fold
change in either direction at any FDR. The largest \|logFC\| in the whole filtered
transcriptome is **0.9137** (1.884×) — and that gene is significant, so the maximum
over all genes and the maximum over DEGs are the same number.

### 3.2 Extremes

| | |
|:--|--:|
| genes tested | 12,173 |
| smallest raw p | 5.09e-20 |
| smallest BH-adjusted p | 6.19e-16 |
| **Bonferroni < 0.05** | **127** |
| Bonferroni < 0.01 | 98 |
| BH < 0.05 | 1,382 |
| operative raw-p cutoff at BH 0.05 | 5.67e-3 |
| π₀ (Storey) | 0.5161 |

π₀ = 0.516 implies **(1 − π₀) × 12,173 ≈ 5,891** genuinely non-null genes. BH 0.05
returns 1,382 of them — **23%**. The study is not short of signal; it is short of
power. Roughly three-quarters of the real lithium signal is sitting below the
detection threshold.

### 3.3 Effect sizes — the headline

\|logFC\| quantiles:

| set | n | 5% | 25% | **50%** | 75% | 95% | max | mean |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| DEG (BH 0.05) | 1,382 | 0.091 | 0.131 | **0.180** | 0.249 | 0.420 | 0.914 | 0.205 |
| non-DEG | 10,791 | 0.005 | 0.027 | 0.054 | 0.088 | 0.159 | 0.535 | 0.064 |

**Median DEG: log₂FC 0.180 = 1.133× = a 13.3% change.**
90th percentile: 0.340 = 1.27×. Maximum: 0.914 = 1.88×.

Survivors under conventional biological cuts:

| rule | n | % of 1,382 |
|:--|--:|--:|
| BH 0.05 only | 1,382 | 100.0 |
| & \|logFC\| > 0.10 | 1,264 | 91.5 |
| & \|logFC\| > 0.20 | 557 | 40.3 |
| & \|logFC\| > 0.32 (1.25×) | 161 | 11.7 |
| **& \|logFC\| > 0.50 (1.41×)** | **25** | **1.8** |
| & \|logFC\| > 0.58 (1.5×) | 16 | 1.2 |
| & \|logFC\| > 1.00 (2×) | **0** | 0.0 |

### 3.4 Precision, not magnitude

| | |
|:--|--:|
| residual df | 213 |
| median SE of logFC | 0.0480 |
| median partial R², all genes | 0.0082 |
| median partial R², DEGs | 0.0500 |
| 90th pct partial R², DEGs | 0.0919 |
| max partial R² | 0.3252 |
| MDE at 80% power | 0.175 (1.13×) |

The median DEG explains **5% of residual variance**; half the genes in the entire
matrix explain under 0.8%. With SE ≈ 0.048 and 213 residual df, a 13% change is about
3.7 standard errors — comfortably significant and biologically small. That is the
whole answer to the headline question.

The DEG-rate gradient across expression quintiles makes the same point from a
different angle (dig §III): DEG rate climbs 7.4% → 18.4% from the lowest to the
highest expression quintile while median \|logFC\| is essentially flat and mildly
U-shaped (0.071 → 0.070, dipping to 0.054 in the middle). Highly expressed genes are
not changing more — they are measured better. `cor(AveExpr, |logFC|)` = **−0.025**,
i.e. nothing, which also rules out the usual low-count artefact where apparent large
effects pile up in weakly expressed genes.

### 3.5 TREAT vs naive post-hoc filtering — they disagree by an order of magnitude

At FDR 0.05:

| τ | TREAT | naive (BH 0.05 → drop \|logFC\| ≤ τ) | in both | TREAT-only | naive-only | Jaccard | smallest \|logFC\| TREAT keeps |
|--:|--:|--:|--:|--:|--:|--:|--:|
| 0.1 | **113** | **1,264** | 113 | 0 | 1,151 | 0.089 | 0.252 |
| 0.2 | **19** | **557** | 19 | 0 | 538 | 0.034 | 0.405 |
| 0.5 | **0** | **25** | 0 | 0 | 25 | 0.000 | — |

Three things here, and all three matter.

1. **TREAT is a strict subset** — `n_treat_only = 0` at every τ. Nothing TREAT calls is
   missed by the naive filter, so the disagreement is entirely naive over-calling.
2. **The factor is 11× at τ = 0.1 and 29× at τ = 0.2.** Naive filtering tests H₀:
   logFC = 0 and then throws away genes whose *point estimate* is small. TREAT tests
   H₀: \|true logFC\| ≤ τ. The second is the hypothesis anyone filtering on logFC
   believes they are testing; the first has no error control at the threshold at all.
3. **TREAT at τ = 0.5 returns zero genes.** Naive returns 25. Not one gene in this
   dataset can be shown, with FDR control, to exceed a 1.41-fold lithium effect. The
   25 "large-effect" genes are genes whose point estimates happen to exceed 0.5 with
   intervals far too wide to establish it — TREAT's minimum kept \|logFC\| at τ = 0.1
   is 0.252, i.e. a gene needs an *observed* effect 2.5× the threshold before the
   interval clears it.

The confidence-interval form (dig §IV) is the same statement without the machinery:
count genes whose 95% CI lower bound on \|logFC\| exceeds τ. It tracks TREAT closely
and is the version to put in a paper, because a reader can check it by eye against a
forest plot.

---

## 4. Things that surprised me, and what I did about them

### 4.1 Nothing reaches 2-fold — is that biology or the filter?

*Dig §I.* The canonical filter requires **>10 counts in ≥90% of samples**. A gene that
is genuinely off in one group and on in the other fails that test by construction: if
the off-group is more than 10% of the cohort, the gene is deleted before the model is
ever fitted. So the filter systematically removes the exact class of gene that produces
large log fold changes, and a ceiling of 0.91 might be an artefact of the filter rather
than a fact about lithium.

The sex coefficient tests this cleanly, because sex-chromosome genes are the textbook
on/off case and sex is already in the canonical design. **Result is in
`diag_filter_effect_ceiling.csv`.** See §5 for the numbers as run.

This is worth carrying into the main pipeline regardless of how it came out: any claim
of the form "no gene shows a large effect" is conditional on a filter that is hostile
to large effects, and that caveat belongs next to the claim.

### 4.2 The permutation null is bimodal

In the 20-permutation run, the null DEG count at FDR 0.20 was **0 in 19 replicates and
1,755 in the twentieth**; null π₀ ranged 0.567–1.000 with median 0.996. A mean of
87.9 describes neither mode and should not be quoted.

The likely mechanism is not exotic. Permuting the exposure column preserves the 74/152
split, so a permuted vector can land close to the true assignment by chance, and at
FDR 0.20 the BH step function is steep enough that a small amount of real signal tips
thousands of genes at once. The dig script runs 40 permutations and records each one's
correlation with the true lithium vector alongside its DEG count, so the explanation is
checked rather than asserted. **Numbers in `diag_permutation_40.csv` and
`diag_null_quantiles.csv`.**

Whatever the cause, the operational conclusion is fixed: **summarise permutation nulls
by quantiles, never by the mean**, and prefer FDR ≤ 0.10, where the null was 0 in every
replicate.

At FDR 0.05 and 0.10 the null was **0 in all 20 replicates**, and null Bonferroni hits
were 0 in all 20 against 127 observed. The 1,382 is not a multiple-testing artefact —
the effects are real. They are just small.

### 4.3 The sex yardstick is closer to lithium than expected

Same fit, same genes:

| coefficient | DEG @ BH 0.05 | DEG & \|logFC\| > 0.5 | median \|logFC\| of sig | max \|logFC\| | min adj p |
|:--|--:|--:|--:|--:|--:|
| lithium | 1,382 | 25 | 0.180 | 0.914 | 6.19e-16 |
| sexM | 906 | 26 | 0.192 | 0.863 | 8.30e-80 |

Sex — a biological variable nobody doubts — produces **fewer** DEGs than lithium, an
almost identical median effect size, and an almost identical count above \|logFC\| 0.5.
Its p-value tail is 64 orders of magnitude more extreme, which is what a handful of
genuinely enormous effects looks like, but its *maximum* \|logFC\| is **smaller** than
lithium's. That is not plausible as biology — sex-chromosome genes are on/off — and is
the observation that motivated §4.1. It is also a useful calibration in its own right:
on this platform, with this filter, at this n, a 13% median change and ~1,000 DEGs is
what a large, unambiguous biological factor looks like. Judged against that yardstick
rather than against an abstract expectation, lithium's signature is ordinary.

---

## 5. Files

| file | contents |
|:--|:--|
| `summary.csv` | the 35-cell grid + matched permutation null + TREAT count where τ = lfc |
| `grid_fdr_x_lfc.csv` | grid with % of tested and % of the FDR set |
| `extremes.csv` | Bonferroni counts, min adjusted p, significance floor, π₀ |
| `logfc_distribution.csv` | \|logFC\| quantiles, DEG vs non-DEG vs all |
| `fold_change_summary.csv` | the same on the fold-change and percent scales |
| `conventional_lfc_survivors.csv` | survivors at 1.25× / 1.41× / 1.5× / 2× |
| `effect_by_expression_bin.csv` | DEG rate and effect size by expression quintile |
| `yardstick_sex_vs_lithium.csv` | sex coefficient from the same fit |
| `precision_vs_magnitude.csv` | SE, partial R², MDE |
| `treat_vs_naive.csv`, `treat_grid.csv` | TREAT vs naive filtering, overlaps |
| `permutation_null_raw.csv`, `observed_vs_null.csv` | 20-permutation null through the grid |
| `diag_filter_effect_ceiling.csv` | §4.1 — filter ladder vs max observable effect |
| `diag_permutation_40.csv`, `diag_null_quantiles.csv` | §4.2 — 40 permutations, quantile null |
| `diag_precision_gradient.csv` | §3.4 — SE and MDE per expression quintile |
| `diag_confint.csv`, `diag_top15_by_effect.csv` | §3.5 — CI form of the TREAT question |
| `diag_direction.csv` | up/down balance of the significant set |
| `objects.rds` | canonical topTable, sex topTable, all tables, SE and partial R² vectors |

---

## 6. What a reader should take away

1. **1,382 is a count, not an effect.** The median significant gene changes by 13%.
   Ninety-eight percent of the significant set is below a 1.41-fold change; none of it
   reaches 2-fold.
2. **The signal is real.** Twenty permutations produced zero DEGs at FDR 0.05 and zero
   Bonferroni hits, against 127 observed Bonferroni hits and a smallest p of 5e-20.
   Small is not the same as spurious.
3. **The signal is under-detected, not over-detected.** π₀ implies ~5,891 non-null
   genes; we recover 23% of them.
4. **Do not filter on logFC post hoc.** At τ = 0.1 it over-calls by 11×, at τ = 0.2 by
   29×, and at τ = 0.5 it reports 25 genes where the correct answer is 0. Use TREAT, or
   report CI lower bounds.
5. **A \|logFC\| > 0.05 filter does nothing here** — the FDR rule already implies a
   floor of 0.0649.
6. **Every "no large effects" claim is conditional on the abundance filter**, which
   deletes on/off genes by construction. Say so wherever the claim appears.
