# norm_sweep — varying normalisation alone on the lithium contrast

**Status: DRAFT — permutation null and split-gene control still running. Numbers below
are computed, not estimated; pending sections are marked.**

Everything in this run holds the cohort, the sample subset, the covariate block, the gene
filter and the eBayes settings at the project baseline and changes exactly one thing: how
counts are scaled across samples.

---

## 0. What was held fixed, and why that matters here

| | |
|:--|:--|
| Contrast | `lithium` — bipolar I only, n = 226 (74 non-users, 152 users) |
| Cohort | the 474-sample prep object (`data/prep.rds`), inherited unchanged |
| Covariates | age, sex, tobacco (`tobacco_imp_01`), RIN, plate, seq PC1–3 |
| Gene filter | > 10 counts in ≥ 90 % of samples → **12,173 genes** |
| Estimator | `voom` + `eBayes`, except where limma-trend is named |
| Seed | 20261218 |

The filter runs on **raw** counts, so it is normalisation-independent. Every configuration
below therefore analyses **the identical 12,173 genes** — asserted in code, not assumed.
That is what makes the comparison clean: any difference in output is attributable to the
scaling factors and to nothing else.

Baseline reproduced exactly: 12,173 genes, 1,382 DEG at FDR 0.05, π₀ = 0.5161.

---

## 1. The headline

Among the four legitimate scaling methods the answer barely moves. Between *normalising*
and *not normalising* it moves enormously — and the entire difference is **one number**.

| configuration | DEG (FDR 0.05) | π₀ | Jaccard vs TMM | Spearman *t* vs TMM | global logFC offset |
|:--|--:|--:|--:|--:|--:|
| voom + TMM | 1382 | 0.516 | 1.000 | 1.000 | 0 (reference) |
| voom + TMMwsp | 1348 | 0.531 | 0.917 | 0.9997 | +0.0039 |
| voom + RLE | 1407 | 0.513 | 0.963 | 1.0000 | −0.0011 |
| voom + upperquartile | 1672 | 0.480 | 0.710 | 0.9968 | −0.0169 |
| **voom + none (CPM)** | **4013** | **0.340** | **0.212** | 0.9617 | **−0.0713** |
| limma-trend + TMM | 1306 | 0.525 | 0.854 | 0.9940 | −0.0013 |
| **limma-trend + none (log-CPM)** | **4018** | **0.331** | **0.203** | 0.9581 | **−0.0719** |

Normalisation factor spread (these multiply library size):

| method | min | median | max | max/min | sd(log2) |
|:--|--:|--:|--:|--:|--:|
| TMM | 0.610 | 1.020 | 1.286 | 2.11 | 0.196 |
| TMMwsp | 0.606 | 1.018 | 1.286 | 2.12 | 0.199 |
| RLE | 0.597 | 1.021 | 1.275 | 2.13 | 0.196 |
| upperquartile | 0.600 | 1.023 | 1.250 | 2.09 | 0.183 |
| none | 1 | 1 | 1 | 1 | 0 |

The three "global" methods agree almost perfectly with one another (pairwise Spearman of
the factor vectors 0.986–0.995); upperquartile is the outlier (0.946–0.955), and it is the
outlier in the DE results too.

---

## 2. The confound: normalisation factors track lithium

This was the question the experiment was asked to check, and the answer is yes.

| quantity | non-user → user | unadjusted *p* | adjusted β | adjusted *p* | partial *r* |
|:--|--:|--:|--:|--:|--:|
| log2 TMM factor | −0.0701 | 0.0106 | −0.0718 | **3.09e-4** | −0.244 |
| log2 TMMwsp factor | −0.0723 | 0.0090 | −0.0757 | 2.50e-4 | −0.247 |
| log2 RLE factor | −0.0691 | 0.0112 | −0.0712 | 4.64e-4 | −0.237 |
| log2 upperquartile factor | −0.0557 | 0.0279 | −0.0554 | 1.53e-3 | −0.215 |
| log2 library size | −0.0942 | 0.279 | −0.0259 | 0.265 | −0.076 |
| granulocyte lineage fraction | +0.0484 | 1.46e-4 | +0.0407 | 1.54e-3 | +0.215 |
| ILR b1 (myeloid vs lymphoid) | +0.192 | 0.0110 | +0.156 | 0.0471 | +0.136 |
| ILR b2 (gran vs mono) | +0.126 | 1.38e-3 | +0.103 | 0.0116 | +0.172 |
| CIBERSORTx RMSE | −0.0258 | 7.25e-5 | −0.0242 | 2.47e-4 | −0.248 |

"Adjusted" uses the **same covariate block as the DE model**, so it is the association that
actually bears on the DE result rather than a marginal curiosity.

Note the contrast with **library size**, which is *not* associated with lithium (p = 0.265).
This is specifically a *composition* effect, not a *depth* effect — exactly the distinction
the experiment brief anticipated.

### Why the adjusted association is stronger than the crude one

| covariates | β | SE | *p* |
|:--|--:|--:|--:|
| none | −0.0701 | 0.0274 | 0.0112 |
| + age | −0.0491 | 0.0272 | 0.0725 |
| + plate | −0.0776 | 0.0234 | 0.00106 |
| + seq PC1–3 | −0.0669 | 0.0210 | 0.00167 |
| + full DE block | −0.0718 | 0.0196 | 3.09e-4 |

Adjustment mostly removes **noise**, not confounding: the SE falls 29 % while β barely
moves. Age is the one genuine partial confounder — lithium users are younger (r = −0.203)
and the factor rises with age (r = +0.272), so age alone absorbs about 30 % of the crude
difference. It does not survive the full block.

### How much of it is the cell mixture?

Product-of-coefficients decomposition of lithium → log2(TMM factor), on the same covariate
block. Descriptive, not causal: no sensitivity analysis for unmeasured mediator–outcome
confounding, and the mediator is CIBERSORTx output estimated from the same matrix, so it
inherits the circularity caveat in CLAUDE.md §4.2.

| mediator | a (lithium→M) | b (M→factor) | indirect | direct | % of total |
|:--|--:|--:|--:|--:|--:|
| granulocyte lineage fraction | +0.0407 | −0.994 | −0.0405 | −0.0314 | **56 %** |
| top-50-gene share of library | +0.0137 | −2.08 | −0.0284 | −0.0435 | 39 % |

Lithium users carry **+4.1 percentage points more granulocytes** (0.368 → 0.417 adjusted),
and that shift alone accounts for 56 % of the movement in the scaling factor. Granulocyte
fraction and top-50 share correlate only r = 0.147, so they are largely *separate* routes
into the factor rather than two names for one thing.

A useful external check: lithium-induced benign neutrophilia is well-established
pharmacology, and a +4 pp granulocyte shift is the right sign and a plausible magnitude.
That the CIBERSORTx fractions recover it is reassuring about the deconvolution, quite apart
from what it does to normalisation.

### The mechanism

TMM's factor is set by the dominant transcripts, so anything that concentrates the library
pushes it down:

- top-50-gene share of the library explains **R² = 0.718** of log2(TMM factor)
- granulocyte fraction alone explains **R² = 0.227** (r = −0.476)
- the 4 ILR balances together explain R² = 0.237
- globin/erythroid transcripts explain only R² = 0.020 — **not** a globin artefact
- library size explains R² = 0.015 — **not** a depth artefact

So the chain is: lithium → more granulocytes → library more concentrated in a few
high-abundance transcripts → TMM factor falls. Every link is measured above.

---

## 3. What "no normalisation" actually does

`norm = "none"` nearly triples the DEG count. It does so by adding **one constant** to
every gene's log fold-change:

| quantity | value |
|:--|--:|
| median(logFC_none − logFC_TMM) | −0.07127 |
| mean | −0.07124 |
| **sd** | **0.00337** |
| IQR | 0.00304 |
| lithium coefficient on log2(TMM factor), adjusted | 0.07183 |
| ratio (shift ÷ coefficient) | **0.992** |
| R² of the shift on a constant | **0.9990** |

The shift is 99.2 % of the lithium-associated difference in the scaling factor, and it is
constant to within sd 0.0034 across all 12,173 genes. Dropping normalisation does not
redistribute signal — it adds a location shift.

The genes this manufactures are exactly what a location shift predicts:

| set | n | mean logCPM | % "down in users" | median \|logFC\| under TMM |
|:--|--:|--:|--:|--:|
| DEG under both | 944 | 5.00 | 45.1 | 0.200 |
| **DEG only when un-normalised** | **3069** | 4.72 | **99.97** | **0.072** |
| never DEG | 8160 | 4.57 | — | — |

3,068 of 3,069 point the same way, and their median effect size is 0.072 — the offset
itself. These are one number wearing 3,069 hats.

### The offset explains everything, across all seven configurations

Regressing each configuration's up/down balance on its global offset:

```
log2(up-genes / down-genes) = 1.140 + 55.42 × offset      Pearson r = 0.9986
```

| configuration | offset | up | down | up/down |
|:--|--:|--:|--:|--:|
| voom + TMMwsp | +0.0039 | 995 | 353 | 2.82 |
| voom + TMM | 0 | 956 | 426 | 2.24 |
| voom + RLE | −0.0011 | 943 | 464 | 2.03 |
| trend + TMM | −0.0013 | 889 | 417 | 2.13 |
| voom + upperquartile | −0.0169 | 842 | 830 | **1.01** |
| voom + none | −0.0713 | 519 | 3494 | 0.15 |
| trend + none | −0.0719 | 488 | 3530 | 0.14 |

One scalar per method reproduces the directional composition of every gene list to
r = 0.999. This is the cleanest statement of the result: in this dataset the normalisation
method is, for practical purposes, a single number.

---

## 4. Distributional normalisation (quantile, cyclic loess)

`de_fit()` only exposes edgeR's scaling methods, so these were run by calling voom/limma
directly with the same design.

| configuration | DEG | π₀ | Jaccard vs voom+TMM | offset |
|:--|--:|--:|--:|--:|
| voom + TMM + quantile | 1604 | 0.514 | 0.681 | −0.0008 |
| voom + none + quantile | 1581 | 0.515 | 0.678 | −0.0019 |
| trend + none + quantile | 1538 | 0.529 | 0.648 | −0.0058 |
| trend + none + cyclic loess | 1471 | 0.548 | 0.656 | +0.0016 |

Two things worth noting. First, **quantile normalisation completely rescues the
un-normalised data** (1581 vs 4013), because forcing identical distributions removes the
location shift whether or not a scaling factor was applied — independent confirmation that
the offset is the whole story. Second, quantile disagrees with TMM about as much as
upperquartile does (Jaccard 0.68 vs 0.71) while producing a similar count, so the
disagreement is about *which* borderline genes, not *how many*.

---

## 5. Can normalisation change a conclusion?

Depends which conclusion, and the honest answer has three parts.

**(a) The ranking: no.** Spearman correlation of the *t*-statistic vector against TMM is
≥ 0.9968 for every scaling method and 0.958 even with no normalisation at all. The top gene
is **TSPAN2 (ENSG00000134198)** under all seven configurations, at p = 5.1e-20 to 1.9e-19.
There are **zero sign flips** among genes significant in both of any pair.

**(b) A thresholded gene list: cosmetically, yes; substantively, no** — as long as you
normalise somehow.

| comparison | sig in all | unstable | % of ever-significant that are unstable |
|:--|--:|--:|--:|
| TMM / TMMwsp / RLE | 1294 | 168 | 11.5 % |
| the 4 scaling methods | 1195 | 633 | 34.6 % |
| TMM vs upperquartile | 1268 | 518 | 29.0 % |
| **TMM vs none** | 944 | **3507** | **78.8 %** |

The 34.6 % figure looks alarming until you look at where those genes sit: **629 of the 633
(99 %) have all four adjusted p-values inside [0.01, 0.20]**, median span 0.0385–0.0717.
They churn because they straddle the 0.05 line, not because the methods disagree about
them. That is threshold arbitrariness, not normalisation sensitivity.

**(c) The directional claim: yes, genuinely.** Under voom+TMM the list is 2.24:1
up-regulated; under upperquartile it is 1.01:1. A sentence like *"lithium predominantly
up-regulates blood transcripts"* is true under TMM/TMMwsp/RLE, false under upperquartile,
and inverted under no normalisation. This is the one conclusion in the sweep that a
defensible change of method actually reverses, and it follows directly from §3 — the
up/down ratio is a deterministic function of the residual offset.

---

## 6. The scaling factor carries the exposure signal

Adding log2(TMM factor) to the design as an ordinary covariate collapses the result:

| added covariate | k | DEG | π₀ | % of baseline | Jaccard | extra R² on lithium | SE inflation |
|:--|--:|--:|--:|--:|--:|--:|--:|
| none (baseline) | 0 | 1382 | 0.516 | 100 | 1.000 | — | 1.000 |
| **log2 TMM factor** | 1 | **158** | 0.932 | **11.4** | 0.105 | 0.018 | 1.031 |
| granulocyte fraction | 1 | 234 | 0.822 | 16.9 | 0.135 | 0.014 | 1.024 |
| CIBERSORTx RMSE | 1 | 198 | 0.867 | 14.3 | 0.103 | 0.018 | 1.032 |
| ILR b1 | 1 | 408 | 0.716 | 29.5 | 0.287 | 0.006 | 1.009 |
| all 4 ILR balances | 4 | 350 | 0.763 | 25.3 | 0.210 | 0.012 | 1.020 |
| top-50 count share | 1 | 1159 | 0.534 | 83.9 | 0.637 | 0.006 | 1.011 |
| expression PC1 *(circular)* | 1 | 470 | 0.684 | 34.0 | 0.246 | 0.009 | 1.015 |
| **log2 library size** | 1 | 1439 | 0.514 | 104.1 | 0.850 | 0.002 | 1.003 |
| **log2 TMM factor, permuted** | 1 | 1469 | 0.509 | 106.3 | 0.926 | 0.001 | 1.001 |
| **random N(0,1)** | 1 | 1383 | 0.515 | 100.1 | 0.995 | 0.000 | 1.000 |

The three negative controls (bold, bottom) land within 6 % of baseline, so this is not
"adding any column kills the result". And crucially the **SE inflation is at most 3 %** and
the covariate explains at most 1.8 % additional variance *in lithium itself* — so the
collapse is not collinearity making the exposure unidentifiable. The covariate is eating
variance in the **outcome**, which is what mediation looks like.

### A correction I had to make to my own comparison

The "expression PC1" row above (470 DEG) looked like the factor comfortably beating the
obvious rival explanation — *"it is just the dominant axis of the matrix"*. That comparison
was **unfair**, and the unfairness is instructive. PC1 there was computed on **TMM-normalised**
log-CPM, so the scale axis had already been removed from it before it was asked to compete
for the scale axis. Its correlation with log2(TMM factor) is **0.014**.

Recomputed on **un-normalised** log-CPM:

| added covariate | DEG | π₀ | % of baseline | Jaccard | \|r\| with log2(TMM factor) |
|:--|--:|--:|--:|--:|--:|
| none (baseline) | 1382 | 0.516 | 100 | 1.000 | — |
| log2 TMM factor | 158 | 0.932 | 11.4 | 0.105 | 1.000 |
| **PC1 of un-normalised log-CPM** | **202** | 0.861 | **14.6** | 0.119 | **0.965** |
| PC1 of TMM-normalised log-CPM | 470 | 0.684 | 34.0 | 0.246 | 0.014 |

So the factor does **not** beat PC1 — the factor *is* PC1. The first principal component of
the un-normalised whole-blood transcriptome correlates 0.965 with log2(TMM factor), and
conditioning on either one removes about the same 85–89 % of the lithium gene list. The
right way to say this is not "the normalisation factor is a mediator" but: **there is one
dominant global axis in this matrix, TMM's scaling factor is a near-perfect index of it,
granulocyte content is a large part of it, and the lithium signal lives on it.**

Two caveats are load-bearing and neither is resolved by this table:

1. **Circularity.** The TMM factor, the CIBERSORTx fractions and expression PC1 are all
   computed *from the matrix that supplies the outcome*. This is the threat flagged in
   CLAUDE.md §4.2. Addressed by the split-gene control in §7.
2. **TMM's own assumption.** TMM assumes most genes are not differentially expressed. If
   lithium genuinely shifts the cell mixture, that assumption is **false here by
   construction**, and TMM is not removing an artefact — it is removing the effect of
   interest. Nothing internal to the data decides this. It needs an absolute anchor (a
   spike-in, or a measured cell count), and this dataset has neither.

---

## 7. Is the factor/lithium association fragile?

Everything above rests on one association, so it was stress-tested before being written
down (`norm_sweep_robust.R`; lm only, no model fitting beyond that).

| check | detail | β | *p* | partial *r* |
|:--|:--|--:|--:|--:|
| baseline | full sample, full covariate block | −0.0718 | 3.09e-4 | −0.244 |
| leave-one-out | **worst** of 226 deletions | −0.0680 | 6.75e-4 | −0.231 |
| trimmed 10 % each tail | n = 180 | −0.0413 | 0.0101 | −0.197 |
| drop one plate | **worst** of 5 plates | −0.0581 | 8.01e-3 | −0.203 |
| bootstrap, 2000 resamples | 95 % CI [−0.1080, −0.0343] | −0.0718 | — | −0.244 |

- **Leave-one-out:** p stays between 1.10e-4 and 6.75e-4 across all 226 single-sample
  deletions. No sample is carrying this. The largest single influence on β is 0.005, i.e.
  7 % of the coefficient.
- **Drop-a-plate:** p ranges 4.98e-5 to 8.01e-3 across the five plates. No plate is
  carrying it either.
- **Bootstrap:** 2000 of 2000 resamples give a negative coefficient.
- **Trimming is the one honest caveat.** Removing the 10 % most extreme scaling factors
  from each tail attenuates β by 42 % (−0.0718 → −0.0413), though it stays significant
  (p = 0.0101). So the association is real across the whole distribution but is somewhat
  concentrated in samples with extreme library concentration — which is what the
  mechanism in §2 would predict, since those are the samples whose cell mixture is most
  unusual.

## 8. Split-gene control

*(pending — S1_gene_split_control.csv)*

## 9. Housekeeping normalisation, and the defensible range

*(pending — S2_housekeeping_norm.csv, S3_defensible_range.csv)*

## 10. Permutation calibration

*(pending — permutation_null_summary.csv)*

---

## Files

*(completed at end of run)*
