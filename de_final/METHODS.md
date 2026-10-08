# Differential expression — methods, explained

The single analysis we report, what each step does, and why each choice was
made. Written to be readable without prior knowledge of RNA-seq statistics.

Run it with `Rscript run_de.R`. Everything below is implemented there.

---

## 0. What we are actually asking

For each of 12,368 genes, one question:

> Is this gene's expression different between two groups of people, once we
> account for the things that differ between those people for reasons other
> than the one we care about?

That is 12,368 separate statistical tests, run four times (two comparisons,
two levels). Almost every complication below exists because of one of two
facts: **the data are counts**, and **there are twelve thousand of them**.

### The four comparisons

| | contrast | n | what it asks |
|:--|:--|--:|:--|
| **WB_LI** | lithium users vs non-users, within bipolar I | 226 (152 vs 74) | what does the drug do? |
| **WB_BPD** | bipolar I off lithium vs healthy controls | 308 (74 vs 234) | what does the illness do, with the drug removed? |
| **CT_LI** | the same lithium contrast, per cell lineage | 226 | which cell type is the drug effect in? |
| **CT_BPD** | the same illness contrast, per cell lineage | 308 | which cell type is the illness effect in? |

WB is whole blood — one measurement per gene per person. CT is cell type — five
measurements per gene per person, one per cell lineage, estimated rather than
observed.

---

## 1. The data, and the two problems with it

The raw matrix holds, for each gene and each person, the number of sequencing
reads assigned to that gene. Bigger number, more of that gene's RNA was
present. Two things stop you comparing those numbers directly.

**Problem 1 — library size.** One person's sample might be sequenced to 10
million reads and another's to 5 million. Every gene in the first looks twice
as expressed, for a reason that has nothing to do with biology.

**Problem 2 — the mean–variance relationship.** With counts, how much a number
bounces around depends on how big it is. A gene averaging 10 reads might vary
between 5 and 20; a gene averaging 10,000 might vary between 9,000 and 11,000.
In absolute terms the second is far more variable, but in relative terms it is
far more stable. Standard linear models assume every observation is equally
precise. Here they plainly are not.

Steps 3 and 4 exist to solve these two problems.

---

## 2. Filtering: which genes to test at all

**What it does.** Keeps genes with more than 10 counts in at least 90% of the
474 samples. 57,773 genes in, **12,368 out**.

**Why bother.** A gene that is switched off in most people carries no
information about group differences, but it still costs you a test — and in
step 6 the price of every test is paid by every other gene. Removing dead
weight makes the surviving genes easier to detect. It is one of the few steps
that both simplifies the analysis and increases power.

**Why this threshold, and not another.** Three reasons, all decided
independently of the answer:

1. **It was fixed in advance**, written into `config.yaml` and committed before
   any analysis was run. A threshold chosen after seeing which one gives the
   most genes is not a test, it is a decision dressed up as one.
2. **It reproduces the source paper.** Krebs et al., who deposited this
   dataset, report 12,344 genes under the same rule; we get 12,368 on our
   474-sample cohort. This is the only external check available on the step.
3. **An independent method agrees.** `edgeR::filterByExpr` is an automatic
   filter we did not choose and did not tune; it keeps 13,925 genes and returns
   a gene count within 2% of ours.

**How much does the choice matter?** We re-ran the entire analysis under **31
different filters**, spanning gene universes from 2,087 to 24,084 genes:

- Across every filter retaining more than ~8,000 genes, the lithium gene count
  ranged from **1,217 to 1,453** — a 19% spread over a threefold change in the
  size of the gene universe.
- **894 of the 1,426** genes we report were called under every filter that
  retained them.
- The median Jaccard overlap of each filter's gene list against ours was
  **0.776**.
- The top gene, TSPAN2, held **rank 1 of 12,368 in all 30 filters that retained
  it**.

The filter is not load-bearing. That is why we report one and not thirty-one.

**One ordering rule that matters:** filter *before* normalising. The scaling
factors in step 3 are computed from the genes you keep, so filtering afterwards
leaves the scaling contaminated by genes you then throw away.

---

## 3. TMM normalisation: making samples comparable

**The naive fix for library size** is to divide each sample by its total count.
That fails when a handful of genes dominate: if one gene is enormous in one
person, dividing by the total shrinks every *other* gene in that person, making
them look under-expressed when nothing happened to them.

**What TMM does instead.** Pick a reference sample. For every other sample,
look at the log-ratio of each gene against the reference. Most genes are not
differentially expressed, so most of those ratios should sit near zero. Trim
away the extremes — the most-changed genes and the most/least expressed — and
take a weighted mean of what's left. That mean is the scaling factor.

The trimming is the point: it measures the shift using the boring majority of
genes rather than letting a few extreme ones set the scale.

**A caveat we report rather than fix.** TMM's scaling factors correlate
**−0.51** with granulocyte proportion across our 474 subjects. Cell composition
therefore leaks into the normalisation itself, meaning the normalisation step
removes a little of the very effect we are studying. We state this; correcting
it is not straightforward and would introduce assumptions of its own.

---

## 4. voom: handling the mean–variance problem

**The problem again.** A linear model assumes every observation is equally
precise. For count data they are not — precision depends on expression level.
Feed counts straight into a linear model and the lowly expressed genes, whose
measurements are least trustworthy, get exactly the same say as the reliable
ones.

**What voom does**, in five steps:

1. Convert counts to log-counts-per-million, which fixes library size using the
   TMM factors.
2. Fit the linear model once, crudely, and record how far each gene's data sits
   from its own fitted line — its residual standard deviation.
3. Plot residual standard deviation against average expression across all
   12,368 genes, and fit a smooth curve through it. This curve *is* the
   mean–variance relationship, measured from this dataset rather than assumed.
4. For every individual observation, read off the curve how variable an
   observation at that expression level should be.
5. Give each observation a **weight** of one over that predicted variance.

Precisely measured observations now count more than noisy ones. The name is
literal: **v**ariance **mo**delling at the **o**bservational level.

**Why voom and not the alternatives.** Two other standard choices exist:
limma-trend and edgeR's quasi-likelihood F test. We report voom because
**it was pre-specified** in `config.yaml` before any analysis ran. That matters
here more than any property of the method, because voom also happens to return
the *most* genes of the three (1,426 against 1,305 and 1,341) — and choosing
the winner afterwards would be exactly the behaviour this project criticises.

**How much does the choice matter?** We ran all three and compared them gene by
gene, not just by headline count:

| | log2FC correlation | p-value rank correlation | genes shared |
|:--|--:|--:|--:|
| voom vs limma-trend | **0.994** | 0.979 | 1,250 (95.8% of the smaller list) |
| voom vs edgeR-QLF | **0.985** | 0.971 | 1,221 (91.1%) |
| limma-trend vs edgeR-QLF | 0.984 | 0.968 | 1,186 (90.9%) |

Of each method's own top 100 genes, **93 are shared by all three**. The
estimator moves the count by about 9% and the gene identities barely at all.

---

## 5. The linear model, and what "adjusting for covariates" means

For each gene, we fit:

```
expression  ~  group + age + sex + tobacco + RIN + plate + seqPC1 + seqPC2 + seqPC3
```

plus `assessment group` for the two case-control comparisons.

**What adjustment actually does.** The coefficient on `group` is not the raw
difference between the two groups. It is the difference *that remains after the
other variables have taken their share*. If lithium users happen to be younger,
and age affects a gene, the raw difference between users and non-users confuses
drug with age. Including age in the model lets it absorb the age-driven part,
leaving a group coefficient that estimates the difference between two people of
the same age, same sex, same smoking status, run on the same plate.

**Why each covariate is there:**

| covariate | why |
|:--|:--|
| **age** | Genuinely imbalanced: within bipolar I, non-users have median age 54.8 vs 48.9 for users (p = 0.0034). Left out, part of the lithium effect would be age. |
| **sex** | Affects expression of many genes; standard. |
| **tobacco** | A strong, well-documented influence on blood expression. Balanced against lithium here (p = 0.70), but not against diagnosis. |
| **RIN** | RNA integrity number. Degraded RNA changes measured expression in ways unrelated to biology. |
| **plate** | Samples processed in five batches; batch effects in RNA-seq are large and routine. |
| **seqPC1–3** | Three principal components of sequencing quality metrics, capturing technical variation not covered by the above. |
| **assessment group** | Controls came from two recruitment groups. Dropped from the lithium contrast because every bipolar I subject is in group A, so it carries no information there and would make the model unsolvable. |

**Verifying the covariates are actually fitted.** Naming a covariate in a
formula is not evidence it was used — this project has already been bitten once
by a function that accepts a covariate argument and silently discards it. So we
check three things, in increasing strength:

1. The design matrix contains the columns. (Necessary, almost worthless.)
2. The design is full rank and `group` is the coefficient being tested. Both
   confirmed for all four comparisons: 226×13 rank 13, and 308×14 rank 14.
3. **Removing the covariates changes the answer.** This is the only check that
   proves they entered the model.

That last test:

| | adjusted vs unadjusted identical? | correlation | genes, adjusted vs unadjusted |
|:--|:--|--:|:--|
| WB_LI | **no** | 0.960 | 1,426 vs 1,500 |
| WB_BPD | **no** | **0.654** | **4 vs 54** |

The bipolar row is a result in its own right. Without covariates the comparison
returns **54 genes**; with them, **4**. Most of the apparent illness signal in
whole blood is age, sex, recruitment group and sequencing batch — not illness.

---

## 6. Empirical Bayes: letting 12,368 genes help each other

**The problem.** Each gene's variance is estimated from 226 or 308 people. That
estimate is itself noisy. By chance, some gene will get an unusually small
variance estimate — and since the test statistic divides by that variance, a
small estimate produces a large t-statistic and a tiny p-value for no good
reason. With 12,368 genes, this is guaranteed to happen many times.

**What limma does.** It assumes the 12,368 gene variances are drawn from a
common distribution, estimates that distribution, and then pulls each gene's
own variance estimate toward it. A gene whose variance looks implausibly small
gets pulled up; one that looks implausibly large gets pulled down. How far each
gene moves depends on how much data supports its own estimate.

The result is a **moderated t-statistic** — the same shape as an ordinary
t-statistic, but with a stabilised denominator. It is both more stable and more
powerful than testing each gene in isolation, and it is the single biggest
reason limma works well on datasets of this size.

---

## 7. Multiple testing: Benjamini–Hochberg at FDR 0.05

**The problem.** Test 12,368 genes at p < 0.05 and, if nothing at all were
going on, you would expect about **618 genes to pass by chance**. Any list
produced that way is mostly noise.

**Two different ways to fix it**, answering two different questions:

- **Family-wise error rate** (Bonferroni, Holm): control the probability of
  making *even one* false claim. Very strict. Appropriate when a single false
  positive is costly.
- **False discovery rate** (Benjamini–Hochberg): control the expected
  *proportion* of false claims among the genes you report.

**BH at 0.05 means:** of the genes on our list, we expect about 5% to be there
by chance. Not "there is a 5% chance of any error" — a weaker, more useful
guarantee.

**How it works.** Sort all 12,368 p-values from smallest to largest. For rank
*k* out of *m* genes, compare p(k) against (k/m) × 0.05. Find the largest *k*
that passes, and call everything up to it significant. The threshold is lenient
for the very smallest p-values and tightens as you go down the list.

**Why FDR and not FWER here.** This is a screen, not a confirmation. The output
is a gene list for follow-up, and the useful property is "most of these are
real", not "none of these could possibly be false". Bonferroni on 12,368 genes
answers a question nobody asked.

**Why 0.05.** Convention, and pre-specified. The number is arbitrary; what is
not arbitrary is choosing it before seeing the results.

**How much does the choice matter?** We computed seven corrections. For the
lithium comparison:

| method | controls | genes |
|:--|:--|--:|
| Holm | family-wise error | 125 |
| Bonferroni | family-wise error | 125 |
| BY | FDR under *any* dependence | 303 |
| **BH** | **FDR** | **1,426** |
| Storey q | FDR, adaptive | 2,568 |

The spread from 125 to 2,568 is not a flaw; it is what choosing an error rate
means. We report BH because it is the field standard for this kind of screen
and was fixed in advance. The other figures are available in `../de_v2/` for
anyone who wants a stricter or looser reading.

---

## 8. The cell-type level

Everything above describes whole blood, where each gene has one measurement per
person. The cell-type level has five, one per lineage, and they are **estimated
rather than measured**.

**Where the cell types come from.** CIBERSORTx estimates what fraction of each
person's blood was each of 22 immune cell types, using the LM22 reference
signature and B-mode batch correction. We aggregate those 22 to **five
lineages** — granulocytes, monocytes, T cells, NK cells, B cells — because at
22-type resolution the agreement between deconvolution methods is poor
(intraclass correlation 0.389) and at lineage level it is acceptable (0.693).
CIBERSORTx is the only deconvolution method used.

**Where the expression estimates come from.** bMIND takes the filtered 12,368
genes as log2(TPM+1), together with the five lineage fractions, and estimates
what each lineage was expressing in each person. It is run **without the
diagnosis labels**, so the deconvolution cannot be shaped by the comparison it
will later be used for.

**Two changes to the model:**

1. **Precision weights.** bMIND reports its own uncertainty — a standard error
   for every gene × lineage × person. We weight each observation by 1/SE², so
   confident estimates count more. Treating the estimates as if they were
   measurements would understate uncertainty and inflate the gene count.
2. **BH within each lineage.** Each lineage is its own family of 12,368 tests,
   which answers "is there signal in granulocytes?". We also report a **pooled**
   BH over all 61,840 gene × lineage tests, because that is the only fair basis
   for comparing the cell-type level against whole blood — otherwise the
   cell-type side gets five independent chances at significance and wins on
   multiplicity alone.

---

## 9. Results

| comparison | n | tests | significant at BH 0.05 | smallest BH q | largest \|log2FC\| |
|:--|--:|--:|--:|--:|--:|
| **WB_LI** | 226 | 12,368 | **1,426** | 7.0×10⁻¹⁶ | 0.911 |
| **WB_BPD** | 308 | 12,368 | **4** | 0.0489 | 0.917 |
| **CT_LI** | 226 | 61,840 | **161** (66 unique genes pooled) | 3.0×10⁻⁸ | 0.829 |
| **CT_BPD** | 308 | 61,840 | **0** | 0.091 | 0.714 |

Cell-type results by lineage:

| | granulocytes | monocytes | T cells | NK cells | B cells |
|:--|--:|--:|--:|--:|--:|
| **CT_LI** | **109** | 20 | 28 | 4 | 0 |
| **CT_BPD** | 0 | 0 | 0 | 0 | 0 |

### Reading these

**WB_LI is the large, robust result.** 1,426 genes, smallest p 5.7×10⁻²⁰, and
it survives every robustness check we ran — three estimators, 31 filters, a
quality screen, and multiple-imputation pooling of the missing covariate.

**WB_BPD's four genes should not be reported as findings.** All four sit on the
boundary at q = 0.0489. None survives a stricter correction, none clears a
minimum fold-change threshold, none is called under all 31 filters, and pooling
the 20 tobacco imputations correctly returns zero.

**CT_LI localises the effect to granulocytes** — 109 of the 161 significant
gene × lineage pairs.

**CT_BPD's zero is not evidence of absence**, and this is the most important
caveat in the analysis. We measured how small an effect each lineage could
detect by planting artificial effects of known size and counting how often they
were recovered:

| lineage | share of blood | smallest detectable effect |
|:--|--:|--:|
| granulocytes | 38.7% | 0.73–0.83 log2FC |
| monocytes | 27.2% | 0.74 |
| T cells | 24.7% | 0.77 |
| NK cells | 7.9% | 4.07 |
| **B cells** | **1.45%** | **26.11** (a 72-million-fold change) |

Every one exceeds the 0.5 log2FC floor we set in advance as the smallest
biologically meaningful change. Whole blood, by contrast, detects **0.203
log2FC** — about a 15% change — so the whole-blood null *is* interpretable and
the cell-type nulls are not.

---

## 10. What was deliberately left out

Earlier work explored variants that are not reported here. They are kept in
`../de_v2/` as the evidence for the choices above, not as parallel results.

| left out | why |
|:--|:--|
| limma-trend, edgeR-QLF | Agree with voom at r > 0.98; reporting three near-identical lists adds length, not information. |
| 30 other gene filters | The filter is not load-bearing (§2). |
| Holm, Bonferroni, Hommel, BY, Storey, permutation FDR | One error rate, chosen in advance (§7). |
| Composition-adjusted fits | See below. |
| Multi-method deconvolution comparison | CIBERSORTx only. |

**On the composition-adjusted fits.** These were dropped as a *parallel set of
gene lists*, which they should be — they are the same hypotheses asked a second
way and they double every table. But the per-gene *comparison* between adjusted
and unadjusted is not a variant; it is the answer to the question this project
exists to ask, and it is reported separately: 20.5% of the 1,426 lithium genes
survive adjustment for cell composition, rising monotonically with signal
strength to **100% of the 125 that pass family-wise error control**. Composition
explains the breadth of the lithium signature but not its core.

---

## 11. Limitations

- **"Off lithium" is not "unmedicated."** The deposit records lithium use and
  nothing else — no antipsychotics, anticonvulsants, antidepressants, dose or
  duration. Every WB_BPD and CT_BPD claim inherits this.
- **No cell-type null is interpretable** (§9).
- **Cell fractions are estimated from the same expression matrix** they are used
  to explain. For the 547 LM22 marker genes this is close to circular.
- **LM22 misfits whole blood in the rarer compartments.** Our estimates put
  monocytes at 26.6% against a clinical reference range of 2–8%, and eosinophils
  at 0.115% — exactly zero in 410 of 474 samples — against 1–6%. A composition
  adjustment built on LM22 cannot remove variance from a compartment the
  signature does not resolve.
- **Five lineages, not 22 cell types**, so a shift *within* a lineage is
  indistinguishable from a change inside its cells.
