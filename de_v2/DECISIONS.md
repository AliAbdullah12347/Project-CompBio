# de_v2 — every decision, and why

The differential expression analysis re-run without verdict bands, with a full
multiple-testing panel. Written as the work was done. Where this version
departs from `../de_analysis/` (v1), both the old choice and the reason for
changing it are kept rather than overwritten.

Companion documents: `README.md` (what is here and how to run it),
`../de_analysis/METHODS.md` (v1's decision log, still the reference for
everything v2 reuses rather than recomputes).

---

## 1. The bands are gone

v1 converted every result into SUPPORTED / INCONCLUSIVE / REFUTED by comparing
the percentage of significant genes against two thresholds, 0.5% and 2%.

Those bands had a real purpose: committed in `config.yaml` before any analysis
ran, they stopped a threshold being chosen after seeing the answer. That is
worth something, and `config.yaml` is kept as the timestamped record of it.

They are dropped here for three reasons.

**They destroy information.** "0.881% of genes in granulocytes" carries strictly
more than "INCONCLUSIVE", which is what the band turned it into. A reader who
wants to judge the size of an effect cannot recover it from a verdict.

**They are an uncorrected second test.** The band is itself a decision rule
applied to data, layered on top of the gene-level testing, and nothing controls
its error rate. There is no sense in which a band has a false-positive
probability — it is a cutoff with no sampling theory attached.

**The verdict wording misleads.** v1's CT_BPD came out "REFUTED", which reads as
*shown absent*. What actually happened is that a prediction was not borne out by
a test that, as `../de_analysis/METHODS.md` §17 shows, had almost no power to
bear it out. The asterisk needed to repair that is a sign the vocabulary was
wrong.

**What replaces them.** Quantities and their uncertainty. Discovery counts are
reported at several thresholds as *columns in a table* — FDR 0.01, 0.05, 0.10,
FWER 0.05 — so the reader sees how the answer moves with the cutoff instead of
being handed one. Nothing in `scripts/` returns a verdict, and `mtc.R`
deliberately has no opinion about 0.05.

---

## 2. What is re-run, and what is reused

**Re-run from scratch:** every model fit. All four comparisons, both
adjustments, eight fits in total, from the count matrix upward. No result file
from v1 is read.

**Reused as input, not as results:**

| file | why it is not recomputed |
|:--|:--|
| `de_analysis/data/prep.rds` | The gene filter, TMM normalisation, TPM, union-of-exons gene lengths, lineage aggregation and ILR balances. None of it depends on a threshold or on the group labels. |
| `de_analysis/results/bmind_profiles.rds` | 3.1 hours of deconvolution. bMIND was run **without** the phenotype, so the cell-type estimates do not depend on the labels at all; re-running it would reproduce the same numbers. |

Reusing the deconvolution is not a shortcut, and it is worth being explicit
about why, because it also licenses the permutation null in §9: since the
labels never entered the deconvolution, permuting them afterwards is a complete
null rather than an approximation.

**Deliberately not revisited:** the per-lineage detection floors from
`de_analysis/scripts/08` and `08b`. They are properties of the design and the
estimator, not of the correction method, so the correction panel here does not
change them. They remain the authority on what a cell-type null can mean, and
§11 restates the consequence.

---

## 3. The four comparisons

| | contrast | n | tests |
|:--|:--|--:|--:|
| **WB_LI** | whole blood; lithium users vs non-users, within bipolar I | 226 (152 vs 74) | 12,368 |
| **WB_BPD** | whole blood; bipolar I off lithium vs healthy controls | 308 (74 vs 234) | 12,368 |
| **CT_LI** | per lineage; same contrast as WB_LI | 226 | 61,840 |
| **CT_BPD** | per lineage; same contrast as WB_BPD | 308 | 61,840 |

Each is fitted twice — unadjusted, and adjusted for the four ILR
cell-composition balances. For whole blood the adjusted fit asks how much of
the effect survives removing the cell mixture. For cell type it is the
circularity control: bMIND's posterior for a lineage depends on that sample's
fraction for that lineage, and both contrasts move those fractions.

**These are not four independent studies.** WB_LI and CT_LI are the same 226
subjects measured two ways; WB_BPD and CT_BPD the same 308. The cell-type and
whole-blood versions of one contrast are near-duplicates. Any reasoning about
error rates across the four has to treat them as strongly positively dependent,
which is what §8 does.

---

## 4. Design and covariates — unchanged from v1, and checked again

Same covariates as v1, identical across levels so that whole blood and cell
type differ only in what is being measured:

- WB_LI / CT_LI: age, sex, tobacco, RIN, plate, seqPC1–3
- WB_BPD / CT_BPD: the above plus assessment group

`assessment group` is dropped from the lithium contrasts because every bipolar I
subject is in group A, making it inestimable there; it varies once controls
enter. Every design matrix is asserted full rank at fit time
(`stopifnot(qr(X)$rank == ncol(X))`), so a rank-deficient design fails loudly
instead of silently dropping a column.

Cell-type fits use **1/SE² precision weights** from bMIND's posterior standard
errors. Treating posterior means as measured data would understate uncertainty
and inflate the cell-type discovery count — biasing exactly the comparison this
project cares about.

---

## 5. Which corrections, and why each one is here

Every method is applied to every analysis, so no result depends on which
correction a reader happens to prefer.

| method | controls | valid under | why it is in the panel |
|:--|:--|:--|:--|
| **BH** | FDR | independence or PRDS | The field standard. The primary number. |
| **Storey q** | FDR | same, plus estimated π₀ | Adaptive BH. Uniformly at least as powerful, and the gain is large when π₀ is well below 1 — which it is here. |
| **BY** | FDR | **arbitrary** dependence | The honest bound if one refuses to assume PRDS for co-expressed genes. Costs a factor of Σ(1/i) ≈ 9.9 at m = 12,368. |
| **Holm** | FWER | arbitrary dependence | For the "which genes would I actually bet on" list. |
| **Hommel** | FWER | PRDS | More powerful than Holm in principle. |
| **Bonferroni** | FWER | arbitrary dependence | Never preferable to Holm, which dominates it. Reported only because reviewers ask for it by name. |
| **Permutation FDR** | FDR | **nothing** | Measures the null instead of assuming it. See §9. |

**Why both π₀ estimators.** Storey's q-value needs the share of true nulls.
Two estimators are reported — the spline extrapolation of Storey & Tibshirani
(2003) and the bootstrap-λ choice of Storey, Taylor & Siegmund (2004) — because
they disagree when the p-value histogram is irregular, and quoting one alone
would hide that. Here they agree closely (within 0.02 everywhere), which is
itself reassuring.

**Hommel is skipped above 20,000 tests.** Its algorithm is quadratic in the
number of tests; on the 61,840-test cell-type families it would run for hours to
gain nothing measurable over Holm. It is computed for the whole-blood families,
where it turns out to equal Holm exactly.

**Independent filtering is reported, not applied.** Bourgon et al. (2010): drop
the genes least likely to be detectable using a statistic independent of the
p-value under the null — mean expression qualifies, since it ignores group — and
power rises for whatever survives. It is reported as a diagnostic rather than
folded into the headline, because our gene filter already removed low-expression
genes. The result confirms that: the optimal additional filter is the 0%
quantile for six of eight analyses, i.e. **further filtering buys nothing**, and
where it helps at all it moves WB_BPD from 4 genes to 5.

---

## 6. Why these are hand-written, and how they are checked

This machine has limma, edgeR, statmod and MIND. It does **not** have `qvalue`,
`IHW`, `fdrtool`, `ashr`, `locfdr` or `swfdr`, and CLAUDE.md §7 says to ask
before adding a dependency. So everything beyond base R's `p.adjust` is
implemented in `scripts/mtc.R` from its published definition.

That is only acceptable if the implementations are **checked rather than
trusted**. A hand-written FDR method that is subtly wrong does not announce
itself; it returns a plausible number of genes. `scripts/00_validate_mtc.R`
pins each one to an identity that must hold exactly, or to a simulation with
known ground truth, and the pipeline will not run if any check fails.

All 18 checks pass. The ones that matter most:

| check | result |
|:--|:--|
| Storey q at π₀ = 1 equals `p.adjust(.,"BH")` **exactly** | max difference 2.2×10⁻¹⁶ |
| π₀ estimators on a pure null (true value 1) | spline 0.9974, bootstrap 0.9939 |
| π₀ estimators on a known 70% null mixture | spline 0.6944, bootstrap 0.7005 |
| Benjamini–Bogomolov with K = 1 reduces to BH | max difference 0 |
| Permutation FDR against a uniform null matches BH | max difference 0.0055 |
| Realised false discovery proportion, 150 simulated datasets, target 0.05 | BH 0.0460, Storey 0.0487 |
| Storey at least as powerful as BH | 0.5020 vs 0.4939 |
| π₀ estimate vs truth (0.95) | 0.9485 |

**Scope of the simulation.** Only hand-written methods need operating-
characteristic validation. BH, BY, Holm, Bonferroni and Hommel come from base R
and are not ours to validate — Part 1 of the script only confirms this file
wires them together in the right order (Bonferroni ≥ Holm ≥ Hommel, Holm ≥ BH,
BY ≥ BH). BH is carried through the simulation purely as the reference Storey
must not beat on error rate.

---

## 7. Two layers of multiplicity that a per-contrast correction misses

Correcting within each analysis is the usual practice and it is not enough here.

**Across lineages.** A cell-type contrast tests every gene five times. Reporting
"significant in any lineage" after correcting inside each lineage gives the
cell-type level five independent chances, so it would beat whole blood on
multiplicity alone. Three treatments are reported, because they answer three
different questions:

- **per-lineage BH** — "is there signal in granulocytes?"
- **pooled BH** over all gene × lineage tests, collapsed to unique genes — "is
  there signal anywhere?", and the only fair basis for comparing the cell-type
  level against whole blood. The percentage uses unique genes out of 12,368, not
  61,840 tests, so the two levels measure the same thing.
- **Benjamini–Bogomolov** — the per-lineage question, but paying for having
  looked at all five and written up the interesting ones.

**Across contrasts.** Four comparisons were run and the interesting ones get
written up. Nothing in a per-contrast correction accounts for that selection.
§8 handles it.

---

## 8. Why Benjamini–Bogomolov rather than pooling everything

Benjamini & Bogomolov (2014): compute one Simes p-value per family, apply BH
across the K families, and inside each **selected** family apply BH at the
reduced level q·R/K where R families were selected. This controls the expected
average FDR over selected families.

The alternative is to pool all 148,416 tests into one BH. That treats a
granulocyte test and a whole-blood test as exchangeable members of one family.
They are not: different power, different null behaviour, different question.

**The data shows exactly why that matters.** Under global BH:

| family | within-family BH | global BH |
|:--|--:|--:|
| WB_LI | 1,426 | 343 |
| WB_BPD | 4 | **16** |
| CT_LI | 105 | 132 |
| CT_BPD | 0 | **22** |

Pooling *increases* the bipolar discovery counts, from 4 to 16 and from 0 to 22.
Nothing was learned about bipolar disorder to justify that. It happens because
BH's threshold depends on the whole sorted p-value vector, so WB_LI's thousands
of strong signals raise the cutoff for everyone — the lithium contrast lends
significance to the bipolar contrast. That is an artifact of the pooling, and it
is the reason global BH is reported as a cautionary comparison rather than used.

Benjamini–Bogomolov across the four families selects **2 of 4** (both lithium
families; neither bipolar family is selected at all), so the within-family level
becomes 0.025, giving WB_LI 866 and CT_LI 75. The bipolar families report
nothing — which is the correct consequence of having looked at four comparisons
and found signal in two.

---

## 9. Permutation FDR — why it is worth three hours

Every method in §5 except one rests on an assumption about dependence. Gene
expression is co-expressed by construction; PRDS is *plausible* for it but it is
an assumption, and nothing in the data checks it. BY avoids the assumption but
costs a factor of ten in power.

A permutation null avoids it at no such cost. Each permutation scrambles the
group labels and leaves the expression matrix — and therefore every gene-gene
correlation — exactly intact. The expected number of false positives at a
threshold is then **measured**:

> FDR(t) = π₀ · E_perm[ #{null p ≤ t} ] / max(1, #{observed p ≤ t})

the Storey–Tibshirani plug-in, the estimator SAM uses.

**What is permuted:** only the group column of the design matrix. Every other
covariate stays attached to its own subject, so the permutation breaks the
group-expression association while preserving the covariate structure and the
full correlation of the outcome.

**Whole blood** repeats the entire voom + limma fit per permutation, because
voom's precision weights depend on the design matrix and reusing them would leak
the real grouping into the null.

**Cell type** does not repeat bMIND — and that is correct rather than a
shortcut, for the reason in §2: the deconvolution never saw the labels, so a
permutation would reproduce it bit for bit. The precision weights depend only on
the posterior SEs and are likewise computed once. Had the deconvolution used the
phenotype, as `bmind_de` does, each permutation would have cost a 3-hour
deconvolution and this analysis would have been impossible.

**B = 100, not v1's 200.** The quantity needed is a *mean* over permutations,
and a mean converges quickly. v1 used 200 because it needed tail quantiles of
the null for a different question.

---

## 10. Thresholds are columns, not verdicts

0.05 appears in this folder only as a column header. `mtc.R` contains no
threshold at all; `01` reports counts at FDR 0.01, 0.05 and 0.10 and FWER 0.05
side by side so the reader can see how much the answer depends on the cutoff.

Where a single level is needed to *define* a procedure — Benjamini–Bogomolov
needs a q to decide which families are selected — it is a named argument with a
stated default, and the selection it produced is reported in full
(`across_contrast_correction.csv` lists every family's Simes p-value and whether
it was selected), so the step is auditable rather than buried.

---

## 11. What v2 does not change

The per-lineage detection floors stand, and they constrain how every null here
may be read (`../de_analysis/METHODS.md` §17–18):

| lineage | share of blood | smallest detectable effect | as a fold change |
|:--|--:|--:|--:|
| gran | 38.7% | 0.73 – 0.83 | 1.7–1.8× |
| mono | 27.2% | 0.74 | 1.7× |
| T | 24.7% | 0.77 | 1.7× |
| NK | 7.9% | 4.07 | 16.8× |
| B | 1.45% | 26.11 | 72,000,000× |

All five exceed the 0.5 log2FC floor set in advance as the smallest
biologically meaningful effect. **No cell-type zero in this folder is evidence
of absence** — including CT_BPD's zero under every correction method. A better
multiple-testing panel cannot repair a power problem; it can only stop us
overstating the positives.

The same applies to π₀ at the cell-type level: `METHODS.md` §19 showed the
permutation null for π₀ there is centred near 1 but so dispersed (5th percentile
down to 0.073) that a low value means nothing. The π₀ figures in
`summary_by_analysis.csv` are used **only** to power the Storey adjustment,
never quoted as evidence of hidden signal.

---

## 12. Incidents

**Hommel on 12,368 genes × 200 simulations did not finish.** The first version
of the validation script put `p.adjust(., "hommel")` inside a 200-iteration
simulation loop at m = 12,368. Hommel's algorithm is quadratic; the job ran past
ten minutes with no output and was stopped. Corrected by removing Hommel from
the simulation entirely — it is base R, so it was never ours to validate — and
reducing the simulation to m = 8,000 and 150 iterations, which is ample for
estimating a false discovery proportion.

**Two leftover temp files were found by the index.** `DIRECTORY.new` and
`DIRECTORY.part2`, fragments of a failed file splice during the v1 write-up,
were still sitting in the project root. The index build flagged them as
UNDESCRIBED, which is the mechanism working as intended. Both were verified as
fragments and deleted.

---

## 13. What the permutation null actually showed

The point of §9 was to stop assuming the null and measure it. Here is the
measurement. `null_inflation` is the measured number of null p-values below a
threshold divided by the number a uniform null would give.

| analysis | π₀ | observed at p ≤ 0.001 | measured null | uniform null | inflation |
|:--|--:|--:|--:|--:|--:|
| WB_LI | 0.487 | 701 | 16.4 | 12.4 | **1.33×** |
| WB_BPD | 0.643 | 67 | 9.7 | 12.4 | **0.79×** |
| CT_LI | 0.773 | 318 | 85.8 | 61.8 | **1.39×** |
| CT_BPD | 0.802 | 147 | 149.9 | 61.8 | **2.42×** |

**The uniform assumption is close but not exact.** For the whole-blood and
cell-type lithium analyses the real null is about a third heavier-tailed than
uniform — correlation and mild model misspecification doing what they do. That
is small enough that BH is not badly wrong, and large enough that it is worth
knowing rather than assuming.

**CT_BPD is the striking one, in two ways.** Its null is 2.4× inflated, so a
correction that assumes uniformity would be materially anti-conservative there.
And its observed count at p ≤ 0.001 is **147 against a measured null of 149.9** —
the observed p-value distribution is indistinguishable from, in fact very
slightly lighter than, its own permutation null. There is nothing in that
comparison at all, under any correction. The zero is not a thresholding artifact.

**WB_BPD's null is *lighter* than uniform (0.79×)**, which is the conservative
direction and means its 4 genes are, if anything, slightly understated by BH.
They still do not survive Rubin pooling (v1 §14), so this does not revive them.

### Discoveries at the 0.05 level, every method

| analysis | tests | Holm (FWER) | BY | BH | Storey | Permutation |
|:--|--:|--:|--:|--:|--:|--:|
| WB_LI | 12,368 | 125 | 303 | 1,426 | 2,568 | 2,757 |
| WB_BPD | 12,368 | 0 | 0 | 4 | 5 | 32 |
| CT_LI | 61,840 | 22 | 33 | 105 | 120 | 215 |
| CT_BPD | 61,840 | 0 | 0 | 0 | 0 | 0 |

The spread across methods is the honest headline: WB_LI ranges from 125 genes
to 2,757 depending only on which error rate you choose to control. That is why
this folder reports all of them instead of one.

**A caveat on the permutation column.** The permutation q-value is a step
function on the 11-point threshold grid: a gene is assigned the FDR of the grid
point at or below its p-value, not the FDR at its p-value exactly. That makes
the per-gene column slightly optimistic and explains why it exceeds Storey
despite resting on a heavier null. **`permutation_fdr.csv` is the primary
output** — it is exact at each grid point — and the per-gene q-values in
`permutation_qvalues.rds` are a convenience with grid granularity.

---

## 14. Two bugs in this folder's own diagnostic code, found and fixed

Both were in code written here, both were caught by reading the output rather
than by a test, and both are recorded because the first printed numbers that
looked entirely plausible.

**π₀ counted twice in the null-inflation diagnostic.** `n_null_if_uniform` was
computed as `pi0 * m * t`. The permutation count it is compared against is the
mean over *all* m tests, i.e. what happens if every gene is null, so the
like-for-like uniform figure is `m * t` with no π₀. Including it made the null
look 2–4× inflated when the true figure is 0.79–2.42×. **The FDR estimates
themselves were never affected** — π₀ belongs in the FDR formula, where the
question is how many of the real nulls fall below t, and it was correct there.
Only the diagnostic column was wrong.

**The per-gene permutation q-value mishandled both ends of the grid.** A p-value
below the smallest grid point (10⁻⁶) was assigned q = 1 — so the *most*
significant genes in the study got the worst possible q-value — and a p-value
above the largest grid point was assigned the FDR at that last point, crediting
genes the grid never calibrated. Fixed to return `fdr[1]` below the grid
(conservative and correct, since FDR is monotone in t) and 1 above it. The
correction moved WB_LI from 2,655 to 2,757 and CT_LI from 189 to 215.

Re-running cost nothing: `03` checkpoints each analysis, so the corrected run
reloaded all four permutation sets from disk and only recomputed the summary.

---

## 15. Expression-filter sweep

31 filters: 30 count/proportion combinations (>C counts in >=P of samples, C in
{0,1,5,10,20,50}, P in {50,75,90,99,100}%) plus `edgeR::filterByExpr` as a
reference we did not choose. Gene universes from **2,087 to 24,084**. The filter
runs before TMM every time, since scaling factors are computed from the genes
that survive.

### Why a DEG count cannot be the headline

Changing the filter moves two things at once: the denominator of any percentage,
and the number of tests the correction divides by. 1,426 of 12,368 (11.5%) and
1,130 of 24,084 (4.7%) are not 2.4-fold apart in signal — they are the same
finding over different universes. The percentage ranges **4.69% to 29.13%**
across the sweep while the count ranges only 608 to 1,453, and neither number is
interpretable alone.

So the result is **stability**: which genes get called, and how far their ranks
move.

### Whole blood, lithium (unadjusted)

| | range across 31 filters |
|:--|:--|
| DEG at BH 0.05 | **608 – 1,453** |
| as % of genes tested | 4.69% – 29.13% |
| π₀ | 0.360 – 0.676 |

The count is **flat from ~8,000 genes upward** (1,217–1,453, a 19% spread over a
3-fold change in universe). It falls only at the two harshest filters
(>50 in 100% of samples, 2,087 genes → 608 DEGs), and that is not instability —
those filters have deleted genes that were genuinely differentially expressed.

π₀ rises monotonically with the gene count (0.36 → 0.68), which is exactly what
should happen: looser filters admit lowly-expressed genes that are mostly null.

`edgeR::filterByExpr` lands at 13,925 genes and 1,452 DEGs, within 2% of our
baseline. Our hand-set filter and the field-standard automatic one agree.

### Per-gene stability — the part worth keeping

| analysis | significant under ≥1 filter | under **every** filter it was tested in | under exactly one |
|:--|--:|--:|--:|
| LI_raw | 2,594 | **1,092 (42.1%)** | 311 |
| LI_ilr | 569 | 309 (54.3%) | 67 |
| BPD_raw | 11 | **0** | 0 |
| BPD_ilr | 3 | **0** | 0 |

**The lithium result has a hard core of 1,092 genes** called under every filter
that retained them. The ten strongest by median rank — TSPAN2, PIGB, ARMC2,
RFX2, LINC00877, FAR2, MOK, MYADM, MIR24-2, SLC31A2 — are the same genes the
baseline analysis put on top, and TSPAN2 is rank 1 in 30 of 30 filters.

**Not one bipolar gene survives every filter.** The best is C6orf163 at 17 of 24
(70.8%), and the four "baseline DEGs" behave the same way: WDR1 18/31, SMIM8
17/25, AC104820.2 17/24. Added to their disappearance under Rubin pooling and
their failure at any fold-change threshold, the filter sweep is a third
independent reason not to report them.

**Read `n_filters_tested` before `frac_sig`.** A gene is only tested under
filters that retain it, so `frac_sig = 1.000` means different things at
different expression levels: TSPAN2 is 30/30, but SPAG6 is 5/5 — perfect, and
far weaker evidence. The robustness table carries both columns for this reason.

### Overlap with the baseline

Median Jaccard of each filter's DEG set against the baseline's: **0.776**.
It exceeds 0.9 for the six filters nearest the baseline in size and falls off
symmetrically in both directions — to 0.19 at the harshest filter (which simply
has 608 genes to offer) and 0.57 at the loosest (where a larger multiple-testing
burden costs borderline genes). The decline is a property of the universes, not
disagreement about the biology.

### Reporting consequence

The baseline filter (>10 counts in ≥90% of samples) is **not a load-bearing
choice** for the lithium result and is not rescuing the bipolar one. Report the
1,092-gene robust core alongside the baseline 1,426, and treat any gene called
under only one filter as provisional — those 311 genes have a median rank of
2,985 and a median |log2FC| of 0.080, which is below the whole-blood detection
floor of 0.203.

---

## 16. Cell-type threshold sweep

Three sweeps, 18 settings. The gene-filter one can only go in one direction:
bMIND's output exists for the baseline 12,368 genes and no others, so looser
filters would need another 3.1-hour deconvolution. Nine nested (stricter)
filters are available, from 2,087 to 12,368 genes.

Cell type gets a second kind of threshold that whole blood does not, because
its values are **estimates rather than measurements**: an across-sample variance
floor, and a posterior-SE quantile.

### A. Gene filter (9 nested settings)

Lithium, unadjusted, genes DE per lineage:

| genes | gran | mono | T | NK | B |
|--:|--:|--:|--:|--:|--:|
| 2,087 | 6 | 3 | 3 | 0 | **0** |
| 4,469 | 3 | 7 | 12 | 1 | **0** |
| 6,560 | 35 | 12 | 17 | 1 | **0** |
| 8,267 | 59 | 16 | 22 | 2 | **0** |
| 10,742 | 92 | 19 | 25 | 4 | **0** |
| **12,368** | **109** | **20** | **28** | **4** | **0** |

Counts scale with the universe, as they must — a stricter filter has fewer genes
to find. The ordering gran > T > mono > NK > B holds at every setting above
6,560 genes. **B cells return zero under all 18 settings in all three sweeps.**

### B. Variance floor — the shrinkage, measured

bMIND shrinks each sample toward a per-gene prior, and §4 of the record showed
the effect is wildly uneven. Imposing a floor on across-sample variance shows
how uneven:

| floor | gran kept | NK kept | B kept |
|:--|--:|--:|--:|
| none | 12,368 | 12,368 | 12,368 |
| > 0.005 | 12,368 | 12,366 | **5,521** |
| > 0.01 | 12,368 | 12,175 | **452** |
| > 0.05 | 12,367 | 1,992 | **0** |
| > 0.1 | 12,339 | 310 | **0** |

A floor that removes **29 of 12,368 granulocyte genes removes every B-cell gene**
and 97.5% of NK genes. The rare lineages have been shrunk nearly flat; there is
almost nothing left varying between subjects for a group difference to live in.

**The granulocyte result is completely insensitive to this.** 109 DEGs at every
floor from 0 to 0.1. The pooled count actually *rises*, 66 to 95 unique genes,
because discarding flattened genes reduces the multiplicity burden.

### C. Posterior-SE quantile — a caution

Keeping only the fraction of genes bMIND estimates most precisely:

| kept | gran | mono | T | NK |
|:--|--:|--:|--:|--:|
| 100% | 109 | 20 | 28 | 4 |
| 90% | 86 | 12 | 16 | 2 |
| 75% | 78 | 5 | 6 | 1 |
| 50% | **35** | 3 | 4 | 1 |

If findings were independent of precision, halving the gene set should retain
roughly half of them — and BH over fewer tests is *more* powerful, so slightly
more than half. We retain **32%** (35 of 109). The cell-type lithium findings
are therefore somewhat **enriched among the genes bMIND is least confident
about**.

Two readings, and the data here does not separate them. The benign one: posterior
SE rises with across-sample variance, and variable genes are exactly where a
group difference is detectable, so the enrichment is expected. The unwelcome one:
some findings are driven by noisy estimates. Worth stating rather than resolving,
and a reason to prefer the variance floor over the SE quantile as a reliability
filter — the variance floor leaves the granulocyte result untouched, while this
one halves it.

### D. Per-gene stability

**128 of 163** gene × lineage pairs significant somewhere are significant under
**every setting that tested them**. By lineage:

| lineage | always significant | sometimes |
|:--|--:|--:|
| gran | **79** | 30 |
| T | 27 | 3 |
| mono | 18 | 2 |
| NK | 4 | 0 |
| **B** | **0** | **0** |

The most robust findings are the same genes the whole-blood sweep put on top —
TSPAN2 (significant in gran, mono and T under all 14 settings each), PIGB, RFX2,
IL8, MIR24-2, LINC00877. TSPAN2 holds rank 1 in granulocytes, monocytes and
T cells simultaneously.

### E. Bipolar

Five gene × lineage pairs significant anywhere across 18 settings, **none under
more than one**, and zero at the baseline. The isolated hits — 3 T-cell genes at
the 4,469-gene filter, 2 at the 2,087-gene filter, 2 monocyte genes at variance
floor 0.1 — appear at one setting each and are what a 5% false discovery rate
produces when run 72 times.

### Conclusion

The cell-type lithium result does not depend on either threshold: granulocytes
give 109 genes at every variance floor and the lineage ordering is stable above
6,560 genes. The cell-type bipolar null does not depend on them either, but
that remains uninterpretable for the reason in §11 — every lineage's detection
floor exceeds the meaningfulness threshold, and no amount of threshold-varying
repairs a power problem.

The B-cell zero now has a mechanism attached to it rather than just a detection
floor: a variance floor that costs granulocytes 29 genes costs B cells all
12,368.

---

## 18. Narrowing to one estimator and one filter

Reporting four parallel variants of everything is unreadable. This section fixes
one of each and records what that choice costs.

### Estimator: voom

**Chosen because it was pre-specified**, in `config.yaml`
(`whole_blood_estimator: voom`) before any analysis ran. That matters more than
any property of the method, because voom also happens to return the *most* genes
of the three (1,426 vs 1,305 and 1,341) and picking the winner after the fact
would be exactly the selection this project criticises. It was not picked after
the fact.

The three estimators agree almost completely on the lithium contrast:

| | log2FC Pearson | p-rank Spearman | DEG shared | % of smaller set |
|:--|--:|--:|--:|--:|
| voom vs limma-trend | **0.9936** | 0.9790 | 1,250 | 95.8% |
| voom vs edgeR-QLF | **0.9845** | 0.9712 | 1,221 | 91.1% |
| limma-trend vs edgeR-QLF | 0.9840 | 0.9682 | 1,186 | 90.9% |

Of each method's own top 100 genes, **93 are shared by all three**; of the top
50, 42; of the top 500, 430. The estimator choice moves the DEG count by ~9%
and the gene identities barely at all.

**The bipolar contrast behaves oppositely, and that is informative.** Effect
estimates still correlate (0.96–0.99) and ranks still correlate (0.91–0.97), but
the called sets barely overlap: voom 4, limma-trend 0, edgeR-QLF 10, with voom
and QLF sharing 3 and only 51 of each top 100 shared by all three. When sets are
this small and this close to the boundary, which genes land inside is decided by
noise. Three methods disagreeing about a 4-gene list is what nothing looks like,
not a discrepancy to reconcile.

### Filter: >10 counts in ≥90% of samples (12,368 genes)

Three independent reasons, none of which is "it gave the most genes":

1. **Pre-specified** in `config.yaml` before any analysis.
2. **Comparability with the source paper.** Krebs et al. report 12,344 genes
   from the same rule; we get 12,368 on our 474-sample cohort. Using a different
   filter would break the one external check available.
3. **Independently corroborated.** `edgeR::filterByExpr` — an automatic,
   design-aware filter we did not choose — keeps 13,925 genes and returns 1,452
   DEGs, within 2% of our 1,426.

Robustness across the 31-filter sweep (gene universes 2,087 to 24,084):

| metric | value |
|:--|:--|
| DEG count range, all 31 filters | 608 – 1,453 |
| DEG count, filters above ~8,000 genes | **1,217 – 1,453** (19% spread over a 3-fold change in universe) |
| Baseline DEGs called under *every* filter that retained them | **894 of 1,426 (62.7%)** |
| Median Jaccard of each filter's DEG set vs baseline | **0.776** |
| TSPAN2 rank | **1 of 12,368 in 30 of 30 filters** |

The count falls below ~1,200 only at the two harshest settings (>50 counts in
100% of samples, 2,087 genes), where the filter has deleted genes that were
genuinely differentially expressed. That is not instability.

### Composition-adjusted results: a partial disagreement

The request was to drop them. Half of that is right and half would remove the
main finding.

**Dropped:** the standalone ILR-adjusted DEG lists (WB_LI_ilr 341 genes,
CT_LI_ilr 130, and the two bipolar zeros). As independent results they add
nothing — they are the same hypotheses asked a second way, and reporting them in
parallel doubles every table for no gain.

**Kept:** the per-gene *comparison* between adjusted and unadjusted. That
comparison is not a variant of the result, it is the answer to the question the
project exists to ask — how much of lithium's apparent transcriptomic effect is
the cell mixture moving. It produced §17's gradient: 20.5% of all 1,426 genes
survive adjustment, rising monotonically to **100% of the 125 that pass
family-wise error control**. Dropping the adjusted fits entirely would delete
that.

So: one DEG table (unadjusted), plus two derived columns per gene —
`logFC_adjusted` and `% effect lost`.

### Covariates: verified by removal, not by inspection

Naming covariates in a formula is not evidence they were fitted. This already
bit the project once: `bmind_de()` accepts a covariate argument and discards it
(§9). So the check is not "are the columns present" but "does removing them
change the answer".

| | LI | BPD |
|:--|:--|:--|
| design | 226 × 13, **full rank** | 308 × 14, **full rank** |
| covariates used | age, sex, tobacco, RIN, plate, seqPC1–3 | the above + assessment group |
| dropped as constant | none | none |
| `grpcase` is the tested coefficient | yes | yes |
| adjusted vs unadjusted log2FC identical | **FALSE** | **FALSE** |
| correlation between them | 0.960 | **0.654** |
| median \|difference\| in log2FC | 0.016 | 0.038 |
| DEG adjusted vs unadjusted | 1,426 vs 1,500 | **4 vs 54** |

Both designs are full rank with every requested covariate present, and removing
the covariates changes every coefficient — so they are genuinely fitted.

**The bipolar row is worth reporting in its own right.** Without covariates the
contrast returns **54 genes**; with them, **4**. The adjusted and unadjusted
effect estimates correlate only 0.654. Most of the apparent bipolar signal in
whole blood is confounding — age, sex, assessment group, plate and sequencing
position — rather than illness. An uncorrected analysis of this contrast would
have produced a publishable-looking 54-gene result that is mostly demographics.
