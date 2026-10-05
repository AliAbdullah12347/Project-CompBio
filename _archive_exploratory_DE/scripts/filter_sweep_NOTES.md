# filter_sweep — how much of the lithium DE result is a property of the filter?

Source: `scripts/filter_sweep.R`, `_deep.R`, `_null.R`, `_xcheck.R`, `_report.R`.
This file is authored at `scripts/filter_sweep_NOTES.md` and placed here by
`write_run()`.

---

## 1. What was varied

Exactly one thing: the gene-level expression filter. The contrast, the covariates,
TMM normalisation, the voom + eBayes estimator and the Benjamini-Hochberg rule at
FDR 0.05 are the shared engine in `scripts/common.R` and were not modified.

Grid: `min_count ∈ {0, 1, 5, 10, 20, 50}` × `min_prop ∈ {0.10, 0.25, 0.50, 0.75,
0.90, 1.00}` on the **lithium** contrast (bipolar I only, n = 226, 74 non-users
vs 152 users). The rule is the engine's own:

```
keep gene if  rowSums(counts > min_count) >= min_prop * n_samples
```

Note the strict inequality: `min_count = 10` means *more than* 10 counts, which is
the project's canonical definition. `min_count = 0` therefore means "detected at
all", not "no filter".

Added to the grid: `edgeR::filterByExpr()` in three parameterisations, and a
no-filter reference that passes all 57,773 genes.

Baseline reproduced exactly before anything else ran: canonical cell
(10, 0.90) → **12,173 genes, 1,382 DEG, pi0 = 0.5161**, matching the established
baseline. The script asserts the gene count, so a silent upstream change cannot
pass unnoticed.

---

## 2. Design decisions, and why each one

**The engine was not modified.** Every fit went through `de_fit()` with the same
arguments except `(min_count, min_prop)`. If two cells disagree, the disagreement
is the filter and not an accident of how the model was assembled.

**filterByExpr runs through the same engine, not alongside it.** `de_fit()` calls
`filter_genes()` internally and has no "use this exact gene list" argument. Rather
than copy `de_fit()`'s body into this script — which would create a second
estimator free to drift from the first — the script temporarily rebinds
`filter_genes` in the global environment to a function returning the precomputed
filterByExpr set, calls `de_fit()`, and restores the original on exit. The fit is
the same code path as a grid cell; only the surviving gene list differs.

**Every cell is cached to `cache/<tag>.rds`.** Nine sibling R processes were
competing for this machine while the sweep ran, so cells took minutes rather than
seconds. Caching also guarantees the sweep table, the decomposition, the CAT
curves and the figure all come from *the same fit* rather than from four re-runs.

**Cells retaining fewer than 20 genes are skipped**, not fitted. `eBayes` borrows
information across genes; with a handful of genes the moderated variance is
meaningless and TMM has nothing to normalise against.

**Storey pi0 with lambda = 0.5**, as in `common.R`. A DEG count conflates signal
with multiple-testing burden; pi0 estimates the *fraction* of tested genes that
are null, so it is comparable across universes of different sizes. It remains an
estimate, and at low counts what it estimates depends on voom's normal
approximation holding — which is what the permutation null checks.

---

## 3. Four analyses the sweep table alone cannot do

### Fixed-universe rescoring (`summary.csv: n_deg_fixedU`)
A DEG count moves when the filter moves for three separable reasons:

1. the universe changes **size**, so the BH threshold moves;
2. the universe changes **membership**, so different genes are eligible;
3. the **estimates** change, because TMM factors and voom's mean-variance trend
   are both fitted on whatever genes are in the matrix.

Each cell is therefore also scored on the intersection of its universe with the
canonical 12,173 genes, with BH re-applied *inside that intersection*, alongside
the canonical fit scored on that identical subset. That holds (1) and (2) fixed;
anything still moving is (3).

### p-value splicing (`decomposition.csv`)
The decomposition is built directly rather than modelled, by gluing p-value
vectors together and re-running BH:

- canonical p-values alone → the published number;
- canonical p-values + uniform draws for the extra genes → what the loose result
  would be if the extra genes were pure noise (500 draws, so it has an interval);
- canonical p-values + the **real** extra genes' p-values → the burden the real
  extra genes impose, with the canonical genes' own statistics frozen;
- the loose fit restricted to the canonical genes → re-estimation alone.

No model is fitted for any of this; every ingredient is an already-computed
p-value from a cached fit.

### Permutation null, exposure only (`null.csv`)
Under a permuted lithium label there is no lithium effect by construction. Only
the exposure column of the design is permuted (`de_fit(permute_exposure=)`), so
every covariate stays attached to its own sample: the null is "lithium carries no
signal", not "these 226 subjects are exchangeable". The same 8 permutations
(seed 20261218) are used in every cell, so a difference between cells is the
filter and not the draw.

Each permuted fit is scored **separately** on the genes the canonical filter keeps
and on the genes it discards. A well-specified test gives pi0 = 1 and exactly
5.00% of p-values below 0.05 *in every stratum*. A stratum that departs from that
is mis-specified, and apparent signal found there is not evidence.

8 permutations, not more, because the machine was nine-way contended. That is
ample for the question asked — "does the null DEG count stay near zero, and are
the strata calibrated?" — and too few for a precise empirical p-value, so no
empirical p-value is leaned on.

### Rank concordance (`cat_curves.csv`)
Jaccard on FDR sets confounds "the ranking changed" with "the cutoff moved": two
cells can rank genes identically and still share few DEGs if one carries a harsher
BH threshold. The concordance-at-the-top curve — fraction of the top *k* genes by
p-value shared with the canonical top *k*, k = 100…5000 — contains no threshold.

---

## 4. Cross-checks (`xcheck.csv`)

**Method.** voom fits one mean-variance trend across whatever genes are in the
matrix, so a loose filter puts tens of thousands of low-count genes into that
lowess and could move the weights of the well-expressed genes. If the filter
effect were largely a voom artefact, limma-trend and edgeR's quasi-likelihood F
test — whose dispersion estimation is less globally coupled — should show a
smaller effect. Run at the loose and canonical cells for all three estimators.

**Contrast.** The lithium contrast is one 226-sample comparison. The same filter
cells were run on the 474-sample case-control contrast. If the behaviour repeats
it is a property of the filter; if it does not it is a property of this contrast.

---

## 5. Files in this directory

| file | what it holds |
|:--|:--|
| `summary.csv` | one row per filter cell, 30+ scored columns |
| `sweep.rds` | the same plus the canonical gene and DEG lists |
| `strata.csv` | signal by expression decile and by detection rate |
| `decomposition.csv` | the p-value-splicing decomposition of the DEG loss |
| `cat_curves.csv` | rank concordance against the canonical ranking |
| `fixed_nominal.csv` | same genes, fixed nominal cutoff, no multiple testing |
| `deg_character.csv` | expression and effect size of the DEGs the filter decides |
| `null.csv` | permutation null, stratified by expression |
| `xcheck.csv` | method × filter, and the case-control contrast |
| `filter_sweep.png` | four-panel summary figure |
| `cache/`, `cache_x/` | one cached fit per cell (resumability, reproducibility) |
| `_log/` | stdout of the downstream scripts |

---

## 6. Caveats

- Nine sibling R processes competed for CPU throughout. The `secs` column in
  `summary.csv` is wall-clock and is **not** comparable between cells.
- `pi0` at low counts is interpreted only against its own permutation null, for
  the reason given in section 3.
- The no-filter reference is a reference point, not a candidate analysis: roughly
  20,000 of its genes are zero in every sample, so voom's trend there is fitted
  partly on structural zeros.
- Nothing outside `experimentation/` was written. Verified by `find` over
  `Implementation/` — the only recently-modified file outside this folder,
  `data/cohort_474/fractions_bmode_474_CIBERSORTx.csv`, has identical mtime and
  btime (2026-09-28 23:52:26) predating this session's first R call, so it was
  created by the cohort build and not touched here.

---

## 7. Results

### 7.1 The grid

`genes tested`
```
     mp0.1 mp0.25 mp0.5 mp0.75 mp0.9   mp1
mc0  37326  28484 23843  20645 18367 12603
mc1  28812  23586 20673  18351 16482 11111
mc5  20132  18126 16411  14871 13588  8377
mc10 17072  15726 14387  13203 12173  6578
mc20 14541  13554 12566  11556 10532  4483
mc50 11771  10988 10082   9137  8002  2091
```
`DEG at FDR 0.05`
```
     mp0.1 mp0.25 mp0.5 mp0.75 mp0.9  mp1
mc0    881   1053  1155   1247  1269 1420
mc1   1029   1136  1237   1290  1369 1367
mc5   1242   1374  1374   1423  1431 1379
mc10  1399   1398  1432   1435  1382 1261
mc20  1455   1417  1459   1393  1323 1059
mc50  1420   1466  1316   1309  1259  605
```
`DEG as % of genes tested`
```
     mp0.1 mp0.25 mp0.5 mp0.75 mp0.9   mp1
mc0   2.36   3.70  4.84   6.04  6.91 11.27
mc1   3.57   4.82  5.98   7.03  8.31 12.30
mc5   6.17   7.58  8.37   9.57 10.53 16.46
mc10  8.19   8.89  9.95  10.87 11.35 19.17
mc20 10.01  10.45 11.61  12.05 12.56 23.62
mc50 12.06  13.34 13.05  14.33 15.73 28.93
```
`pi0`
```
     mp0.1 mp0.25 mp0.5 mp0.75 mp0.9   mp1
mc0  0.772  0.734 0.679  0.640 0.612 0.522
mc1  0.756  0.677 0.645  0.611 0.574 0.501
mc5  0.648  0.593 0.569  0.548 0.530 0.460
mc10 0.575  0.563 0.546  0.521 0.516 0.437
mc20 0.545  0.529 0.516  0.514 0.504 0.426
mc50 0.511  0.492 0.502  0.481 0.465 0.406
```

### 7.2 The answer to the key question

**The DEG count is not driven by how many genes are tested — and the DEG
percentage is driven by almost nothing else.**

- universe size spans **17.9x** (2,091 to 37,326 genes)
- the DEG count spans **2.42x** (605 to 1,466), and only **1.66x** (881 to 1,466)
  once the degenerate 2,091-gene corner is set aside
- the DEG *percentage* spans **12.3x** (2.36% to 28.93%)

If the count were bookkeeping — more tests, proportionally more hits — then
scaling the canonical result (12,173 genes, 1,382 DEG) to the loosest cell would
predict 1,382 x (37,326/12,173) = **4,238 DEG**. The observed value is **881**,
about one fifth of that, and *lower* than the canonical count despite testing
three times as many genes. The relationship is not merely sublinear, it is
inverted over most of the range.

The count is an **inverted U** in universe size, peaking near 11,000-15,000 genes.
Loosening past that adds genes that almost never clear FDR while raising the
Benjamini-Hochberg burden for the genes that would have; tightening past it starts
deleting genes that carry real signal.

**Practical consequence.** The canonical filter (10, 0.90) returns 1,382 DEG, within
**6%** of the best value found anywhere on a 36-cell grid (1,466 at min_count 50 /
min_prop 0.25). The project's filter choice is not costing it anything.

**Reporting consequence.** Quoting "X% of genes tested were differentially
expressed" is close to meaningless across studies: that number moved 12-fold here
while the underlying data never changed. The absolute count is the far more stable
summary, which is the opposite of the usual intuition.

### 7.3 Relation to the published analysis

Krebs et al. filtered to genes "expressed at >10 counts in 90% of samples" — that
is *exactly* the canonical cell of this grid, so the project is using the published
filter. Their reported lithium DEG count is **976** of 12,344 genes.

976 does lie inside the range this sweep spans (605-1,466), so a filter change
*could* in principle produce a count like it. But because the published filter is
the canonical one, the gap between 976 and this project's 1,382/1,459 is **not** a
filtering effect and must come from the sample set or the model. Two prior
reproductions stored under `Implementation/Krebs Paper/results/` make the same
point from the other side: on the 444-sample QC cohort they give 1,023 DEG
(`DE_lithium_all444`) and 435 DEG (`DE_lithium_casesonly`) from the same 12,353
genes — a 2.4x spread with the filter held completely fixed, wider than most of
what this sweep achieves by moving the filter. (Those files were read, not
re-derived, and their cohort is 444 rather than the 474 used here.)

### 7.4 filterByExpr

| parameterisation | MinSampleSize | as % of 226 | genes | DEG | pi0 | DEG Jaccard vs canonical |
|:--|--:|--:|--:|--:|--:|--:|
| `design =` | 3.93 | 1.7% | 23,246 | 1,120 | 0.735 | 0.601 |
| `group =` lithium | 54.8 | 24.2% | 15,978 | 1,394 | 0.573 | 0.796 |
| neither | 161.2 | 71.3% | 13,789 | 1,459 | 0.529 | 0.879 |
| canonical (10, 0.90) | — | 90% | 12,173 | 1,382 | 0.516 | 1.000 |

All three retain all 1,382 canonical DEGs *in the universe*, so the differences are
burden and re-estimation, not genes made untestable.

**The `design =` form — the one edgeR's documentation steers you towards when you
have covariates — is the worst of the three here, costing 262 DEG (19%) against
the canonical filter.** The mechanism is worth stating because it is not obvious:
`filterByExpr` turns `min.count = 10` into a CPM cutoff using the median library
size, then demands that CPM in `MinSampleSize` samples, and with `design=` it sets

```
MinSampleSize = 1 / max(hat(design))
```

Here `max(hat) = 0.2544` against a mean of `13/226 = 0.0575` — 4.4x the average —
so `MinSampleSize = 3.93`. **A gene need only be expressed in 4 of 226 samples to
be tested.** The two highest-leverage samples are the RIN = 2.3 sample (the cohort
minimum, |z| = 6.79) and a seqPC2 outlier (|z| = 5.91).

So a single degraded-RNA sample silently sets the gene filter for the whole cohort.
This connects directly to the project briefing's note that the deposited QC flag
does not screen RNA quality and that four QC-passing samples have RIN < 5: that
unscreened sample does not merely add noise, it propagates into how many genes get
tested at all. Recommendation: do not use `filterByExpr(design=)` on this cohort
before an RIN screen is applied.
