# Methods and decision log

Every choice made in this analysis, and the reason for it. Written as the work
was done, in order. Where a decision was changed, the original and the reason
for changing it are both kept rather than overwritten.

The pre-specification is `config.yaml`, committed at **2026-09-30 18:05:55**
(`958c1a5`) before any analysis in this folder existed.

---

## 0. Scope and honesty declaration

Exploratory differential expression was done earlier and is archived in
`../_archive_exploratory_DE/`. Nothing here reads it. But the people who wrote
`config.yaml` had seen it, so the four hypotheses are labelled by how much
prior knowledge they carry:

| | Status |
|:--|:--|
| **WB_LI** | **Positive control.** Result already known (~11%). Not a discovery. |
| **WB_BPD** | Confirmatory-of-exploratory. Direction was informed by the archived work. |
| **CT_LI** | **Genuinely pre-specified.** No cell-type analysis had been run. |
| **CT_BPD** | **Genuinely pre-specified.** |

Pretending all four were equally blind would be the exact practice this
project criticises.

---

## 1. Why WB_LI is a control, not a hypothesis

Boltz et al. reported zero cell-type lithium genes and could not distinguish
"no effect" from "no power." A null is only evidence if the same pipeline
demonstrably finds a known effect. WB_LI is that demonstration. If it had
failed, every null below would have been uninterpretable and the correct action
would have been to stop.

---

## 2. Thresholds in percentages, not gene counts

The archived exploratory work resampled 200 of 226 subjects and watched the DEG
count move across an **11-fold range** (251 to 2,894). A decision rule built on
a raw count is therefore fragile.

Two consequences, both written into `config.yaml`:

- Bands are **percentages of genes tested** — which also makes whole blood
  (12,368 genes) comparable with per-lineage sets of other sizes.
- Every count criterion is paired with **π₀**, which uses the whole p-value
  distribution instead of only the genes crossing a line, and is far more
  stable.

Bands are **0.5%** and **2%**. Justification is empirical: every null-like
result on record sits below 0.1% of genes (Krebs' lithium-adjusted
case-control, 6/12,353; Boltz's cell-type lithium, 0) and every clear-signal
result above 10%. The 0.5–2% window is empty space, so the classification does
not depend on where exactly in that window the line falls.

---

## 3. One gene filter, applied once

The filter (>10 counts in ≥90% of samples) runs **once on all 474 samples**,
giving 12,368 genes, and that set is used for every contrast.

This differs from the archived exploratory work, where the filter was re-run per
subset. The change is deliberate: the hypotheses compare percentages **across**
contrasts and levels, and a percentage computed over a different gene universe
each time is not comparable.

---

## 4. Gene lengths by union of exons

TPM needs a per-gene length. Exons of different transcripts overlap, so their
lengths cannot be summed — the intervals are merged first. Union is also the
correct denominator for HTSeq union mode, which generated these counts.

The first implementation grew a per-gene matrix inside the read loop, which is
quadratic in the number of exon records and did not finish. Rewritten to
accumulate flat vectors and split once: 753,811 exon records, 12,368 genes,
3 minutes. TPM columns sum to exactly 1e6.

---

## 5. Five lineages, not 22 LM22 types

Evidence-based, not convenience. At 22-type resolution the deconvolution's
between-method intraclass correlation is **0.389**, and the compositional
distance between two methods exceeds the spread between subjects. Aggregated to
five lineages, median ICC is **0.693** and the myeloid/lymphoid balance reaches
**0.907**.

bMIND takes proportions as *input*, so unreliable fine-grained proportions
propagate into every cell-type expression estimate.

---

## 6. Covariate adjustment — verified, not assumed

Every design was printed and checked rather than trusted
(`scripts/VERIFY_designs.R`):

| Fit | n | Design | Rank | Full rank |
|:--|--:|--:|--:|:--|
| WB_LI_raw | 226 | 226 × 13 | 13 | yes |
| WB_LI_ilr | 226 | 226 × 17 | 17 | yes |
| WB_BPD_raw | 308 | 308 × 14 | 14 | yes |
| WB_BPD_ilr | 308 | 308 × 18 | 18 | yes |

All eight requested covariates appear, plate correctly expands to four dummies,
nothing was dropped as constant or aliased, and the tested coefficient is
`grpcase` in every case.

`assessment group` appears in BPD but not LI, which is correct: all bipolar I
subjects are group A, so it is inestimable in the lithium contrast, but it
varies once controls are included.

**The covariates change the lithium answer very little** — 1,500 unadjusted vs
1,426 adjusted. That is not a failure: within bipolar I, tobacco (p = 0.70) and
plate (p = 0.89) are genuinely balanced against lithium. Age is the one real
imbalance and is adjusted for.

---

## 7. Two uncertainty estimates for π₀, and which to believe

`config.yaml` specifies a **bootstrap over genes**, which is standard. It is
also wrong in a known direction: it assumes genes are independent, and
co-expression is the premise of the whole field. Its interval is therefore too
narrow, which biases the decision rule toward "the CI excludes 1" and hence
toward INCONCLUSIVE.

A **permutation null** does not have that defect. Each permutation scrambles the
group labels while leaving the expression matrix — and every gene-gene
correlation — intact.

Both are computed. The gap is large and in the predicted direction:

| | bootstrap CI | permutation null range |
|:--|:--|:--|
| π₀, LI contrast | [0.493, 0.523] | **[0.606, 1.000]** |

Where they disagree, **the permutation is the one to trust**, and the
pre-specified bootstrap rule should be read as the conservative version.

---

## 8. Minimum detectable effect by spike-in, not formula

A power formula needs an assumed variance. The spike-in uses the actual
variance of this matrix at this sample size under this design. Labels are
permuted before spiking so no real signal contaminates the calibration.

Whole blood: **MDE ≈ 0.20 log2FC at 80% power** for both contrasts — about a
15% change in expression, comfortably below the 0.5 meaningfulness floor. So
whole-blood nulls are interpretable.

---

## 9. bmind_de abandoned — two verified defects

`bmind_de()` was the planned cell-type estimator. It is the function Boltz et
al. used. Four runs completed before two defects were found, both confirmed by
direct test rather than inferred.

### Defect 1 — a hard p-value floor

Minimum p-value across all four runs: **2×10⁻⁴**. Only ~600 distinct p-values
appear across 61,840 tests.

BH at FDR 0.05 over 12,368 genes requires the smallest p-value below
**4.0×10⁻⁶**. The floor is **50× too large**, so no gene can reach significance
however strong the true effect is. The zero DEGs was mathematically forced.

Raising MCMC iterations tenfold (1,300 → 13,000) moved the floor from 0.086 to
0.073 on a test set. It is structural in how the function forms p-values, not a
sampling-resolution problem more iterations fix.

### Defect 2 — covariates silently ignored

`bmind_de(covariate = NULL)`, `bmind_de(covariate = CV)` and
`bmind_de(covariate_bulk = CV)` return **byte-identical** p-values. The argument
is accepted and discarded.

This was first noticed because the `raw` and `ilr` runs — differing by four ILR
covariates — produced identical results to machine precision.

Consequence: every bmind_de result is unadjusted for age, sex, tobacco, RIN,
plate, sequencing PCs **and** the cell-composition balances.

### Implication beyond this project

Boltz et al. used `bmind_de` and reported zero cell-type lithium DEGs across a
1,730-sample cohort. Defect 1 means that null was very likely **unobtainable
rather than observed**. This is worth reporting regardless of what our own
analysis finds.

---

## 10. The replacement

Run `bMIND()` **without** the phenotype to get per-sample, per-lineage
expression estimates, then test with limma in a design we control.

| | Gains |
|:--|:--|
| p-values | continuous, no floor |
| covariates | actually enter the model |
| estimator | same family as whole blood, so the two levels are comparable |
| circularity | the phenotype never informs the deconvolution |

bMIND returns a posterior SE per gene × lineage × sample. These are used as
limma **precision weights** (1/SE²), so the deconvolution's own uncertainty is
propagated rather than discarded. Treating posterior means as measured data
would understate uncertainty and inflate the cell-type DEG count — biasing
CT_BPD toward "supported" exactly where no thumb on the scale is wanted.

One run serves both contrasts. Cost: ~7.6 min per 400-gene chunk at n=474,
≈ 4 hours for 12,368 genes.

---

## 11. Cell-type MDE must be per lineage

Bulk expression is the abundance-weighted average of its cell types, so a
within-lineage effect of size *L* in lineage *c* moves the bulk by only
*L × p_c*.

Granulocytes are 38.7% of this mixture; B cells are 1.45%. The same biological
effect inside a B cell is **~27× harder to detect**, before any statistical
consideration. Predicted MDE from abundance alone:

| Lineage | Mean fraction | Predicted MDE (bulk MDE ÷ fraction) |
|:--|--:|--:|
| gran | 0.387 | 0.52 log2FC |
| mono | 0.272 | 0.75 |
| T | 0.247 | 0.82 |
| NK | 0.079 | 2.59 |
| **B** | **0.015** | **14.0** |

A B-cell MDE of 14 log2FC is a 16,000-fold change — not a biologically
reachable effect. **"No B-cell genes" and "no granulocyte genes" are therefore
completely different statements** and must not be reported as one number.

`08_celltype_mde.R` measures this empirically by spiking the effect where it
actually lives — in the bulk, scaled by each sample's lineage fraction — and
running the whole deconvolution-and-test pipeline on the spiked data. Spiking
bMIND's *output* instead would skip the deconvolution and measure nothing but
limma.

---

## 12. Two multiplicity corrections, for two different questions

- **Per-lineage BH** answers "is there signal in granulocytes?" Each lineage is
  its own family of ~12,368 tests.
- **Pooled BH** over all 5 × 12,368 = 61,840 tests answers "is there signal
  anywhere?", and is the only fair basis for comparing the cell-type level
  against whole blood — otherwise the cell-type side gets five independent
  chances and wins on multiplicity alone.

The cell-type percentage is computed over **unique genes out of 12,368**, not
over the 61,840 gene × lineage tests, so that it measures the same thing as the
whole-blood percentage.

---

## 13. Incidents

**Killed the wrong process.** A `taskkill` aimed at the bmind_de loop also
stopped the running permutation job. LI had completed (200 permutations,
retained); BPD was restarted.

**The bmind_de wrapper outlived its R process.** An earlier `taskkill` killed
the R process but not the bash loop driving it, which proceeded to the next job.
Three of four runs completed unattended as a result. They were later discarded
anyway under §9.
