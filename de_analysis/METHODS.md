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
correlation — intact. 200 permutations per contrast.

Both are computed. The gap is large and in the predicted direction: the
bootstrap interval is **7–10× narrower** than the distance from the null's
centre to its 5th percentile.

| | bootstrap CI width | null median → 5th pct | ratio |
|:--|--:|--:|--:|
| π₀, LI | 0.0296 | 0.3081 | **10.4×** |
| π₀, BPD | 0.0351 | 0.2530 | **7.2×** |

### Correction — the first summary used the wrong tail

`04` judged π₀ by whether it fell inside the null's **central 95% interval**.
That is two-sided and wrong. Signal can only push π₀ *down*; there is no
alternative under which π₀ is too large. The correct reference is the null's
**5th** percentile. `04b` recomputes it. For BPD this changed the reading:
π₀ = 0.7091 sits above the 2.5th percentile (0.6901) and so appeared "inside
the null", but only 7 of 200 permutations reach a π₀ that low, giving a
one-sided p of **0.0398**.

### The null is bimodal — never quote its mean

195 of 200 permutations give **exactly zero** DEGs; the remaining five give up
to 1,683. The mean (17.1 for LI) describes neither mode. Medians and exact
tail counts are reported instead. The five heavy draws are permutations that
happened to land near a real axis of variation; they are legitimate and are
kept.

### Result

| | obs DEG | null median | p(DEG) | obs π₀ | null π₀ 5th pct | p(π₀) |
|:--|--:|--:|--:|--:|--:|--:|
| **LI** | 1,426 | 0 | 0.0100 | 0.5084 | 0.6919 | 0.0149 |
| **BPD** | 4 | 0 | 0.0100 | 0.7091 | 0.7470 | 0.0398 |

**Both contrasts carry signal distinguishable from chance.** The permutation
therefore does *not* overturn the pre-specified bootstrap verdict for WB_BPD —
it agrees that π₀ is genuinely low, and only disagrees about how dramatically.
WB_BPD remains **INCONCLUSIVE**.

What separates the two contrasts is not whether signal exists but its **size**:
11.53% of genes versus 0.032%, and nothing in BPD survives TREAT at
|log2FC| ≥ 0.1. The coherent reading of BPD is a *weak, diffuse,
transcriptome-wide shift* — thousands of genes moving by less than the 0.203
log2FC detection floor, which is exactly the signature that produces a low π₀
with almost no FDR-significant genes.

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

## 13. bMIND clips its output to the range of whatever matrix it is given

Found by reading the MIND source, not documented in the package. The internal
`bmind()` ends with:

```r
res$A[res$A < min(X)] = min(X)
res$A[res$A > max(X)] = max(X)
```

The per-gene model itself is independent across genes — `profile[i, ]` and
`covariance[i, , ]` are built from gene *i* alone and each gene gets its own
MCMC — but this final clip is computed across **all genes in the call**. Since
`06` ran in chunks of 400 genes for memory reasons, the clip bounds are
chunk-specific. The chunk maximum ranged **9.52 to 18.91**.

Measured (`06b`), over all 29,312,160 posterior means:

| | share | what it is |
|:--|--:|:--|
| at the lower bound | 0.6521% | the bound is exactly 0 in **every** chunk, i.e. a non-negativity floor. Not a chunking artifact — an unchunked run applies the same floor. |
| at the upper bound | **0.0353%** | this one *is* chunk-dependent. This is the real cost of chunking. |

0.035% is negligible, so chunking did not materially change the estimates.
Clipping also concentrates in the **abundant** lineages (gran 1.34%, mono 1.21%)
rather than the rare ones (B 0.03%, NK 0.01%), which is the opposite of the
worry.

This finding is also what licenses `08` to test three effect sizes in one
bMIND run using disjoint gene blocks — the only cross-gene coupling is a clip
that almost never binds. That made it affordable to cover all five lineages
instead of three.

**Interpretive consequence beyond this project:** a cell-type expression value
is never free to lie outside the observed bulk range. A lineage at 1.45% of the
mixture can legitimately sit above or below anything seen in bulk, and bMIND
cannot report that.

---

## 14. Both pre-declared sensitivity analyses, and what they changed

`config.yaml` declared four. `03` ran the estimator swap, `04` the permutation
null. The remaining two are `09`.

### A. Tobacco multiple imputation, pooled by Rubin's rules

20 of 474 subjects had no recorded tobacco status. Analysing only column 1
treats a guess as a measurement. Rubin's rules add the between-imputation
variance to the average within-imputation variance; degrees of freedom use the
Barnard–Rubin correction, because the naive Rubin df can exceed the
complete-data df, which is incoherent.

| | main (column 1) | Rubin-pooled | spread across the 20 |
|:--|--:|--:|:--|
| WB_LI_raw | 1,426 | **1,426** | 1,426 – 1,426 |
| WB_LI_ilr | 341 | **341** | 341 – 341 |
| WB_BPD_raw | 4 | **0** | 0 – 5 |
| WB_BPD_ilr | 0 | **0** | 0 – 0 |

Median fraction of missing information is **0.0000** for the lithium contrast
and 0.0030 for the case-control one; median SE inflation 1.0000× and 1.0015×.
So the imputation is immaterial — as it should be, since tobacco is balanced
against lithium within bipolar I (p = 0.70).

**But it is not immaterial to the four BPD genes.** They sit exactly on the
FDR boundary (adj.P = 0.0489), the count moves 0–5 across single imputations,
and correct pooling returns **zero**. Those four genes are not a robust finding
and are not reported as one.

### B. Depth / RIN screened subset

Krebs et al. state twice that no samples were removed for quality, so the main
analysis keeps all 474 and records the flags only. 461 pass both screens.
The gene filter is deliberately **not** re-run — the hypotheses are percentages
and a percentage over a different gene universe is not comparable.

| | full | screened | band |
|:--|--:|--:|:--|
| WB_LI_raw | 11.530% | 8.692% | many → many |
| WB_LI_ilr | 2.757% | 2.264% | many → many |
| WB_BPD_raw | 0.032% | 0.024% | few → few |
| WB_BPD_ilr | 0% | 0% | few → few |

Every band identical. No verdict depends on the quality screen.

---

## 15. Is the lithium signature a composition signature? Tested, not asserted

The top lithium genes *look* like granulocyte genes. "Look like" is not
evidence, and naming genes from memory is the kind of unchecked claim this
project exists to avoid. The checkable version: if the signature is driven by a
shift in the cell mixture, the genes it flags should be the genes whose
expression tracks the granulocyte fraction across subjects. That is computed
from our own matrix and fractions over all 474 subjects (`11`), Spearman so a
few extreme subjects cannot manufacture it.

- **92.8%** of the 1,426 lithium genes have `sign(logFC)` equal to
  `sign(ρ with granulocyte fraction)`. Chance is 50%; binomial p ≈ 2×10⁻²⁷⁰.
- Median |ρ| is 0.253 for DEGs vs 0.205 for non-DEGs; 33.1% vs 18.0% exceed
  |ρ| = 0.3.
- Adjusting for the four ILR balances shrinks |log2FC| by a median of **40.8%**
  (IQR 27.9–55.0%), and only **292 of 1,426** genes survive.

The shrinkage is a **lower bound** on the compositional share: the balances are
estimated from the same matrix, are only four coordinates on five lineages, and
cannot absorb composition LM22 never resolved.

The four BPD genes show no such pattern (ρ = −0.05 to −0.32), which is a small
point in favour of them not being composition-driven — and also of there being
very little there at all.

---

## 16. Where the effect lives once cell types are separated

`07` showed the same genes (TSPAN2, IL8, PIGB, RFX2, MIR24-2) topping the
significant list in granulocytes, monocytes, T and NK cells, with 100% of the
monocyte, T and NK hits also whole-blood hits. The obvious worry is that bMIND
is not separating cell types at all — that each estimate is a rescaled copy of
the bulk, which would make every cell-type result here meaningless.

That worry is testable, and it is **wrong**. If the lineages were copies their
per-gene log2FCs would be near-perfectly correlated. They are not (`12`):

| LI contrast | gran | mono | T | NK | B | bulk |
|:--|--:|--:|--:|--:|--:|--:|
| gran | 1.000 | 0.210 | 0.074 | −0.020 | 0.575 | 0.650 |
| mono | | 1.000 | 0.087 | 0.000 | 0.212 | 0.166 |
| T | | | 1.000 | 0.072 | 0.044 | 0.074 |
| NK | | | | 1.000 | −0.032 | 0.011 |
| B | | | | | 1.000 | 0.427 |

Maximum between-lineage |r| is 0.575. The lineages carry distinct information.

What is actually happening is that the effect is **concentrated in one
lineage**. Over the 1,426 whole-blood lithium genes (median |log2FC| 0.176):

| lineage | fraction | r with bulk | median \|log2FC\| | ratio to bulk |
|:--|--:|--:|--:|--:|
| **gran** | 0.387 | **0.650** | 0.0941 | **0.262** |
| mono | 0.272 | 0.166 | 0.0156 | 0.008 |
| T | 0.247 | 0.074 | 0.0175 | 0.003 |
| NK | 0.079 | 0.011 | 0.0050 | 0.000 |
| B | 0.015 | 0.427 | 0.0076 | −0.005 |

Monocytes, T and NK have a log2FC ratio of essentially **zero** for genes that
clearly move in bulk: the bulk movement did not originate inside them. The
shared gene names are simply the largest bulk signals surfacing wherever there
is any power at all.

**Caveat that must travel with this result.** "Within granulocytes" is only as
fine-grained as our 5-lineage aggregation. The gran lineage pools neutrophils,
eosinophils and mast cells, so a shift in the *mix* of those — more immature
neutrophils, say — appears here as a within-granulocyte expression change. This
analysis cannot separate those two, and the earlier observation that
eosinophil/basophil markers (IL5RA, CLC) carry part of the signal makes that a
live possibility rather than a formality.

The negative B–T (−0.671) and B–NK (−0.525) correlations in the BPD contrast
are a sum-constraint artifact: in a decomposition that must reproduce one bulk
measurement, one compartment rising forces another down.

---

## 17. Cell-type detection floors, measured — and what they force

`08` spikes an effect of size *L* into lineage *c* where it actually lives — in the
bulk, scaled by that sample's fraction, `bulk += L * p_ci * y_i` — then runs the
whole deconvolution-and-test pipeline on the spiked bulk. Labels are permuted
first. Spiking bMIND's *output* instead would skip the deconvolution and measure
nothing but limma.

Three effect sizes share one bMIND run per lineage, spiked into **disjoint gene
blocks**. This is safe because bMIND's model is per-gene and the only cross-gene
coupling is the clip measured in §13, which almost never binds. It cut the cost
threefold and is what made all five lineages affordable; the previous version
covered three and would have left mono and NK without a floor, i.e.
uninterpretable.

| Lineage | Share | MDE at 80% power | As a fold change | Predicted from share alone | Bulk-equivalent at MDE |
|:--|--:|--:|--:|--:|--:|
| gran | 38.7% | **0.73 – 0.83** | 1.7–1.8× | 0.51 | 0.333 |
| mono | 27.2% | 0.74 | 1.7× | 0.75 | 0.201 |
| T | 24.7% | 0.77 | 1.7× | 0.84 | 0.186 |
| NK | 7.9% | 4.07 | 16.8× | 2.79 | 0.297 |
| B | 1.45% | **26.11** | 72,000,000× | 14.10 | 0.376 |

Whole blood, for comparison: **0.203** log2FC (LI), 0.208 (BPD) — a 15% change.

### Consequence 1 — no cell-type null is interpretable

`config.yaml` fixed 0.5 log2FC as the smallest meaningful effect. **Every lineage
exceeds it.** So "0 genes in B cells" is not evidence of absence; it is evidence
that a B-cell analysis at this sample size is blind. The B-cell figure is not a
near miss — a 72-million-fold change is not a reachable effect.

This applies to our own **CT_BPD = REFUTED** verdict, which must be read as "the
prediction was not borne out", never as "shown absent".

### Consequence 2 — the deconvolution buys almost nothing

Convert each floor back into how much it would move the bulk (last column): the
range is **0.186 to 0.376**, against a whole-blood floor of **0.203**. For a
cell-type analysis to see an effect, that effect must move the bulk by roughly as
much as a plain bulk analysis would have needed anyway. The deconvolution does not
recover hidden signal; it re-expresses the same detection limit in per-lineage
units, and for the rare lineages it is strictly worse.

mono (0.74 measured vs 0.75 predicted) and T (0.77 vs 0.84) match the abundance
arithmetic closely, which is the check that the spike-in is measuring what it
claims. gran, NK and B come out 1.5–1.9× worse than abundance alone predicts —
the deconvolution adds noise on top of the dilution.

---

## 18. Our own spike-in design inflated false positives — found and corrected

`08` returned false-positive counts among its **unspiked control** genes of:

| gran | mono | T | NK | B |
|--:|--:|--:|--:|--:|
| 26/200 (**13.0%**) | 2/200 (1.0%) | 0/200 (0%) | 2/200 (1.0%) | 0/200 (0%) |

13% where 5% was intended, in exactly the lineage whose result matters most.

**The cause is our design, not bMIND.** `08` spikes 600 of 800 genes (75%), so
limma's variance prior — fitted across genes, and with `trend = TRUE` fitted
against average expression — is estimated mostly from genes carrying a real
effect. The whole-blood MDE in `03` did not have this problem: 1,200 of 12,368,
under 10%.

The mechanism is worth recording because it is specific to cell-type spike-ins. A
within-lineage effect reaches the bulk as `L * p_ci * y_i`, which **varies between
subjects** because the fraction does. A plain group indicator absorbs only its
mean, so the remainder stays in the residual and inflates the residual variance of
spiked genes. That tilts the variance-versus-expression trend and makes genuinely
null genes at other expression levels look more significant than they are.

`08b` re-measures granulocytes with the ratio inverted — 200 spiked against 600
controls, one run per effect size so blocks stay disjoint:

| spiked share | MDE at 80% power | mean FP rate |
|:--|--:|--:|
| 75% (`08`) | 0.832 | 13.00% |
| **25% (`08b`)** | **0.731** | **7.63%** |

The false-positive rate roughly halves. The floor moves by 0.10 log2FC, and both
values sit above the 0.5 meaningfulness threshold, so **§17's conclusion does not
depend on which is used**. The residual error runs in a known direction —
anticonservative, i.e. optimistic — so the true floors are if anything higher
than reported, which strengthens rather than weakens the conclusion.

This distortion belongs to the **calibration**, not to the real analysis, where
under 1% of genes carry a large effect and the variance trend is undistorted.

---

## 19. π₀ is not usable at the cell-type level; and the mechanism test

### The problem `07` raised

`07` produced per-lineage π₀ values that are not internally coherent:

| | π₀ | DEGs at FDR 0.05 | smallest p |
|:--|--:|--:|--:|
| BPD, B cells | 0.4395 | 0 | 8×10⁻⁵ |
| LI, B cells | 0.4166 | 0 | 1.3×10⁻⁴ |
| BPD, T cells | 0.5268 | 0 | 3×10⁻⁵ |

π₀ = 0.44 asserts that 56% of genes carry signal. A gene set with 56% true effects
at n = 308 should not produce **zero** significant genes. Reporting those numbers
at face value would mean claiming widespread hidden signal in B cells.

### The test

Because `06` ran bMIND **without** the phenotype, the cell-type estimates do not
depend on the labels at all. Permuting labels afterwards and re-running only the
limma step is therefore a *complete* null — the deconvolution would be identical
and does not need repeating. (With `bmind_de`, which takes y, each permutation
would have cost a 3-hour deconvolution. This is a concrete payoff of §10's design
choice.) 200 permutations per contrast, on a fixed 800-gene subsample drawn once
with the config seed.

### The answer

The null is **centred correctly** — every per-lineage null median π₀ is 0.986 to
1.000, so there is no systematic bias. It is not **stable**: the null's 5th
percentile runs down to **0.073**. Shuffled labels routinely produce π₀ of 0.3 or
lower by chance.

| Contrast | Lineage | obs π₀ | null median | null 5th pct | p |
|:--|:--|--:|--:|--:|--:|
| LI | gran | 0.335 | 0.991 | 0.318 | 0.065 |
| LI | B | 0.425 | 1.000 | 0.073 | 0.244 |
| BPD | T | 0.522 | 1.000 | 0.423 | 0.085 |
| BPD | B | 0.430 | 1.000 | 0.085 | 0.169 |

No lineage in either contrast reaches significance. The striking B-cell figure is
ordinary noise. **π₀ is therefore not reported as evidence at the cell-type
level** — only at the whole-blood level, where §7's permutation showed it behaves.

Gene **counts** behave far better: across both contrasts the 95th percentile of
the null DEG count is 0 or 1 gene. A detection like the 109 granulocyte genes is
not something shuffled labels produce. Detections are trustworthy here; π₀ is not.

### The pre-specified mechanism test for CT_BPD

`config.yaml` declared it, so running it is a commitment kept rather than a
post-hoc search. CT_BPD is otherwise decided on a count comparison, and a count
cannot demonstrate a mechanism. Masking predicts something sharper and checkable:
a gene hidden in bulk *because* two lineages move it in opposite directions should
show opposite-sign effects at nominal p < 0.05.

| | |
|:--|--:|
| observed discordance, 12,364 bulk-null genes | **0.000%** |
| permutation null, median | 0.000% |
| permutation null, 95% range | 0.000% – 0.125% |
| empirical p | 1.000 |
| pre-specified threshold | ≥ 10% |
| **verdict** | **NOT SUPPORTED** |

Not one gene showed the pattern. Bipolar signal is not being hidden in whole blood
by lineages pulling in opposite directions — at least not at a size this study
could see, which §17 shows is not a strong constraint.

---

## 20. Incidents

**Killed the wrong process.** A `taskkill` aimed at the bmind_de loop also
stopped the running permutation job. LI had completed (200 permutations,
retained); BPD was restarted.

**The bmind_de wrapper outlived its R process.** An earlier `taskkill` killed
the R process but not the bash loop driving it, which proceeded to the next job.
Three of four runs completed unattended as a result. They were later discarded
anyway under §9.

**Ran `08` with too many workers and paged to disk.** At `ncore = 6` on a 7.7 GB
machine the run went into the pagefile: 1.5 GB paged, and only ~1.2 of 6 cores
actually busy, making it roughly 5× slower than fewer workers would have been. It
was stopped and restarted at `ncore = 4` with per-lineage checkpointing, which
recovered to 3.96 effective cores and ~5 min per lineage. More parallelism is not
free when each PSOCK worker carries its own copy of the data.

**Two scripts died on R's top-level `else`.** `06b` and `07` both wrote
`x <- if (a) u` followed by `else v` on the next line at top level, which R
parses as a complete expression followed by a syntax error. Both printed their
results and then failed before writing any CSV. Fixed by wrapping each in braces.
The results were unaffected; only the files were missing.

**A wrong summary was published to the user before being caught.** The first
permutation summary judged π₀ against the null's *central 95%* interval, which is
two-sided and wrong for a statistic that signal can only push downward. It was
reported as flipping WB_BPD to SUPPORTED. `04b` corrected it to the one-sided
test; WB_BPD remains INCONCLUSIVE. See §7.
