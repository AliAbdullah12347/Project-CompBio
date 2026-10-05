# cellcomp — is the lithium expression signature a cell-composition shift?

Everything below was computed in this experiment. No number is carried forward
from a paper or from a briefing note without being recomputed here.

Scripts, in the order they must run:

| script | what it does |
|:--|:--|
| `scripts/cellcomp.R` | the six models, the permuted-ILR null, path a, the loss-by-loading breakdown, specificity |
| `scripts/cellcomp_power.R` | collinearity, first spike-in, estimator sensitivity, survivor list |
| `scripts/cellcomp_spike2.R` | effect-size sweep for the spike-in (and the diagnosis of why the first COMPOS arm failed) |
| `scripts/cellcomp_spike3.R` | the corrected, variance-preserving spike-in — this is the one to quote |
| `scripts/cellcomp_mediator.R` | interrogates the mediator itself: which cell types, deconvolution bias, lithium vs illness |
| `scripts/cellcomp_annot.R` | gene identity — direction, marker enrichment, how far the lost genes fell, named survivors |

`data/gene_annot.csv` is built once from the GENCODE v19 GTF in
`../cibersortx/data/` (read-only). It has 57,820 rows, which is exactly the
GENCODE v19 gene count, so the annotation is pinned to the same build the
counts were generated against.

---

## 1. The question and the design

Lithium raises neutrophil proportion. Bulk blood RNA is a weighted average over
cell types, so moving the weights moves thousands of genes without any single
cell changing what it is doing. The lithium contrast (bipolar I only, n = 226,
74 non-users vs 152 users) was fitted six ways. Every argument to `de_fit()`
except `adjust` was held at the canonical value, so each row differs from the
canonical row by the composition term and nothing else.

| id | composition term | rationale |
|:--|:--|:--|
| a | none | canonical baseline |
| b | 4 ILR balances | the pre-specified compositional model |
| c | 4 of the 5 lineage proportions | one must be dropped: they sum to 1 |
| d | first 5 PCs of the 22 fractions | 95.1% of fraction variance |
| e | ILR balance 1 alone | myeloid vs lymphoid |
| e2 | logit(granulocyte share) | the number a routine blood count reports |

Choices worth recording:

- **PCA is centred but not scaled.** Four LM22 types are exactly zero in all
  474 samples of this CIBERSORTx b-mode run, so scaling divides by zero. Those
  four columns were dropped first, leaving 18.
  *Note:* the four all-zero types here are T follicular helper, NK activated,
  macrophages M2 and mast activated. The project briefing records a different
  set of four structural zeros (T follicular helper, T regulatory, macrophages
  M1, mast activated) for the **deposited classic CIBERSORT** fractions. These
  are two different deconvolution runs, so a different zero set is expected,
  but the briefing's list should not be quoted for the CIBERSORTx b-mode data.
- **B cells are the omitted lineage in model c**, chosen as the smallest.
- All models run on the identical gene set (12,173) and sample set (226); this
  is asserted in code, because "how many canonical DEGs survive" is meaningless
  if the two models tested different genes.

---

## 2. Headline result

12,173 genes, 226 samples, FDR 0.05. The canonical model reproduces the stated
baseline exactly (12,173 genes, 1,382 DEG, π₀ = 0.516).

| model | n DEG | π₀ | canonical survive | canonical lost | **prop lost** | new | r(logFC) all | shrink survivors |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| a canonical | 1382 | 0.516 | 1382 | 0 | — | 0 | 1.000 | 1.000 |
| b ILR ×4 | 350 | 0.763 | 301 | 1081 | **0.782** | 49 | 0.801 | 0.696 |
| c lineage ×4 | 323 | 0.765 | 255 | 1127 | 0.816 | 68 | 0.758 | 0.699 |
| d fraction PC ×5 | 312 | 0.754 | 243 | 1139 | 0.824 | 69 | 0.743 | 0.723 |
| e ILR1 only | 408 | 0.716 | 399 | 983 | 0.711 | 9 | 0.936 | 0.775 |
| e2 logit(gran) | 260 | 0.814 | 229 | 1153 | **0.834** | 31 | 0.802 | 0.710 |

**78.2% of the canonical lithium DEGs do not survive adjustment for cell
composition**, and the ones that do survive keep only ~70% of their effect.
The median |logFC| over the canonical DEG set falls from 0.180 to 0.096.

π₀ rises from 0.516 to 0.76 — the signal is not merely pushed below a
threshold, the whole p-value distribution flattens towards the null.

Only 2 genes out of 1,382 change sign (in model d), so this is attenuation
towards zero, not a reversal.

---

## 3. The result that makes it interpretable

Losing DEGs on adjustment is not by itself evidence of confounding, because
adjustment also costs power. Three controls separate the two.

### 3.1 Collinearity is negligible

Regressing lithium on the ILR balances given the other covariates:

- R² = **0.0384** — composition explains under 4% of lithium's variance
- variance inflation 1.040, so **standard errors grow by only 2.0%**
- residual df falls from 213 to 209

A 2% increase in standard errors cannot remove 78% of a DEG list.

### 3.2 A permuted-ILR null prices the degrees of freedom

50 refits in which the rows of the ILR matrix are shuffled and relabelled. Same
four columns, same marginal distribution, same 4 degrees of freedom, no true
link to any sample. Whatever this costs is the price of the degrees of freedom
alone.

| | permuted ILR (n = 50) | real ILR |
|:--|--:|--:|
| DEGs | 1379 ± 189 (range 884–1952) | 350 |
| canonical DEGs surviving | 1279.9 | 301 |
| proportion lost | **0.074** | **0.782** |
| r(logFC) vs canonical | 0.998 | 0.801 |
| shrinkage | 0.9998 | 0.696 |

**Four degrees of freedom of pure noise cost 7.4% of the DEG list. The real
composition costs 78.2%.** The excess — 979 genes, 70.8% of the canonical list
— is attributable to genuine confounding. z = 9.0; no permutation of 50 came
close (permutation p = 0.020, which is the floor at 50 draws).

Note also that the permuted null shrinks logFC by 0.02% and leaves the
correlation at 0.998. Adjustment per se does not attenuate anything.

### 3.3 Spike-in calibration: the decisive control

Known effects were planted in the real expression matrix, drawn from the
observed canonical DEG |logFC| distribution (median 0.180, IQR 0.131–0.249) so
the synthetic signal is the same size as the real one.

- **TRANS** — a purely transcriptional lithium effect, carrying no composition
  component. This is the ceiling: what adjustment costs an honest signal.
- **COMPOS** — an effect entering only through the myeloid/lymphoid balance,
  so the gene has no lithium term at all and its apparent lithium association
  is entirely inherited. This is the floor.

TRANS results (`spikein_calibration_curve.csv`), 1,200 genes per multiplier:

| multiplier | detected canonically | power | survive adjustment | survival | shrinkage |
|--:|--:|--:|--:|--:|--:|
| 1.0 | 879 | 0.73 | 797 | 0.907 | 1.013 |
| 1.5 | 1085 | 0.90 | 1046 | 0.964 | 1.007 |
| 2.0 | 1150 | 0.96 | 1124 | 0.977 | 1.019 |
| 3.0 | 1195 | 1.00 | 1190 | 0.996 | 0.999 |
| 4.0 | 1197 | 1.00 | 1197 | 1.000 | 0.994 |
| 6.0 | 1198 | 1.00 | 1197 | 0.999 | 1.000 |

**A genuine transcriptional effect of the same size as the observed one
survives composition adjustment 91–100% of the time, and is not attenuated at
all** (shrinkage 0.99–1.02 throughout). The real signature survives 21.8% of
the time and shrinks to 0.70.

Corrected, variance-preserving version (`cellcomp_spike3.R`, see §3.4 for why):
TRANS ceiling **0.944** over 4,746 detections, mean shrinkage **0.998**.

Reading the real signature against the anchors:

| construction | COMPOS floor | REAL | TRANS ceiling | implied % compositional |
|:--|--:|--:|--:|--:|
| additive (`spike2`) | 0.071 | 0.218 | 0.977 | **84%** |
| variance-preserving (`spike3`) | 0.000 | 0.218 | 0.944 | **77%** |

Two different constructions, agreeing to within 7 points. The COMPOS floor is
thinly estimated in both (14 and 1 detections respectively — see §3.4), so the
TRANS ceiling is the number to lean on: it is estimated from thousands of genes
and it says adjustment costs an honest signal 3–9%, not 78%.

### 3.4 A mis-specification worth recording

The first COMPOS arm (`cellcomp_power.R`, and the sweep in
`cellcomp_spike2.R`) detected almost nothing — 2 to 4 genes out of 1,200 at
every multiplier tried, up to 6×. That is not low power to be fixed with a
bigger effect; it is a construction error, and an easy one to make.

Injecting an effect of size *d* through the balance needs a coefficient of
*d*/a₁ = *d*/0.156 = 6.4*d* on the balance. Adding that term **on top of** the
gene's existing residual inflates the gene's noise far more than it inflates
its lithium signal, so the simulated gene is simply a noisy gene. A real
composition-driven gene is not like that: its expression *is* the
mixture-weighted average, so composition explains variance the gene already
had rather than adding new variance.

`cellcomp_spike3.R` corrects this by taking the injected signal out of the
gene's own residual budget, holding total residual variance fixed, and applies
the same treatment to both arms so the comparison stays symmetric. The TRANS
arm behaves well under both constructions (0.914 vs 0.907 survival at the
matched 1× multiplier), which is reassuring.

**But the corrected COMPOS arm degenerates too, and more informatively.** At
1× it had to cap 1,019 of 1,200 genes: the variance the balance would have to
explain *exceeds the gene's entire residual variance*. Working the arithmetic
through — a1 = 0.156, so an induced logFC of 0.18 needs a balance coefficient
of 1.16; sd(b1) = 0.536, so that implies a composition-explained variance of
1.16² × 0.536² = 0.386, against a median gene residual variance of 0.323² =
0.104, nearly four times larger.

This is not a failed simulation so much as a statement about the data: **a
gene cannot have a median-sized lithium effect running entirely through b1
unless composition dominates its variance.** Which raises the obvious question
— how did 1,382 genes get detected in the real data, then? The answer must be
that the real casualties are exactly the genes where composition *does*
dominate, unlike the randomly chosen genes the spike-in used. That is
measurable rather than simulable, and §3.5 measures it.

So: spike2 and spike3 are both kept. The TRANS arms are the usable result; the
COMPOS arms' failure is itself the diagnosis, and it pointed at §3.5.

---

### 3.5 The drop in logFC *is* the mediated path — checked, not assumed

For a linear model the difference-of-coefficients and product-of-coefficients
definitions of mediation are algebraically the same quantity. So for every gene
the drop in the lithium coefficient should equal the sum over balances of
(lithium → balance) × (balance → gene):

> **cor(Σ aₖ·b_gk, logFC_canonical − logFC_adjusted) = 0.9965, slope 1.015**

They are not forced to agree, because voom reweights between the two fits, so
0.9965 is a real check on both. It means the 78% is not a thresholding
coincidence: the amount each gene's effect falls is, gene by gene, exactly the
size of its composition-mediated path.

### 3.6 Where the losses fall — and an honest negative

Canonical DEGs binned by how strongly they load on the myeloid/lymphoid
balance (quintiles of |b1 coefficient|, 2,435 genes each):

| quintile | median abs(b1) | canonical DEGs | lost | % lost |
|:--|--:|--:|--:|--:|
| Q1 (weakest) | 0.050 | 121 | 84 | 69.4 |
| Q2 | 0.152 | 193 | 148 | 76.7 |
| Q3 | 0.260 | 231 | 194 | 84.0 |
| Q4 | 0.386 | 313 | 254 | 81.2 |
| Q5 (strongest) | 0.622 | 524 | 401 | 76.5 |

Two things, and the second is a negative result I did not expect.

- **Where DEGs live is strongly compositional.** Each quintile holds 2,435
  genes, so a flat signature would put ~276 DEGs in each. Instead Q5 holds 524
  and Q1 holds 121 — a 4.3-fold enrichment of lithium DEGs among the genes
  most tied to the myeloid/lymphoid balance. That is the composition signature
  in its rawest form.
- **But the loss rate conditional on being a DEG is essentially flat** (69–84%,
  no monotone trend). I expected a gradient and there is not one. The
  explanation is that b1 is one of four balances and not the most relevant one
  (see §6) — a gene with a weak b1 loading can have a strong b2 loading, and
  the adjustment uses all four. Quintiles of |b1| are therefore a poor proxy
  for total composition dependence. §3.7 measures that properly instead.

## 4. Mediation path a — lithium does move the blood

*(`cellcomp_mediator.R`; same 226 samples, same covariates as the DE models)*

### 4.1 It is neutrophils, and essentially only neutrophils

Of the 18 LM22 types with non-zero variance, **exactly one clears FDR 0.05**:

| cell type | non-users | users | difference | t | p | FDR |
|:--|--:|--:|--:|--:|--:|--:|
| **Neutrophils** | 0.346 | 0.392 | **+4.64 pct pts** | 3.16 | 0.0018 | **0.032** |
| T CD4 memory resting | 0.116 | 0.100 | −1.57 | −2.63 | 0.0093 | 0.084 |
| NK resting | 0.081 | 0.069 | −1.15 | −2.01 | 0.045 | 0.268 |
| Monocytes | 0.273 | 0.261 | −1.25 | −1.46 | 0.146 | 0.380 |

This is the strongest possible form of the result: the shift is not a diffuse
smear across unrelated cell types, which is what a deconvolution artifact would
look like. It is specifically neutrophilia, which is exactly lithium's
documented haematological effect.

At the coordinate level — all 14 candidate mediators, BH-corrected across the
14 (`path_a_mediator_models.csv`), sorted by evidence:

| mediator | β | β in SD | t | p | **FDR** | partial R² |
|:--|--:|--:|--:|--:|--:|--:|
| granulocyte proportion | 0.0407 | 0.456 | 3.21 | 0.0015 | **0.013** | **0.046** |
| neutrophil fraction (raw) | 0.0394 | 0.450 | 3.16 | 0.0018 | **0.013** | 0.045 |
| logit(granulocyte) | 0.1720 | 0.420 | 2.94 | 0.0036 | **0.014** | 0.039 |
| log neutrophil/lymphocyte | 0.1922 | 0.416 | 2.89 | 0.0042 | **0.014** | 0.038 |
| lymphocyte fraction (raw) | −0.0302 | −0.410 | −2.85 | 0.0049 | **0.014** | 0.037 |
| **b2 (gran vs mono)** | 0.1026 | 0.360 | 2.55 | 0.0116 | **0.027** | 0.030 |
| T proportion | −0.0206 | −0.340 | −2.34 | 0.0201 | **0.040** | 0.025 |
| NK proportion | −0.0077 | −0.277 | −2.01 | 0.0452 | 0.073 | 0.019 |
| **b1 (myeloid vs lymphoid)** | 0.1558 | 0.291 | 2.00 | 0.0471 | **0.073** | 0.018 |
| monocyte proportion | −0.0105 | −0.210 | −1.49 | 0.137 | 0.186 | 0.010 |
| B proportion | −0.0019 | −0.177 | −1.19 | 0.236 | 0.275 | 0.007 |
| b3 (T vs NK+B) | 0.0354 | 0.082 | 0.56 | 0.578 | 0.588 | 0.002 |
| b4 (NK vs B) | 0.0629 | 0.081 | 0.54 | 0.588 | 0.588 | 0.001 |

**Path a is significant — the mechanism is real.** But note which row is on
top and which is not: the pre-specified primary balance b1 does **not** clear
FDR 0.05 across the 14 mediators (FDR 0.073), while granulocyte proportion
clears it at 0.013 with 2.5× the partial R². See §6.

Lithium moves the blood by roughly 0.45 SD on the granulocyte axis. The
lymphocyte side moves the mirror image (−0.41 SD), as a closed composition
must.

### 4.2 The composition shift is lithium, not illness

| contrast | granulocyte proportion | p |
|:--|--:|--:|
| lithium within BP1 (n = 226) | +0.0407 | **0.0015** |
| case vs control (n = 474) | +0.0185 | 0.069 |
| **unmedicated BP1 (n = 74) vs controls** | **−0.0123** | **0.352** |

The third row is the clean test: bipolar I patients who are *not* on lithium
are indistinguishable from controls on every compositional coordinate (all
p ≥ 0.138). Whatever moves the blood composition in this cohort is the drug,
not the diagnosis.

---

## 5. Two caveats that must travel with this result

### 5.1 Circularity

The cell fractions are estimated **from the same expression matrix** that
supplies the outcome. Conditioning on them therefore conditions on a summary of
the data itself, and for any gene that contributes to the deconvolution the
relationship is partly mechanical. The TRANS spike-in bounds how bad this can
be — a signal orthogonal to composition passes through untouched at 91–100% —
but it does not eliminate the concern, because a *real* global transcriptional
effect of lithium would itself perturb the estimated fractions and so be
partially absorbed. This is the "estimated mediator" problem the project
briefing flags for Arm 2, and the number above is a **lower bound on the
transcriptional component / upper bound on the compositional one**.

### 5.2 The mediator is measured differently in the two exposure groups

CIBERSORTx returns a per-sample fit quality. Regressed on lithium with the same
covariates:

| metric | non-users | users | t | p |
|:--|--:|--:|--:|--:|
| Correlation | 0.9064 | 0.9121 | 2.76 | 0.0063 |
| RMSE | 0.5146 | 0.4888 | −3.73 | **0.00025** |

Lithium users are fitted **better**, not worse. Library size is balanced
(p = 0.230), so this is not a depth artifact. It is differential measurement
error in the mediator and it should be reported; it does not have an obvious
direction of bias, but it means the fractions are not equally reliable across
the contrast.

---

## 6. Why the pre-specified balance b1 is the wrong single coordinate

This is a concrete recommendation for the main pipeline, not just an
observation about this experiment.

The sequential binary partition in `00_prep.R` puts granulocytes **and**
monocytes together in the numerator of b1. But lithium moves them in opposite
directions — neutrophils +4.64 percentage points, monocytes −1.25 — so lumping
them dilutes the very signal b1 is supposed to carry. The evidence:

- as a mediator, b1 is the weakest of the three myeloid-facing coordinates
  (p = 0.047, partial R² 0.018) while raw granulocyte proportion is the
  strongest (p = 0.0015, partial R² 0.046);
- as a DE adjustment, **logit(granulocyte) alone removes 83.4% of the canonical
  DEGs — more than all four ILR balances together (78.2%) — using one degree of
  freedom instead of four.**

Whether to revise the SBP is a call for the main project, but b1 should not be
quoted as "the" composition coordinate without this qualification.

**The practical implication is larger than the methodological one.** The single
most predictive composition variable is the one a haematology analyser prints
for free. A routine blood count, with no sequencing at all, accounts for more
of the lithium expression signature than a 22-type deconvolution does. That is
directly relevant to Arm 1's question of whether the transcriptome adds
anything over a blood count.

---

## 7. Files

| file | contents |
|:--|:--|
| `summary.csv` | the six-model comparison in §2 |
| `gene_level.csv` | per-gene canonical vs adjusted logFC, b1 loading, predicted indirect effect, status |
| `path_a_mediator_models.csv` | lithium → each compositional coordinate |
| `lithium_vs_celltypes.csv` | lithium → each of the 18 testable LM22 types |
| `composition_lithium_vs_illness.csv`, `composition_unmedicated_vs_control.csv` | §4.2 |
| `deconv_qc_vs_lithium.csv` | §5.2 |
| `collinearity.csv` | §3.1 |
| `permuted_ILR_null.csv`, `permutation_draws.csv` | §3.2 |
| `spikein_calibration.csv`, `spikein_calibration_curve.csv` | §3.3 (additive; TRANS arm valid, COMPOS arm mis-specified) |
| `spikein_varpreserving_curve.csv`, `spikein_varpreserving_anchors.csv` | §3.3 corrected |
| `estimator_sensitivity.csv` | voom / voomWQW / trend / QLF |
| `loss_by_b1_loading.csv` | DEG loss by strength of composition loading |
| `deg_by_lfc_floor.csv` | DEG counts at |logFC| floors 0 to 0.5 |
| `specificity_casecon.csv` | the same adjustment on the case/control contrast |
| `attenuation_depth.csv` | how far the lost genes actually fell |
| `marker_enrichment.csv`, `direction_check.csv`, `biotype_by_status.csv` | gene identity |
| `survivor_genes.csv`, `survivor_genes_named.csv` | candidate composition-independent lithium genes |
| `top_composition_driven.csv` | the clearest casualties |
| `logs/` | full stdout of every script |
