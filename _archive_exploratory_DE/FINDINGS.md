# Findings

> **Verification status.** All agent-sourced numbers were independently
> re-computed by a different route (`scripts/VERIFY_agent_claims.R`,
> `runs/VERIFY/verdict.csv`). **15 of 16 claims reproduce exactly or within
> tolerance.** The one apparent failure — the compositional spike-in survival
> rate — was traced to a flaw in the *verification*, not the original: my spike
> was absorbed exactly by the adjustment (1.5e-15), and 100% of the apparent
> survivors were genes already differentially expressed before any spike was
> added. Corrected, survival is 0.000, matching. The only genuine discrepancy is
> split-half reproducibility (agent 0.32, mine 0.615), which is a difference in
> the denominator definition, not in the data; sign agreement is 0.97–0.99 both ways.

**Exploratory.** Nothing here is pre-specified. Thresholds were varied on
purpose, several analyses were chosen after seeing earlier ones, and no p-value
is corrected for the number of analyses run. Anything promoted to the paper has
to be re-declared in `config.yaml` and re-run on that basis.

Canonical settings throughout unless stated: filter > 10 counts in ≥ 90% of
samples, TMM, voom, FDR 0.05.

| Contrast | n | Genes | DEG | π₀ |
|:--|--:|--:|--:|--:|
| `lithium` (BP1, 74 vs 152) | 226 | 12,173 | 1,382 | 0.516 |
| `casecon` (BD vs control) | 474 | 12,368 | 1,816 | 0.546 |
| `lithium_krebs` (240 cases) | 240 | 12,207 | 1,459 | 0.503 |

---

## 1. The bipolar signature is a lithium signature — three independent routes

| Route | Case-control DEGs |
|:--|--:|
| Ignore lithium entirely | **1,816** |
| **Adjust** for lithium (Krebs' published spec, n=444) | **6** |
| **Stratify** to bipolar I patients not on lithium (74 vs 234) | **4** |

Two methods that share almost nothing — covariate adjustment and sample
stratification — converge on the same answer. The published spec reproduces
Krebs' own **6** DEGs exactly: 5 of their 6 genes recovered (C6orf163, COG4,
PVT1, DOCK3, BBS9, RP11-53B2.2), 100% sign concordance.

By contrast the on-lithium comparison gives **2,603**. The published
case-control logFC vector correlates **0.959** with the *on-lithium* contrast
and **0.130** with the *off-lithium* one.

### Not a power artefact

Subsampling on-lithium cases to exactly 74 against the same 234 controls, 40 draws:

| | n cases | DEG |
|:--|--:|:--|
| off lithium (observed) | 74 | **4** |
| on lithium, subsampled to 74 | 74 | median **1,472**, range 505–2,784 |
| **control vs control (null calibration)** | 74 | **median 0, max 0** |

0 of 40 draws came within two orders of magnitude of 4. The control-vs-control
null returns exactly zero, so the pipeline is calibrated.

### Not an estimator artefact

| Estimator | off lithium | on lithium | ratio |
|:--|--:|--:|--:|
| voom | 4 | 2,603 | 651× |
| limma-trend | 0 | 2,625 | >2000× |
| edgeR-QLF | 11 | 2,730 | 248× |

### Confirmed on the original authors' gene list

Krebs' published 3,031-gene lithium DEG list used as a probe set:

| | median p over their genes | background | DEG |
|:--|--:|--:|--:|
| off lithium vs control | **0.347** | 0.327 | 1 |
| on lithium vs control | **0.0018** | 0.160 | 2,477 |

Off lithium, their gene set sits at background.

> **The caveat that governs this entire section.** "Off lithium" is **not**
> "unmedicated." The dataset records lithium use and nothing else — no
> antipsychotics, anticonvulsants, antidepressants, dose or duration. Those 74
> subjects may well be on other psychotropics. The defensible claim is that the
> signature tracks **lithium specifically**, not that unmedicated bipolar blood
> is indistinguishable from healthy blood.

---

## 2. Roughly three quarters to nine tenths of the lithium signature is compositional

| | DEG before | after | removed | mean \|logFC\| shrinkage |
|:--|--:|--:|--:|--:|
| **Ours** (4 ILR balances) | 1,459 | 328 | **77.5%** | 34% |
| **Krebs** (their own correction, published) | 3,031 | 233 | **92.3%** | 46% |

Two independent groups, different adjustments, same dataset, same direction.

### The adjustment does what it claims — spike-in validated

A synthetic-signal calibration (agent-run, `runs/cellcomp/`) spiked 1,200 genes
with either a transcriptional effect or a compositional one and asked whether
the balance adjustment could tell them apart:

| Spiked effect | Detected | Survives adjustment | logFC shrinkage |
|:--|--:|--:|--:|
| **Transcriptional** | 879/1,200 | **91%** | **1.3%** (none) |
| **Compositional** | 2/1,200 | **0%** | **83%** |

So the 77.5% figure is not the adjustment over-fitting. It specifically removes
compositional signal and specifically preserves transcriptional signal.

### The closed-form version, and where it breaks

Under pure mediation, `β_lithium = γ_balance × δ` for every gene.

- **δ (path a)** = +0.156 (SE 0.078, p = 0.047, partial R² = 0.017) — real but weak
- β vs γ across 12,173 genes: Pearson **+0.695**; among the 1,382 DEGs, **+0.835**
- Fitted slope **+0.309** [0.303, 0.315] vs δ = +0.156 — **δ outside the interval**
- Composition accounts for **36.5%** of the variance in the lithium effect

The slope exceeding δ looks like proof of a direct effect, but it is not clean:
γ is a coefficient on an *estimated* balance, so measurement error attenuates
it, which inflates the β/γ slope by exactly this kind of factor. δ is unbiased
(the balance is the outcome there). **These numbers cannot separate a genuine
direct effect from attenuation in γ** — which is what the pre-proposal's
regression-calibration step exists to fix.

---

## 3. Part of the "direct effect" is composition the deconvolution cannot see

408 genes keep a lithium association after adjustment. The two strongest are
**IL5RA** (β = +0.73) and **CLC/galectin-10** (+0.61) — both eosinophil/basophil
genes. LM22 puts eosinophils at 0.1% of leukocytes against a 1–6% clinical
reference and **has no basophil type at all**.

Marker panels fixed from the literature before looking:

| Panel | Composition explains | % "direct" |
|:--|--:|--:|
| Neutrophil + T cell — LM22 **sees** these | **48%** | 5% |
| Eosinophil — LM22 is **blind** | **6%** | **50%** |

Background direct rate 3.4%; hypergeometric **p = 6.9 × 10⁻⁴**. The basophil
panel had *zero* genes pass the expression filter. Gene by gene: IL5RA, CLC,
CYSLTR2 (+0.43), IL1RL1/ST2 (+0.38) — a coherent type-2/eosinophil axis, all
rising with lithium, all essentially unexplained by the balances.

### Abundance, not induction

| Panel | mean pairwise r | r with LM22's own estimate |
|:--|--:|--:|
| **Eosinophil** | **+0.719** | +0.123 (nonzero in 14% of samples) |
| Neutrophil (abundance benchmark) | +0.486 | +0.403 |
| T cell | +0.467 | **−0.214** |

The eosinophil genes co-vary *more tightly than the neutrophil benchmark* — one
latent cell-number factor, not independent induction. So the mediated
proportion is a lower bound **for a second, separate reason**: the mediator is
blind to some of the composition it should be capturing.

Lithium's effect on the eosinophil panel score is d = +0.26, p = 0.063 on four
genes — suggestive, not significant. A measured eosinophil count would settle it.

> **Incidental red flag.** LM22's T-cell fraction *anti-correlates* (−0.21) with
> canonical T-cell markers (CD3E, CD3D, CD2, IL7R, LCK, ZAP70).

---

## 4. Cell composition dominates this transcriptome; lithium is a minor axis

Each variable used as its own exposure, same 226 subjects, same covariate set:

| Exposure | DEG | π₀ |
|:--|--:|--:|
| **myeloid/lymphoid balance (b1)** | **9,165** | **0.153** |
| sequencing metric PC3 | 6,945 | 0.276 |
| sequencing plate | 5,842 | 0.327 |
| RIN | 3,174 | 0.488 |
| **lithium** | **1,382** | 0.516 |
| sex | 906 | 0.703 |
| age | 106 | 0.802 |
| tobacco | 88 | 0.884 |

Cell composition is associated with **75% of the transcriptome** and carries
**6.6× more signal than lithium**. Three technical variables also beat lithium.

Consistent with this, residual PCA after the full covariate model:

| Residual PC | var. explained | r with granulocyte % | r with lithium |
|:--|--:|--:|--:|
| **PC1** | **20.2%** | **0.688** | 0.001 |
| PC2 | 9.6% | 0.417 | 0.003 |

The single largest unmodelled axis of variation is cell composition, and it is
orthogonal to lithium. Full-matrix PC1 (20%) and PC2 (11%) correlate 0.55 and
0.58 with granulocyte proportion.

> **Do not read the surrogate-variable DEG counts as signal.** Adding residual
> PCs raises the lithium DEG count 1,382 → 5,286, but mean |logFC| is unchanged
> (0.0799 → 0.0799). Residual PCs are orthogonal to the design by construction,
> so they shrink the error term without moving any effect estimate. That is an
> artefact of the naive construction, not added power — and a trap worth
> documenting, since a real `sva` run guards against it.

---

## 5. Normalisation choice changes the answer threefold — and it is not innocent

| Normalisation | DEG | π₀ | Jaccard vs TMM |
|:--|--:|--:|--:|
| TMM (baseline) | 1,382 | 0.516 | 1.00 |
| RLE | 1,407 | 0.513 | 0.96 |
| TMMwsp | 1,348 | 0.531 | 0.92 |
| upper quartile | 1,672 | 0.480 | 0.71 |
| **none (plain CPM)** | **4,013** | 0.340 | **0.21** |

TMM, RLE and TMMwsp agree. Skipping normalisation inflates the count **2.9×**
and shares a fifth of its gene list.

The reason matters here. TMM's own normalisation factors are **explained by cell
composition**: r² = 0.237 against the four ILR balances, r² = 0.227 against
granulocyte proportion, r = −0.476 with granulocytes. TMM corrects for shifts in
library composition, and a neutrophil shift *is* such a shift — so the
normalisation step is partly removing the effect under study. This needs stating
explicitly in the write-up rather than left as a generic methods sentence.

---

## 6. Effect sizes are small, which is itself consistent with a compositional cause

| FDR | \|logFC\| ≥ 0 | ≥ 0.2 | ≥ 0.5 | ≥ 1.0 |
|--:|--:|--:|--:|--:|
| 0.001 | 156 | 137 | 19 | **0** |
| 0.01 | 495 | 320 | 23 | **0** |
| 0.05 | 1,382 | — | — | — |

Nothing reaches a two-fold change. Many genes moving slightly and coherently is
what a mixture shift produces; a few genes moving hard is what transcriptional
regulation produces.

---

## 7. The signal is real, but DEG *counts* are not stable quantities

**Permutation null** (lithium labels permuted, covariates left attached):
mean and max **0 DEGs** at FDR 0.05; null π₀ = 0.971 against an observed 0.516.
Enrichment 312× at FDR 0.001. The signal is not an artefact of the pipeline.

**But the count is fragile.** 100 random subsamples at n=200 — dropping just 26
of 226 subjects — give a median of 1,011 DEGs with **range 251–2,894**. An 11×
spread from sampling alone.

**Split-half reproducibility** is 0.32: only about a third of DEGs recur in a
held-out half, though sign agreement is 0.967. Directions replicate; individual
gene calls do not.

Practical consequence: report effect-size distributions, π₀ and enrichment, not
headline DEG counts. The 4-vs-2,603 result survives because its gap (126× below
the minimum of 40 matched draws) is far outside this sampling variability.

---

## 8. Robustness checks that came back clean

**Exposure × mediator interaction gate — PASSES.** No lithium × ILR, × age or
× sex interaction survives BH at FDR 0.05, so the pre-specified
product-of-coefficients decomposition remains valid. Worth noting honestly: the
lithium × b1 term shows mild inflation (10.2% of genes at nominal p < 0.05
against 5% expected, π₀ = 0.80). Not significant, but not a clean null either —
it may become significant with more power.

**Tobacco imputation barely matters.** Across all 20 completed columns: DEG
median 1,764, range 1,654–1,826, SD 38.5. 1,610 genes DEG in all 20.
Rubin-pooled 1,787 vs naive 1,792 — SEs understated by a median of **0.08%**
if imputation uncertainty is ignored. Doing it properly was cheap and correct;
it changes nothing. That is a reassuring negative result, and consistent with
tobacco carrying only 88 DEGs on its own.

**p-value histograms are well shaped** for all three contrasts: a ~9–10× spike
in the first bin over the mid-range, roughly flat elsewhere, far-tail ratio
0.79–0.85. All three sit slightly below 1, similarly so; that mild dip is most
likely inter-gene correlation rather than misspecification. No contrast shows a
U shape or a spike at 1.

**Model choice is not driving anything.** Across voom, voomWQW, limma-trend,
edgeR-QLF and edgeR-LRT, DEG counts stay in the same range and test statistics
correlate highly.

---

## 9. Things that matter for the project's decisions

**The depth/RIN decision does change the answer.** Screening to depth ≥ 3M and
RIN ≥ 5 drops the lithium contrast from 1,382 to 1,117 DEGs with a Jaccard of
only **0.65** — a third of the gene list changes. The decision to follow Krebs
and adjust rather than exclude is defensible, but it is not free and should be
reported as a sensitivity analysis.

**Krebs' published lithium count is internally inconsistent.** Their Results
text says 976 DEGs, their MAGMA methods text says 897, and their own deposited
pre-correction sheet has 3,031 at FDR 0.05. Our reproduction of their stated
spec gives **1,023** — within 5% of the 976 in the text. The 3,031 in their
supplementary file is not the number they report.

**Assessment group is inestimable in any cases-only design.** All 240 cases are
group A (confirmed in both the 474 and 444 cohorts). The published lithium model
lists ascertainment group as a covariate, but it cannot have been estimated.
Our engine detects and drops it; anyone reproducing the published model needs to
know this.

---

## Caveats applying to everything above

1. **"Off lithium" ≠ "unmedicated"** — no other medication is recorded.
2. **The mediator is derived from the outcome matrix**, so composition
   adjustment is partly self-referential. §3 quantifies one consequence.
3. **The mediator is blind to eosinophils and basophils**, so the mediated
   proportion is a lower bound for that reason too.
4. **Many thresholds were tried on purpose** and nothing is corrected for that.
5. **DEG counts are unstable** at this sample size (§7). Prefer effect sizes,
   π₀ and enrichment.
6. **This is one cohort.** Every result is within-GSE124326 and unreplicated
   externally.
