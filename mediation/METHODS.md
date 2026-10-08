# Mediation arm — methods

Everything this arm did, why, and what went wrong on the way.
Companion documents: `PRESPEC.md` (what was fixed in advance), `RESULTS.md` (what came out).

---

## 1. The question

Lithium is known to raise neutrophil proportions. Blood is a mixture, so changing the
mixture changes bulk gene expression arithmetically — without any cell behaving
differently. The differential-expression arm had already found lithium's whole-blood
signature to be large (1,426 genes) and largely compositional.

This arm puts a number on it:

> **How much of lithium's effect on whole-blood gene expression runs *through* the cell
> mixture, rather than directly?**

| | |
|:--|:--|
| exposure | lithium use, 0/1 |
| mediator | the ILR balances of the cell mixture, **as one vector** |
| outcome | gene expression, **one model per gene** |
| sample | the 226 bipolar I subjects — 152 users, 74 non-users |

Lithium use varies only within bipolar I, so nothing is gained by pooling and a great deal
is risked. Assessment group is constant in these 226 subjects and is therefore dropped;
`01_prep.R` asserts this rather than trusting it.

---

## 2. Two mediators, run side by side

The written proposal and the implemented pipeline disagreed about the cell-type tree, so
**both were run end to end** rather than one being chosen.

| | **lineage** (pipeline) | **proposal** (proposal text) |
|:--|:--|:--|
| parts | 5 lineages: gran, mono, T, NK, B | 6: neut, mono, CD4-naive, CD4-memory, NK, B |
| balances | 4 | 5 |
| cell mass used | 100 % | 92.4 % |
| samples needing zero replacement | **0** | 71 of 474 (15.0 %) |
| exact zeros among parts | 1 of 2,370 (0.04 %) | 72 of 2,844 (2.5 %) |

For reference, the raw 22 LM22 types are 46.1 % exact zeros, which is why neither tree uses
them directly: a log-ratio built on a part that is zero in 100 % of samples (NK activated)
or 98 % (B memory) measures nothing.

Both trees are **sequential binary partitions** — each balance splits a group formed by an
earlier split — and `01b_prep_proposal.R` verifies this rather than assuming it, checks the
contrast matrix is orthonormal to 2.2e-16, and confirms the transform is lossless (the
balances reconstruct the proportions to 2.8e-16).

Zeros are replaced **multiplicatively** (δ = 0.65 × the smallest observed value in that
part, other parts rescaled so the composition still sums to one). Additive replacement
would distort exactly the ratios we then take logs of.

---

## 3. The estimator

### 3.1 The models

Mediator models, k = 1..K, sharing one design:

```
M_k = b0_k + b1_k·X + b2_k'·C + e_k
```

Outcome model, per gene, **with exposure–mediator interaction**:

```
Y = t0 + t1·X + Σ_k t2_k·M_k + Σ_k t3_k·(X·M_k) + t4'·C + e
```

Covariates `C`: age, sex, tobacco, RIN, plate, seqPC1–3 — the same set as the DE arm, so
the two arms adjust identically.

`M` and `C` are centred at their sample means. This makes the mediator intercept `b0_k`
*be* the quantity the direct-effect formula needs, and removes the near-collinearity
between `X` and `X·M` that raw ILR values (means 1.88, 0.24, 1.76, 1.37) would create.
Centring provably changes nothing — the design spans the same column space — and
`00_validate_med.R` checks that numerically rather than asserting it.

### 3.2 The estimands

```
NDE = t1 + Σ_k t3_k·b0_k                 natural direct effect
NIE = Σ_k (t2_k + t3_k)·b1_k             natural indirect effect   ← PRIMARY
TE  = NDE + NIE
PM  = NIE / TE                           proportion mediated
```

These are Pearl's (2001) natural effects, after Robins & Greenland (1992), in the
regression form of VanderWeele & Vansteelandt (2009, 2010), extended to a **vector**
mediator by VanderWeele & Vansteelandt (2014) — the reference already in the proposal.

Also reported, the **four-way decomposition** (VanderWeele 2014):

```
TE = CDE + INT_ref + INT_med + PIE
```

the parts of the total effect needing neither mediation nor interaction, interaction only,
both, and mediation only.

**A limit worth stating.** Because the mediator is a vector, this is the **joint** indirect
effect through the whole cell mixture. Splitting it into path-specific effects through
*individual* balances would need an assumed causal ordering among the balances, which we do
not have. Per-balance terms are reported as a **descriptive decomposition**, never as four
separate causal effects.

### 3.3 Why not product-of-coefficients, which the proposal's flowchart specifies

Because the feedback on the pre-proposal said not to:

> "product of coefficients assumes no lithium-by-cell-mix interaction. **Test that
> interaction before reporting a proportion mediated.**"

The estimator above *contains* product-of-coefficients as the special case `t3 = 0`, so
nothing is abandoned. We additionally run a joint test of `H0: t3 = 0` per gene (moderated
F on K df, BH-corrected) and report **both estimators side by side**, so the cost of the
original assumption is a number rather than an argument. The proposal was updated to match
and now says so.

### 3.4 Standard errors

Delta method, with an **exactly block-diagonal** joint covariance:

```
Var(g) = a'·V_theta·a + b'·V_beta·b        (no cross term)
```

This is not a convenience approximation. The outcome model's estimating equation is driven
by `e_Y = Y − E[Y|X,M,C]`, the mediator models' by `e_M = M − E[M|X,C]`. Since
`E[e_Y|X,M,C] = 0` and `e_M` is a function of `(M,X,C)`, iterated expectation gives
`E[e_Y·e_M·f(X,M,C)] = 0`. The two scores are uncorrelated.

The K mediator models share one design, so this is a seemingly-unrelated regression in
which equation-by-equation OLS is fully efficient and

```
V_beta = kronecker( (D'D)^-1 [1:2,1:2] , Sigma )
```

carries the **full** across-mediator correlation. Nothing is assumed independent that is
not.

Per-gene residual variances are shrunk with limma's empirical Bayes (`squeezeVar`, the
engine inside `eBayes`), for consistency with the DE arm, which used limma-voom with eBayes
throughout.

---

## 4. The identity that looked like a bug, and was not

In a linear model `total = direct + indirect` is an exact algebraic identity, so `TE`
should equal the lithium coefficient from the plain model `Y ~ lithium + C` to machine
precision. **It did not** — the gap reached 4.3e-02.

The identity requires both models to use the *same projection*. Ours deliberately do not:

* **outcome models are weighted** (voom precision weights), for comparability with the DE arm;
* **mediator models are unweighted**, because voom weights describe count-level noise in
  *gene expression*. They say nothing about how precisely a subject's cell fractions were
  estimated, and the lithium→cell-mix effect is a single gene-independent quantity.
  Weighting it per gene would claim lithium shifts the myeloid/lymphoid balance by 0.1558
  when you look at one gene and 0.1949 when you look at another — not a coherent statement
  about blood.

`02c_identity.R` proves this is the whole explanation: refit the mediator model with each
gene's weights so the projections match, and the gap returns to **2.5e-16**. The unweighted
choice is kept, and **S12** confirms the conclusions do not depend on it (per-gene NIE
correlation 0.9984; median PM 0.396 weighted vs 0.400 unweighted).

---

## 5. Validation — `00_validate_med.R`, 25 checks, 0 failures

Nothing checks the estimator against my own algebra. Every test checks it against something
independent.

| # | Check | Result |
|:--|:--|:--|
| 1–3 | **Brute-force Monte Carlo of the counterfactuals.** Potential outcomes `Y(1,M(0))`, `Y(0,M(0))`, `Y(1,M(1))` generated directly from known truth over 4 M draws; natural effects taken as plain averages, no formula reused | matches the closed form to 1.4e-4 |
| 4–5 | estimator recovers known truth at n = 200,000 | NDE 0.7230 vs 0.7250; NIE 0.6247 vs 0.6200 |
| 6–7 | single-mediator case vs an independent hand computation with `lm()` | agrees to 1e-9 |
| 8–9 | four-way decomposition sums to TE; `CDE+INT_ref = NDE`, `INT_med+PIE = NIE` | exact (diff 0) |
| 10 | `NIE − product-of-coefficients = INT_med` | exact to 1e-12 |
| 11 | with no true interaction, `INT_med` is null | z = 0.54 |
| 12–13 | centring invariance | identical to 1e-8 |
| 14–15 | design full rank; voom weights unchanged by centring | rank 21/21; max diff 5.3e-12 |
| 16–18 | **delta SE vs a 4,000-replicate bootstrap** | within 1.9–6.3 % |
| 19–20 | **interval coverage over 500 simulated datasets** | 93.0 %, 94.4 % vs nominal 95 % |
| 21–23 | mediator covariance equals the kronecker form; SE matches `lm()`; cross-mediator covariance carried | exact |
| 24–25 | degenerate inputs (no mediation → NIE null; no direct effect → PM ≈ 1) | z = 1.27; PM = 0.935 |

**Finding that changed the plan:** the delta method is mildly **anticonservative** —
about 5 % low in simulation, about 10 % low against a direct bootstrap on the real data
(`02b_diagnose.R`). That is why the bootstrap is not decoration here.

---

## 6. Resampling

| | |
|:--|:--|
| **Bootstrap** | B = 1000, resampling **subjects** with replacement |
| **Permutation** | B = 1000, **lithium labels shuffled** across subjects |
| refitted inside every replicate | TMM normalisation, **voom weights**, mediator models, centring |
| cores | 4 (six thrashed a 7.7 GB machine earlier in this project) |

Refitting voom inside every replicate costs about 4.5 s of the ~9.5 s per replicate. It is
paid because the mean–variance trend is *estimated*; holding the weights fixed would treat
an estimated quantity as known and understate the variance.

The resampling unit is the subject, who carries exposure, balances, covariates and whole
expression profile together. The permutation breaks lithium's link to **both** mediator and
outcome at once, so the null is "lithium does nothing", under which NDE, NIE and TE are all
zero. Permutation p-values are **one-sided**: signal can only add discoveries.

Index sets and shuffles are drawn once in the master process from seed 481, so the run is
reproducible regardless of how work lands on workers. `med_replicate` run on the identity
index reproduces the primary estimates to **5.4e-15**.

---

## 7. Sensitivity and secondary analyses

All declared in `PRESPEC.md` **before** any estimate existed.

| # | Analysis | Script |
|:--|:--|:--|
| S1 | 350 LM22 marker genes analysed separately — the circularity measurement | `05` |
| S2 | second, independent deconvolution (Krebs' deposited CIBERSORT) | `05` |
| S3 | regression calibration for the estimated mediator (Valeri, Lin & VanderWeele 2014) | `06` |
| S4 | E-value (VanderWeele & Ding 2017) | `06` |
| S5 | all 20 tobacco imputations, pooled by Rubin's rules | `05` |
| S6 | zero-replacement δ swept over {0.1, 0.3, 0.5, 0.65, 0.8} | `05` |
| S7 | correction panel: BH, BY, Bonferroni, Holm, Storey q | `02` |
| S8 | product-of-coefficients reported alongside | `02` |
| S9 | unadjusted and age-only-adjusted | `05` |
| S10 | power / minimum detectable NIE | `06` |
| S11 | the proposal's own 6-part tree, as a full parallel primary | `01b`, `02` |
| S12 | unweighted outcome models, where the identity is exact | `05` |
| T1 | **exploratory** — age as exposure, per 10 years | `07` |
| T2 | **exploratory** — moderated mediation, with its power | `07` |

### Notes on individual analyses

**S3 regression calibration.** Classical error in a *regressor* attenuates its coefficient,
so naive `t2`/`t3` are too small and the mediated effect is understated. The multivariate
correction is `t2_true = (Σ_W^-1 Σ_M)^-1 · t2_naive`, with `Σ_U` built from the per-balance
ICCs (b1 0.720, b2 0.585, b3 0.642, b4 0.516) measured in `00_reliability.R`. `b1` is *not*
corrected: error in a mediator model's *outcome* inflates its standard error but leaves the
slope unbiased. **The ICCs themselves overstate reliability**, because both deconvolutions
share the LM22 signature and the same algorithm family — so the correction is too small and
the corrected figure is *also* a lower bound. Every caveat in this arm points the same way.

**S5 Rubin pooling is a no-op here, and that is a finding.** The 30 subjects with imputed
tobacco status are all controls or bipolar II; **0 of the 226 bipolar I subjects** in this
contrast had tobacco imputed. The 20 imputations are therefore identical on this sample
(fraction of missing information 0.000). Reported rather than quietly dropped.

**S4 E-value — one thing is incomplete.** The proposal specifies benchmarking against the
obesity→neutrophil association (Furuncuoğlu et al. 2016), because bodyweight is not recorded
in GSE124326. **That paper is not in the project library, so its numbers were not used
rather than guessed.** Instead the E-values are benchmarked against the measured covariates
in this very dataset, on the same scale. Extracting the Furuncuoğlu figure is an open action
item for the write-up.

---

## 8. Assumptions, stated plainly

Natural effects are identified under (Pearl 2001; Imai, Keele & Yamamoto 2010):

1. no unmeasured exposure–outcome confounding given C;
2. no unmeasured mediator–outcome confounding given C and X;
3. no unmeasured exposure–mediator confounding given C;
4. no mediator–outcome confounder itself affected by the exposure.

**(2) and (4) are not guaranteed and cannot be made so.** Illness severity, mood state and
other psychotropics plausibly affect lithium prescription, cell mix and expression together,
and **none is in this dataset**. There is no lithium dose. BMI is not recorded — checked, it
does not exist in GSE124326 — so it is bounded rather than adjusted. This is why S4 is not
optional, and why every estimate here is **associational under a causal interpretation**.

A further assumption specific to this design: the mediator is **estimated from the same
expression matrix** that supplies the outcome. S1 measures the damage; the 547 LM22
signature genes (350 after filtering) are held out of the primary set entirely.

---

## 9. Incidents — mistakes made in this arm and how they were caught

Recorded because a methods section that reports only successes is not a methods section.

1. **A validation test was wrong, not the code.** The check "with no true interaction, NIE
   and product-of-coefficients agree closely" failed at tolerance 0.02. The gap (0.0217) was
   about one standard error: with true `t3 = 0` the *estimated* `t3` is still noise, so the
   two estimators must differ by sampling error. Replaced with the correct pair of tests —
   an exact algebraic identity (`NIE − prod = INT_med`) and a z-test on `INT_med`.

2. **The total-effect identity appeared to fail.** Traced to the deliberate
   weighted-outcome/unweighted-mediator mismatch, proved in `02c_identity.R`, and shown
   immaterial by S12. See section 4.

3. **A heuristic was reported as a hard bound.** `02b_diagnose.R` originally flagged the
   proposal variant for "violating" a root-sum-square ceiling on the NIE z-statistic. That
   bound is only valid for mutually orthogonal balances; ours correlate up to |r| = 0.86.
   The verdict text was corrected and the **real** ceiling is now computed properly in S10b
   by zeroing the outcome-model variance.

4. **An apparent contradiction between the two mediator trees was real, and is the
   headline caveat.** Per-gene NIE estimates correlate at 0.980 across trees, yet BH gave 0
   significant genes for one and 3,647 for the other. Not a bug: BH is a step-up procedure,
   so a diffuse shift in the whole p-value distribution flips it between almost nothing and
   thousands. The discovery *count* is therefore not a stable quantity here, which is
   precisely what the permutation null was built to calibrate.

5. **An independent adversarial review found four real defects, after the 25-check
   validation suite had already passed.** Five reviewers with distinct lenses (estimand
   algebra, variance derivation, code, statistical interpretation, proposal compliance)
   examined the arm; thirteen candidate findings were put to adversarial verifiers
   instructed to refute them, and six survived. Four were genuine and are fixed:

   * **The "product-of-coefficients" comparison was tautological.** `NIE_prod` was coded as
     `T2 %*% b1` using `t2` from the model that *contains* the interaction columns. That is
     byte-identical to `PIE`, so `NIE - NIE_prod` was `INTmed` by construction and the
     reported "cost of the no-interaction assumption" was a within-model decomposition term
     wearing a disguise. Zeroing `t3` in the *formula* is not the same fit as removing the
     interaction columns from the *design* — for a binary exposure the no-interaction slope
     is a within-group weighted average, `t2_ni = t2 + w1*t3`. `prod_of_coef()` now refits
     properly. **The published figures changed: 25.3 % to 8.8 %, and sign disagreement 6.6 %
     to 2.30 %.** The correction *strengthens* the conclusion that the proposal's original
     estimator was defensible. Two new validation checks now guard it.
   * **The regression calibration was inadmissible and nobody stopped.** The reliabilities
     imply more measurement error than the observed balance covariance can carry: the
     implied true-score covariance had two negative eigenvalues. The script detected this,
     printed it, applied a fallback that was *also* not positive definite, printed that too,
     and carried on to produce 0.429. It now shrinks the error covariance to the largest
     admissible multiple (c = 0.56), asserts positive definiteness, and — because the matrix
     correction turns out to be ill-conditioned at that boundary (condition number 92,
     giving an artefactual 0.877) — reports the **stable per-balance correction (0.587)**
     plus the full sweep instead of a point on a cliff edge.
   * **The calibration let the total effect drift 18.5 %.** `t2` and `t3` were corrected but
     `t1` was not, so `TE = NDE + NIE` moved. TE is the coefficient of a model containing no
     mediator, so measurement error in the mediator cannot change it; a correction may only
     redistribute effect between the two paths. The direct effect is now obtained by
     subtraction and the identity is asserted (`max |TEc - TE| = 0`).
   * **The bootstrap never produced an interval for the headline quantity.** It saved
     per-gene intervals but summarised the proportion mediated as an unrestricted *mean*
     over all 12,018 genes — a rule `PRESPEC.md` §1.4 explicitly forbids — and no bootstrap
     number reached the report at all. It now computes the pre-declared estimand per
     replicate (median PM over the genes with a significant total effect, gene set fixed
     from the primary analysis) and reports its percentile interval.

6. **One confirmed finding was rejected after testing it, and the test mattered.** The
   review argued that `SE(INTref)` is inflated because the centred mediator model forces
   `b0 = -xbar * b1` identically, so `Var(b0)` should be `xbar^2 Var(b1)` rather than the
   textbook `(D'D)^-1[1,1] * Sigma` — the two differ by exactly `1/n`, which is true, and we
   verified the identity holds to 1.9e-16. But the claim is only correct *conditional on the
   sample centring*. Our inference is unconditional: the primary intervals come from
   resampling subjects, under which the sample means are themselves random. Arbitrated
   against a 500-replicate subject bootstrap, the textbook form is **closer** to the truth
   and the proposed fix makes it **worse**:

   | | current (textbook) | proposed fix |
   |:--|--:|--:|
   | SE(NDE) / bootstrap SE | **0.960** | 0.955 |
   | SE(INTref) / bootstrap SE | **0.805** | 0.760 |
   | SE(TE) / bootstrap SE | **0.989** | 0.984 |

   The fix was **not applied**. What the exchange did establish, and what is reported, is
   that `SE(INTref)` is the least trustworthy standard error in the arm (about 20 % below
   the bootstrap), which is one more reason the bootstrap is the primary interval.

7. **Attribution error on my part.** I initially cited the 39 PDF comments as the team's own
   review notes after a correction, then reverted to Prof. Ay once Ali confirmed they are
   his. The annotation author field reads `ahmetay` throughout.

---

## 10. Files

| Path | Holds |
|:--|:--|
| `PRESPEC.md` | the protocol, fixed before any estimate existed |
| `METHODS.md` | this document |
| `RESULTS.md` | the findings |
| `scripts/med.R` | the estimator — one implementation, used everywhere |
| `scripts/00_validate_med.R` | 25 correctness checks |
| `scripts/00_reliability.R` | per-balance ICCs from two deconvolutions |
| `scripts/00b_ilr_trace.R` | the ILR construction, traced step by step |
| `scripts/01_prep.R` / `01b_prep_proposal.R` | the two mediator datasets |
| `scripts/02_effects.R` | the primary result, both trees |
| `scripts/02b_diagnose.R` / `02c_identity.R` | why the two trees differ; the identity proof |
| `scripts/02d_vs_de.R` | external validation against the differential-expression arm |
| `scripts/03_bootstrap.R` / `04_permutation.R` | B = 1000 each, both trees |
| `scripts/04b_perm_pm.R` | permutation null for the headline quantity, on a fixed gene set |
| `scripts/05_sensitivity.R` / `06_calibration.R` / `07_secondary.R` | S1–S12, T1–T2 |
| `scripts/08_report.R` | assembles `RESULTS.md` |
| `results/primary/` | per-gene effect tables, both trees |
| `results/boot/`, `results/perm/` | resampling output |
| `results/sens/` | every sensitivity analysis |
| `logs/` | run logs |
