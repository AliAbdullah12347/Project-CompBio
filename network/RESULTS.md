# Arm 3 — Results

**Updated:** 2026-10-09 session 2
**Runs completed:** 0 of 48 queued (base-001 RUNNING)

## The question

Do gene co-expression networks differ across healthy / bipolar-off-lithium /
bipolar-on-lithium, and does any difference survive removing cell composition?

Exploratory. A well-documented negative or inconclusive result is a complete
deliverable.

---

## Scale-free topology

*Confirmed before any network was built.*
Scale-free fit fails at all powers 1:20 in the control group (n=234).
Max signed-R² = 0.955 at power 2 (positive slope — not scale-free at any power).
WGCNA FAQ for signed networks, n>40: use power 12 when no power satisfies criterion.
The failure of scale-free topology in this data is itself a result.

Source: `runs/base-001/sft_r2_curve.csv`

| group | n | powerEstimate | max signed-R² | power used |
|:------|--:|:-------------|:-------------|:----------|
| control | 234 | 1 (function failing) | 0.955 | **12** |
| bp_nolith | 74 | TBD (pow-009) | TBD | **12** |
| bp_lith (74 subsample) | 74 | TBD (pow-009) | TBD | **12** |

---

## Reference network — base-001

| parameter | value |
|:----------|:------|
| genes | 12,368 |
| samples | 234 controls |
| power | 12 (signed, bicor) |
| minModuleSize | 30 |
| mergeCutHeight | 0.25 |
| n_modules | **building** |

---

## Module preservation — base-001 (raw expression)

Status: **RUNNING** — reference TOM building.

### bp_nolith (74 samples, 500 permutations)

*(pending)*

### bp_lith (20 draws of 74 from 152, 50 perms/draw)

*(pending)*

---

## Module preservation — base-002-ceil (split-half ceiling)

*(pending — script written, awaiting base-001 completion)*

160-control reference vs 74 held-out controls. Gives the n=74 Zsummary ceiling:
without this, we cannot distinguish low bp_lith Zsummary from ordinary small-n noise.

---

## KEY QUESTION: composition residualised (inp-022)

*(pending — script written)*

Reference and test expression regressed for ILR b1-b4 (cell composition proxies).
If base-001 finds reduced Zsummary in bp_lith and it persists here: disruption is
not explained by composition. If Zsummary recovers: composition was the driver.

---

## KEY QUESTION: technical+composition residualised (inp-023)

*(pending — script written)*

Most conservative: age + sex + rin + ILR b1-b4. Any remaining preservation after
doubly residualised expression cannot be attributed to technical or compositional
confounding.

---

## Circularity check — minus LM22 genes (gene-027)

*(pending — script written)*

547 LM22 signature genes excluded before network construction. CIBERSORTx used
LM22 to estimate fractions. If preservation results are similar with and without
these genes, the findings do not depend on the fraction-estimation method.

---

## Technical residualised (inp-021)

*(pending — script written)*

Full technical set: age + sex + rin + plate + seqpc1-3. Contrasts with inp-023
(age + sex + rin + ILR, no plate/seqPCs).

---

## Summary table (to be populated)

| id | family | comparison | Zsummary bp_nolith | Zsummary bp_lith (mean±sd) | verdict |
|:---|:-------|:-----------|:-------------------|:---------------------------|:--------|
| base-001 | baseline | raw | — | — | running |
| base-002-ceil | baseline | ctrl ceiling | — | — | pending |
| inp-022 | input | ILR residualised | — | — | pending |
| inp-023 | input | tech+ILR residualised | — | — | pending |
| inp-021 | input | full-technical residualised | — | — | pending |
| gene-027 | geneset | minus LM22 | — | — | pending |

---

## Results file index

| file | experiment | content |
|:-----|:-----------|:--------|
| `runs/base-001/sft_r2_curve.csv` | base-001 | Scale-free R² by power, n=234 controls |
| `cache/ref_control_all_bicor_signed_p12.rds` | base-001 | Reference modules (building) |
| `runs/base-001/nolith_preservation.rds` | base-001 | Pending |
| `runs/base-001/results_nolith.csv` | base-001 | Pending |
| `runs/base-001/results_lith_alldraws.csv` | base-001 | Pending |
| `runs/base-001/results_lith_summary.csv` | base-001 | Pending |
| `runs/base-001/results_combined.csv` | base-001 | Pending |
