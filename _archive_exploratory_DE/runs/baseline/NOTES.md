# baseline — reference fits, and a replication of Krebs et al. (2020)

Everything here was computed from the data. No number is carried forward from the
paper or from the project notes without being recomputed; where a published figure
is quoted it is labelled as such and sits next to the figure we obtained.

Scripts, in the order they run:

| script | what it does |
|:--|:--|
| `scripts/ref_krebs.R` | loads the deposited Krebs DE tables and the comparison functions |
| `scripts/baseline.R` | canonical fits, replication, specification search, design diagnostics |
| `scripts/baseline_gap.R` | case/control decomposition, one-subject sensitivity, estimator variants |
| `scripts/baseline_confound.R` | how much of the case/control signature is lithium |
| `scripts/baseline_bp2.R` | why 13 bipolar-II subjects change the answer by 2.5× |

Fits are cached outside the project (scratchpad) and keyed by specification name, so
a rerun is cheap and no stale object can be analysed by accident. Nothing outside
`experimentation/` was written to. `data/gene_symbols.csv` was derived from the
GENCODE v19 GTF and is used only to attach readable symbols, never to join.

---

## 1. Canonical fits

voom + TMM, filter `>10 counts in ≥90% of the samples in that subset`, FDR by
Benjamini–Hochberg. These reproduce the baselines given in the task brief exactly
(lithium 12,173 genes / 1,382 DEG / π₀ 0.516; casecon 12,368 / 1,816 / 0.546;
lithium_krebs 12,207 / 1,459 / 0.503), which is the check that the engine is being
driven the same way here as elsewhere.

| contrast | n | genes | resid df | FDR .10 | **FDR .05** | FDR .01 | π₀ | % up | mean\|logFC\| (DEG) | max\|logFC\| |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `lithium` (BP1 only) | 226 | 12,173 | 213 | 2,447 | **1,382** | 495 | 0.516 | 69.2 | 0.205 | 0.914 |
| `casecon` (BD v control) | 474 | 12,368 | 460 | 2,721 | **1,816** | 762 | 0.546 | 41.5 | 0.142 | 0.728 |
| `lithium_krebs` (cases) | 240 | 12,207 | 226 | 2,578 | **1,459** | 499 | 0.503 | 67.7 | 0.203 | 0.919 |

Adding a fold-change floor thins these a lot, which is worth knowing before any
downstream list is taken at face value — at FDR .05 **and** |logFC| ≥ 0.2 the lithium
contrast keeps 557 of 1,382 and casecon keeps 249 of 1,816. The effects are small:
median |logFC| across *all* genes is 0.062 (lithium) and 0.051 (casecon).

The logFC distributions are near-symmetric and tight (`logFC_histogram.csv`,
quantiles in `canonical_summary.csv`): lithium q25/median/q75 = −0.062 / −0.009 /
0.062, casecon = −0.051 / −0.003 / 0.052. Long right tails in both.

π₀ near 0.5 implies ~5,900 non-null genes in the lithium contrast. That is a large
claim and it is **not** established here that all of it is lithium: any unmodelled
structure that tracks the exposure inflates the same statistic. Treat π₀ as an upper
bound on how much signal is present, not as a count of lithium-responsive genes.

Top genes are in `top30_per_contrast.csv`. The lithium top ten is TSPAN2, PIGB,
ARMC2, RFX2, RP4-666F24.3, LINC00877, FAR2, MYADM, MIR24-2, SLC31A2 — the same
leading genes, in nearly the same order, as the deposited Krebs table.

---

## 2. What the deposited Krebs tables actually contain

Both deposited DE tables carry **12,353 genes** with versioned Ensembl IDs. Two things
follow immediately:

- The "map symbols to Ensembl via the GTF" step the task anticipated is not needed.
  After stripping version suffixes, **12,353 of 12,353 match our row names exactly**.
  The GTF was still read, and independently confirms the annotation: it contains
  57,820 genes, exactly the deposited matrix row count.
- **12,353 is our number, not the paper's.** The paper's text says 12,344 filtered
  genes. Our filter on the Krebs 444 returns 12,353 — identical to their own deposited
  table. So the 9-gene discrepancy recorded in the project notes is a discrepancy
  *within the paper*, between its text and its supplementary file, not between the
  paper and us.
- And the match is on the **gene set, not merely the count** (`filter_gene_set_identity.csv`).
  Running our filter on the Krebs 444 gives a set identical to their deposited table:
  **0 genes in ours not theirs, 0 in theirs not ours.** Matching a count can happen by
  luck; matching all 12,353 identifiers cannot. The filtering step is exactly reproduced.
  It stops matching as soon as the sample set changes, which is the expected behaviour
  of a cross-sample filter (cases only: 213 of theirs missing; our 474 cohort: 24 extra).

DEG counts in the deposited files:

| file | FDR .05 | FDR .01 | FDR .10 | max\|logFC\| |
|:--|--:|--:|--:|--:|
| `Li_DEGs_PreCellTypeCorrection` | 3,031 | 1,578 | 3,907 | 1.026 |
| `Li_DEGs_PostCellTypeCorrection` | **233** | 131 | 348 | 0.673 |

The POST file's 233 matches the paper's stated 233 exactly, so that file is the
analysis it claims to be. The PRE file's 3,031 matches nothing in the paper (which
says 976 in Results and 897 in the MAGMA methods). Its max |logFC| of 1.026 rules out
the paper's 976-gene set being any subset of it, since the paper caps that set at 0.82.
This agrees with the conclusion already recorded in `../Krebs Paper/REPRODUCIBILITY.md`
— that the PRE sheet holds a different contrast than its label claims — and we reached
it independently from the file alone.

**Sign concordance is 100% for every specification tested**, across 668–1,393 genes
called significant by both us and the deposited sheet — recomputed from scratch with
disagreements counted explicitly, because a clean 100% is usually a bug:
**zero disagreements in total** (`sign_concordance.csv`). Among *all* shared genes,
not just the co-significant ones, agreement is 72–76%, which is the honest figure and
is what mostly-null genes should give.

---

## 3. The specification search — the main methodological result

### 3.1 `assessment group` is constant among cases, and that settles the model

Verified directly (`design_diagnostics.csv`):

```
          group A   group B
Control       125       109
BP1           226         0
BP2            14         0
```

- Among all 240 cases, `assessment group` takes **one value (A)**. Same on the Krebs
  444 subset (239 cases, all A).
- Group B contains **109 controls and zero cases**. Group is nested inside diagnosis.
- `lithium = 1` occurs only within BP1 (152 of 226); every control and every BP2
  subject is a non-user.

So the engine is right to drop `group` from `lithium_krebs`: within cases it carries no
information and makes the design singular. But the paper states it retained assessment
group in the lithium model. **Both statements can only be true if the published lithium
model was not restricted to cases.** The supplementary methods confirm it in a sentence
that is only meaningful if controls were present:

> "Although there were no controls being treated with lithium, diagnosis was included
> in the lithium-use comparison to account for BD-effects within the non-lithium using
> group."

**Implication: the engine's `lithium_krebs` is misnamed. It does not reproduce the
published contrast.** The published lithium fit ran on the whole cohort with diagnosis
as a covariate. Because lithium use is nested within BP1, the lithium coefficient is
still a within-BP1 comparison, but the gene filter, the normalisation, the mean–variance
trend, the covariate estimates, the residual degrees of freedom (428 vs 226) and the
empirical-Bayes moderation all come from 444 subjects rather than 240.

### 3.2 How diagnosis is coded matters more than anything else

Fitting the published specification on the Krebs 444, varying only how "BD diagnosis"
enters the model (`replication_vs_krebs.csv`):

| specification | n | genes | **DEG .05** | mean\|FC\| | max\|FC\| | sd\|FC\| | % up | carry-over of the 233 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| **Krebs, as published** | 444 | 12,344* | **976** | 0.20 | **0.82** | 0.10 | 77.3 | 194 (83.2%) |
| full cohort, diagnosis **2-level** (case/control) | 444 | 12,353 | **1,023** | **0.195** | **0.815** | **0.100** | 71.2 | **197 (84.5%)** |
| full cohort, diagnosis **3-level** (Control/BP1/BP2) | 444 | 12,353 | **2,652** | 0.179 | 0.952 | 0.095 | 56.6 | 196 (84.1%) |
| full cohort minus BP2 (either coding) | 431 | 12,377 | 2,591 | 0.180 | 0.952 | 0.096 | 58.6 | 196 (84.1%) |
| cases only, 444 (`lithium_krebs` on their set) | 239 | 12,151 | 1,450 | 0.204 | 0.921 | 0.106 | 68.4 | 193 (82.8%) |
| cases only, 474 (**engine `lithium_krebs`**) | 240 | 12,207 | **1,459** | 0.203 | 0.919 | 0.106 | 67.7 | 194 (83.3%) |
| BP1 only, 474 (**engine `lithium`**) | 226 | 12,173 | 1,382 | 0.205 | 0.914 | 0.107 | 69.2 | 192 (82.4%) |

\*paper's text; its own deposited file has 12,353.

The 2-level row matches the published fingerprint on four statistics at once —
n (1,023 v 976), mean |FC| (0.195 v 0.20), **max |FC| (0.815 v 0.82)** and sd (0.100 v
0.10) — and recovers 20/20 of the Anand replication genes and 130/134 of the Breen
genes. The 3-level row is decisively inconsistent with the published max of 0.82.
**The published lithium analysis coded diagnosis as a binary case/control indicator.**

### 3.3 Why that costs a factor of 2.5 — and it is 13 people

This is the result worth carrying forward. With a 3-level diagnosis the design is
saturated — one free parameter per cell — so the lithium coefficient is exactly
BP1-users minus BP1-non-users. With a 2-level indicator, the 13 bipolar-II subjects
(all lithium-free) carry the same covariate value as the 74 BP1 non-users, so the model
is forced to give both groups one shared mean. Least squares sets it to the
count-weighted average, and the lithium coefficient becomes

```
b_2level = y_use − (74·y_non + 13·y_BP2)/87
         = b_3level − (13/87)·(y_BP2 − y_non)
```

so the published lithium effect is measured against a reference group that is **15%
bipolar-II patients**. Removing those 13 subjects makes the two codings identical
(2,591 DEGs each), which is the proof that they are the whole difference.

`baseline_bp2.R` tests the two consequences of that algebra directly; see §7.

---

## 4. The case/control arm replicates essentially perfectly

The engine's `casecon` omits lithium use from the covariates; the published BD model
includes it. Fitting both (`casecon_decomposition.csv`):

| specification | n | FDR .10 | **FDR .05** | FDR .01 | π₀ | overlap with their 6 |
|:--|--:|--:|--:|--:|--:|--:|
| 474, lithium **not** a covariate (engine `casecon`) | 474 | 2,721 | **1,816** | 762 | 0.546 | 6/6 |
| 444, lithium **not** a covariate | 444 | 2,867 | **1,912** | 746 | 0.538 | 6/6 |
| 474, lithium a covariate | 474 | 13 | **1** | 1 | 0.757 | 1/6 |
| **444, lithium a covariate — the published BD spec** | 444 | 35 | **6** | 0 | 0.735 | 5/6 |

The published specification returns **exactly 6 DEGs**, as published, and their six genes
come back with almost the same fold changes (`bd_six_genes.csv`):

| gene | symbol | Krebs logFC | our logFC | our adj.P |
|:--|:--|--:|--:|--:|
| ENSG00000203872 | C6orf163 | 0.2582 | 0.2628 | 0.0190 |
| ENSG00000103051 | COG4 | −0.1163 | −0.1034 | 0.0536 |
| ENSG00000249859 | PVT1 | 0.2835 | 0.2816 | 0.0392 |
| ENSG00000088538 | DOCK3 | 0.4370 | 0.4262 | 0.0392 |
| ENSG00000267702 | RP11-53B2.2 | 0.3122 | 0.3139 | 0.0456 |
| ENSG00000122507 | BBS9 | 0.1606 | 0.1661 | 0.0392 |

All six sign-concordant; five clear FDR .05 and COG4 sits just outside at 0.0536. Our
sixth is RP1-102H19.8. This is as close to an exact reproduction as the deposit allows.

**The finding that matters for this project:** on the same 444 subjects, adding lithium
use as a covariate takes the case/control result from **1,912 DEGs to 6** — a 99.7%
collapse — and π₀ rises from 0.538 to 0.735, so the whole evidence distribution moves,
not just the genes near the threshold. The sample set is almost irrelevant by comparison
(1,816 v 1,912). The paper's headline contrast between a large lithium signature and a
near-absent diagnosis signature is produced by that one covariate.

---

## 5. Choices made here, and why

- **Contrasts added in memory only.** `li_case`, `li_dx`, `casecon_li`, `casecon_nol`
  and `cells` are attached to the in-memory copy of `prep`; `data/prep.rds` is never
  rewritten, so the other experiments running against it are unaffected.
- **Krebs' own QC flag is used as `samples =`**, not re-derived, so "the 444" means
  exactly the subjects they analysed.
- **Tobacco: `tobacco_imp_01` throughout, and it is inert where it matters.** Verified:
  all 444 QC-passing subjects have *observed* tobacco, and all 30 subjects with missing
  tobacco fall outside the 444. So no 444-based comparison in this document depends on
  the imputation, and the single-completion shortcut cannot explain any gap against
  Krebs. (Rubin-rules pooling across the 20 completions is still the right thing for the
  474-based fits; not done here, and flagged as a limitation.)
- **Sign concordance is reported two ways** — among genes both lists call significant,
  and among all shared genes — because the first is the generous reading and the second
  is the honest one.
- **Rank correlation uses signed −log10 p**, the same statistic Krebs use in their own
  BMI sensitivity analysis, so the scale is theirs.
- **The estimator was held fixed at voom + TMM** for every specification comparison, so
  that a difference between two rows is a difference in the model and never in the
  estimator. Estimator variation is tested separately, on the published spec only (§6),
  and deliberately not swept broadly — a sibling experiment owns that axis.
- **Seed** `20261218`, set once in `baseline_gap.R`, used only by the resampling in §6.

---

## 6. The residual gap, and estimator sensitivity

*(filled in by `baseline_gap.R` parts B and C — see below)*

---

## 7. Bipolar-II and the lithium axis

*(filled in by `baseline_bp2.R` — see below)*

---

## 8. What this does not establish

- π₀ ≈ 0.5 is an upper bound on signal, not a count of lithium-responsive genes.
- Nothing here separates a lithium effect from a cell-composition shift; that is the
  job of the deconvolution and mediation arms. The 1,912 → 6 collapse in §4 shows only
  that diagnosis and lithium are near-inseparable in this cohort's *expression* data.
- The one-subject deposit/publication discrepancy (Krebs count 240 cases / 204 controls,
  the deposit has 239 / 205) is characterised in §6 but cannot be resolved: the subject
  cannot be identified from public data.
- The 897 figure in the paper's MAGMA section matches no threshold on either deposited
  table and remains unexplained.
