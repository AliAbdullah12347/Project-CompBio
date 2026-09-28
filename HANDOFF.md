# Handoff — session of 19–28 Sep 2026

Context for the next session. Excludes the Krebs reproduction work, which is
self-documenting in `Krebs Paper/STATUS.md`, `README.md` and `REPRODUCIBILITY.md`.

---

## 1. The hard rule

**All code and changes go in `Implementation/`. Do not read, write, or modify anything
outside it without explicit written permission from the user — this holds even when bypass
permissions is toggled on.**

Permission has been granted, file by file, for exactly three files outside it:
- `AliAva_PreProposal_Feedback.pdf`
- `Ali, Ava Project Preproposal.docx`
- `Papers- Annotated Bibliography/Whole-blood-transcriptome-analysis-...pdf`

Ask before touching anything else, including directory listings.

---

## 2. Pre-proposal revision — done

Professor Ay returned the pre-proposal with **~35 inline comments plus two summary notes**
(PDF annotation author: `ahmetay`). His verdict: structure passed, the dataset section was
praised as "essentially complete", and the plan was explicitly credited for handling data
leakage and compositional data. Four named defects and one existential question.

**Deliverable:** `preproposal/Ali, Ava Project Preproposal (revised).docx`
- Original preserved beside it as `original.docx` — the revision was written to a new file,
  not over the original.
- Reproducible: edit prose in `preproposal/content.py`, re-run `build.py`. `verify.py`
  checks all 51 feedback items are present and passes 51/51.
- Grew ~1,900 → ~9,100 words (user said to ignore the page limit). Seven numbered sections,
  31 references, a pre-specification parameter table, four captioned figures.

### What he cared about most (repeated concerns)
1. **Novelty vs Boltz et al. 2024** — raised three times. "Just checking mediation is not a
   full project." Answered in a new §1.1 with four specific contributions.
2. **Circularity** — raised three times. The cell fractions are estimated from the same
   expression matrix that is the mediation outcome. Answered with marker exclusion,
   stratified reporting, cross-method sensitivity, and an honest statement of the residual.
3. **Prediction ≠ mediation** — he asked directly "which result would separate the two?"
   Answered by adding a third Arm 1 comparator: expression residualised on the ILR coordinates.
4. **Pre-specification with actual numbers.** He supplied two himself — ΔAUC = 0.05 and
   Zsummary < 2 — and both were adopted verbatim.

### Factual corrections made to our own draft
- Heritability redefined as share of variance in liability, not of individual risk.
- **Boltz was misstated in our draft.** They adjusted for cell-type proportions and *still*
  detected lithium-associated genes. That makes the question quantitative, not binary, and
  means a partially direct effect is the honest expectation.
- 12,343 → 12,353 (we had three different numbers across the document).
- Imputation sentence withdrawn: all covariates are complete post-QC.
- Dropped PPIXpress (wrong tool — protein interaction networks, not co-expression centrality;
  also miscited as "Will & Helms" while the reference list had the authors' names swapped).
- Dropped the borderline-personality methylation reference. Added BayesPrism (Chu et al. 2022).

### Two judgment calls made without the user (both stated explicitly in the text)
- **The 32 recoverable controls** — primary analysis stays at 444. Recovering them via MICE is
  a pre-specified *sensitivity analysis only*, with his three rules written in: never impute
  exposure or outcome, fit imputation models inside training folds only, and the four sample
  mix-ups stay excluded permanently.
- **H4 (age × lithium moderated mediation)** — kept, but labelled "exploratory, declared
  underpowered", with the power argument spelled out. Cutting it silently would have looked
  like hiding it.

Either can be reversed; they are flagged as decisions, not facts.

### The one thing not done
**The four flowchart PNGs still contain the errors he flagged.** Raster images can't be edited
and matplotlib is not installed. Handled by making the text authoritative and listing each
diagram's corrections in its caption — the missing inner tuning loop in Fig 3, the "500-100
draws" typo in Fig 4, the withdrawn imputation step and "use the best method" box in Fig 1,
the missing interaction test in Fig 2. **An offer to redraw all four is outstanding** and
would need matplotlib or graphviz installed.

---

## 3. Toolchain — installed once, machine-wide

R lives in the user profile, not in any project, so it is installed once and shared by
everything under `Implementation/`.

- **R 4.6.1** at `C:\Users\hp\R\R-4.6.1`, added to the **user PATH permanently**.
- Library at `C:\Users\hp\R\win-library\4.6` (also set as `R_LIBS_USER`):
  limma 3.68.5, edgeR 4.10.5, WGCNA 1.74, statmod 1.5.2, MASS 7.3.65,
  dynamicTreeCut 1.63.1, fastcluster 1.3.0. Bioconductor 3.23.
- Load it with `source Implementation/env.sh` (bash) or `. Implementation\env.ps1` (PowerShell).

**Python 3.13** has numpy 2.4.1, pandas 2.3.3, networkx 3.6.1, pypdf 6.18.0.
**Missing and likely needed: `scipy`, `statsmodels`, `scikit-learn`, `matplotlib`** (~300 MB).
pip reaches PyPI fine. Not installed — CLAUDE.md requires asking first, and the user has not
yet approved (see §5).

---

## 4. Environment quirks worth knowing

These each cost time to discover:

- **The Read tool falsely reports PDFs as password-protected.** Both PDFs in this project
  trigger it and neither is encrypted. Extract with `pypdf` instead.
- **Chocolatey needs admin and silently reports success while installing nothing** — it
  returned exit 0 having installed 0/0 packages. Don't retry it; install to the user profile.
- **Not available:** pandoc, LibreOffice/soffice, openpyxl, defusedxml. The bundled docx
  skill's `merge_runs.py` and `validate.py` both fail on the missing defusedxml, so docx and
  xlsx work was done by parsing OOXML directly with `zipfile` + `ElementTree`.
- **Bash heredocs hit a command-length limit** (`ENAMETOOLONG`) well before a long file fits.
  Use the Write tool for anything substantial.
- **OOXML enforces a strict child order inside `<w:pPr>`** (`pStyle → keepNext → keepLines →
  spacing → ind → jc → outlineLvl → rPr`). Emitting out of order produces a file Word may
  open but strict readers reject.
- **Every sample title in GSE124326 starts with a digit**, so R's `make.names` prefixes each
  with `X`. Always read with `check.names = FALSE` and join on the title key, never on column
  order. The metadata `characteristics` lines are also ragged — parse each `field: value` pair
  by its own key, never by row index.

---

## 5. Open decisions

1. **Language architecture.** CLAUDE.md names Python modules (`load.py`, `arm1_prediction.py`)
   but the reference tools are split: limma-voom and WGCNA are R, while elastic net / random
   forest / nested CV are much better in scikit-learn. **Recommended: hybrid** — R for voom,
   WGCNA and BayesPrism; Python for the ML and mediation. Needs the four Python packages in §3.
   *Awaiting user approval.*
2. **Whether to redraw the four flowcharts** (§2).

---

## 6. Recommended next steps

The user's stated priority is getting real results fast. Recommended order:

1. **`config.yaml` first, before anything else.** Not the interesting work, but it is what
   makes every later result legitimate — Ay's sharpest demand was that a margin chosen after
   the results is not a test. Committing ΔAUC = 0.05, Zsummary < 2, proportion mediated > 0.5,
   depth ≥ 3M, RIN ≥ 5, the zero-replacement rule and the seed to version control, with a
   timestamp that provably precedes any result, converts the project from "we promise we
   pre-specified" to "here is the commit". Cannot be retrofitted.
2. **`load.py`** — roughly 80% already written and validated in `Krebs Paper/src/00_prepare_data.py`.
   Port it, then add the two things it lacks for our project: our own depth/RIN screen, and
   BP1-only subsetting for the lithium contrast.
3. **`tests/test_facts.py`** — 17 dataset facts already verified against the real data. Write
   them as assertions before any analysis code, per CLAUDE.md's ordering.
4. **`deconvolve.py` through to ILR balances** — the genuinely new work (structural zeros,
   zero replacement, lineage tree), and where all three arms converge, so it unblocks the most.
5. **Arm 2 (mediation)** then gives a real result quickly. CLAUDE.md calls it the cheapest real
   result and that is correct.

**Do not look at Arm 1 numbers before `config.yaml` is committed**, or the claim that the
margin was pre-specified is lost.

---

## 7. Corrections owed to CLAUDE.md

Three stated facts were checked against the data and are wrong. They should be fixed in
`CLAUDE.md` so they stop propagating:

1. **The gene filter gives 12,353, and there is no annotation-level gap.** CLAUDE.md says
   "Krebs report 12,344. A ~9-gene gap is expected and is annotation-level, not a pipeline
   error." In fact our filter reproduces the source paper's own deposited DE table gene-for-gene,
   12,353 of 12,353, with zero difference. Their 12,344 is a *later* gene count, after a
   further nine genes drop at the network stage. No annotation gap exists.
2. **The deposited cell fractions do not sum to 1.000000.** They miss by up to 5.8 × 10⁻⁵ from
   rounding. A 10⁻⁶ tolerance wrongly flags 106 of 480 samples; 10⁻⁴ is the correct check.
3. **The source paper's lithium contrast was all 240 cases (152 vs 88), not BP1-only**, and its
   model retained `ascertainment group`. That variable is constant within cases, so our
   BP1-only contrast (74 vs 152) is a different and smaller comparison and **must drop that
   covariate or the design matrix is singular.** CLAUDE.md already warns about the singularity;
   what it does not say is that the published contrast is a different one from ours.

Everything else in CLAUDE.md §2 that was checkable was verified and holds — 17/17, including
the 444 QC count, the 74/152 BP1 lithium split, the Fisher p = 3.04 × 10⁻¹⁰ on the
non-random QC exclusion, the 47 ENSGR rows totalling 1,518 counts, the four structural-zero
cell types, and zero covariate missingness post-QC.

---

## 8. Scientific warning carried forward

The source paper found neutrophils predict lithium use with **adjusted R² = 0.028** (β = 0.63,
p = 0.024) — a weak effect, and we reproduced it almost exactly. If five ILR coordinates carry
about that much signal while the gene matrix has over a thousand genes at FDR < 0.05, the gene
set may beat the coordinates by more than the pre-declared 0.05 AUC margin, meaning **H1 is
rejected.**

That is a perfectly publishable answer — the transcriptome carries information a blood count
does not — but it should be expected rather than come as a surprise, and it is a reason to get
`config.yaml` committed before anyone looks.
