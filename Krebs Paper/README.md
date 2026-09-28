# Krebs et al. 2020 — reference recreation

A from-scratch reimplementation of

> Krebs CE, Ori APS, Vreeker A, Wu T, Cantor RM, Boks MPM, Kahn RS, Olde Loohuis LM, Ophoff RA
> (2020). **Whole blood transcriptome analysis in bipolar disorder reveals strong lithium
> effect.** *Psychological Medicine* 50(15), 2575–2586.
> [doi:10.1017/S0033291719002745](https://doi.org/10.1017/S0033291719002745)

The authors released data but **not code**. This directory rebuilds their analysis from the
published methods and checks every stage against their own deposited result tables, so that
our project has a trustworthy reference for what this dataset does and does not support.

Read [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) before relying on any result here — several of
the paper's analyses cannot be reproduced from public data, and one of its headline numbers
does not match its own deposit.

---

## Layout

```
Krebs Paper/
├── README.md               ← this file
├── REPRODUCIBILITY.md      ← what is and is not reproducible, and why
├── download_data.sh        ← fetch all inputs (~29 MB)
├── run_all.sh              ← run the whole pipeline
├── paper/
│   └── Krebs2020_PsychologicalMedicine.pdf
├── docs/
│   ├── paper_fulltext.txt          ← extracted, for grepping the methods
│   └── supplementary_methods.txt   ← the exact covariate models
├── data/
│   ├── raw/                ← GEO downloads, never edited (+ SHA-256 checksums)
│   ├── supplementary/      ← Cambridge supplementary files (+ checksums)
│   ├── reference/          ← the authors' published tables, as CSV — the validation target
│   └── processed/          ← derived; regenerable from src/
├── src/
│   ├── 00_prepare_data.py          ← parse GEO, QC, gene filter, fact checks
│   ├── 01_extract_ground_truth.py  ← supplementary xlsx → CSV
│   ├── 02_differential_expression.R
│   ├── 02b_model_sweep.R           ← diagnostic: which lithium model was deposited?
│   ├── 02c_presheet_resolution.R   ← the full search behind that question
│   ├── 03_wgcna.R
│   └── 04_celltype.R
└── results/                ← all output tables
```

## Running it

R 4.6.1 with limma, edgeR, WGCNA, statmod and MASS is installed once per machine under
`C:\Users\hp\R\` and shared by every project anywhere under `Comp Biology Research/`.
`~/.Renviron` pins that library path, so R finds the packages from any working directory
with nothing sourced. To also put `Rscript` on `PATH`:

```bash
source "../../env.sh"      # or: . "..\..\env.ps1"  from PowerShell
```

Then:

```bash
./download_data.sh   # ~29 MB, idempotent, verifies checksums
./run_all.sh
```

Total on-disk footprint is about 79 MB including all derived data.

## The dataset in one paragraph

GSE124326 is whole-blood RNA-seq from 480 subjects in a Dutch bipolar cohort: raw HTSeq union
counts for 57,820 genes (exactly the GENCODE v19 total). Forty-seven `ENSGR` rows are
pseudoautosomal Y duplicates and are dropped, leaving 57,773. The deposited quality flag keeps
444 samples. Lithium use is recorded only for cases, and the 22 CIBERSORT LM22 cell fractions
are deposited in the series matrix, so CIBERSORT does not need re-running.

Two traps worth knowing about, both handled in `00_prepare_data.py`:

- **Every sample title starts with a digit**, so R's `make.names` silently prefixes each with
  `X`. Always read with `check.names = FALSE` and join on the title key, never on column order.
- **The characteristics lines are ragged** — a field does not occupy the same line index for
  every sample (`Sex` is present for 478 of 480). Parse each `field: value` pair by its own
  key.

## What we found that the paper does not say

Three things surfaced from checking against the authors' own deposit:

1. **The gene filter yields 12,353 genes, not 12,344.** Our filter reproduces their deposited
   DE table gene-for-gene, 12,353 out of 12,353, with zero difference. The paper's stated
   12,344 is the *WGCNA* gene count — nine genes are dropped between the DE and network stages.
   Earlier notes of ours describing this as an unexplained annotation-level gap were wrong;
   there is no gap.

2. **The deposited lithium DEG sheet does not match the documented lithium model.** The paper
   says 976 DEGs, its MAGMA section says 897, and the deposited sheet has 3,031 at FDR < 0.05.
   Details in `REPRODUCIBILITY.md`.

3. **The deposited cell fractions do not sum to 1 exactly** — they miss by up to 5.8 × 10⁻⁵
   from rounding. A 1e-6 tolerance wrongly flags 106 of 480 samples; 1e-4 is the right check.

## Relevance to our project

- The lithium contrast Krebs ran is **all 240 cases, 152 users vs 88 non-users**, and their
  model keeps `ascertainment group` — which is constant within cases. Our project's BP1-only
  contrast (74 vs 152) is a *different and smaller* comparison, and must drop that covariate
  or the design is singular.
- Their neutrophil result (elevated in lithium users) is the positive control our deconvolution
  has to reproduce before anything downstream is trustworthy.
- The LM22 misfit to whole blood — monocytes at roughly three times their reference upper
  bound, four cell types structurally zero — is visible in their own deposited fractions, not
  an artifact introduced by us.
