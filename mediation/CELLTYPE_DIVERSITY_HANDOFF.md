# HANDOFF — alpha and beta diversity of the cell-type-specific estimates

Paste everything below the line into Claude Code on the supercomputer, inside a checkout of
`https://github.com/AliAbdullah12347/Project-CompBio` (the `Implementation/` folder is the repo root).

---

You are running one pre-specified analysis for an undergraduate computational-biology project.
**The design is fixed. Do not change any analytic decision**; if something seems wrong, record it
in `mediation/results/celltype_diversity/ISSUES.md` and carry on with the declared design.

## 1. Get the code and data

```bash
git fetch origin
git checkout mediation-circularity-diversity   # skip if it has been merged into main
git pull
```

Everything needed is tracked — no download, no `.rds` cache:

| Input | Path |
|:--|:--|
| bMIND cell-type expression, log2(TPM+1), 12,368 genes × 474 samples, one file per lineage | `data/analysis_matrices/bmind/bmind_logtpm_{gran,mono,T,NK,B}.csv.gz` (checksums in `CHECKSUMS.sha256`) |
| bMIND precision per sample × lineage | `data/analysis_matrices/bmind/bmind_mean_posterior_se.csv` |
| counts, covariates, ILR balances, fractions | loaded by `network/00_load.R` from `data/analysis_matrices/` |
| gene biotypes, frozen MSigDB 2026.1.Hs gene sets | `mediation/data/gencode_v19_gene_biotypes.tsv`, `mediation/data/genesets_msigdb_2026.1.Hs.csv` |
| the estimator | `mediation/scripts/med.R` |
| **the script you run** | `mediation/scripts/16_celltype_diversity.R` |
| the protocol it implements | `mediation/PRESPEC.md`, **section 9** — read it first |

Verify the data: `cd data/analysis_matrices/bmind && sha256sum -c CHECKSUMS.sha256` — all must say OK.

**Never commit** `.env`, `env.sh`, `env.ps1` or any `kubeconfig*` file. They are in `.gitignore`;
do not override it, and do not `git add -f` anything.

## 2. Environment

R ≥ 4.3 with `vegan`, `permute`, `limma` (Bioconductor) and `edgeR`. Install them if missing; nothing else.

## 3. Run

The first thing the script does is reproduce the published whole-blood result. If it prints
`GATE FAILED` it stops — **do not work around it**; report what differs.

**Smoke test first** (about 3 minutes):

```bash
LINEAGES=gran NBOOT=50 NPERM=50 NPERM_V=99 NBOOT_V=20 OUT=/tmp/ct_smoke \
  Rscript mediation/scripts/16_celltype_diversity.R
```

It must print `GATE PASSED` and finish lineage `gran`.

**Full run.** The five lineages are independent. Run them as five parallel jobs (one per lineage,
any scheduler), each with the default settings:

```bash
LINEAGES=gran NCORES=4 Rscript mediation/scripts/16_celltype_diversity.R   # likewise mono, T, NK, B
```

then merge and apply the pre-declared Holm correction:

```bash
COMBINE_ONLY=1 Rscript mediation/scripts/16_celltype_diversity.R
```

Defaults are the declared values: 10,000 bootstrap and 10,000 permutations for alpha, 9,999
permutations and 1,000 bootstrap replicates for beta. **Do not lower them for the real run.**
Expect roughly 20–40 minutes per lineage on one core; the Jensen–Shannon distance and the vegan
permutation tests dominate. `NCORES` speeds up the permutation tests.

## 4. Outputs (all in `mediation/results/celltype_diversity/`)

| File | Holds |
|:--|:--|
| `alpha_<lineage>_<exposure>.csv` | per lineage and exposure: total effect on Shannon (and Hill 2, markers-in, Hallmark pathway Shannon), the same adjusted for bMIND precision, the indirect effect through the cell mix with bootstrap CI and permutation p, interaction p |
| `beta_<lineage>_<exposure>.csv` | PERMANOVA (B1), the same adjusted for precision, PERMDISP raw and residualised (B2), direct effect given the cell mix (B3), shared fraction with CI, path b; Aitchison and Jensen–Shannon layers |
| `ALPHA_ALL.csv`, `BETA_ALL.csv` | everything merged, with Holm-adjusted p for the primary tests |
| the run logs | save stdout of every job as `log_<lineage>.txt` in the same folder |

## 5. When it finishes

1. Commit the results folder and logs to the same branch with a clear message, and push.
2. Write `mediation/results/celltype_diversity/SUMMARY.md` containing, without interpretation
   beyond section 6: the two primary tables the combine step prints; the gate line; wall-clock
   per lineage; any warnings or issues.
3. Report back in this chat with the two primary tables.

## 6. How to read the results — fixed in advance (PRESPEC section 9.3)

| Result | Reading |
|:--|:--|
| no effect in gran, mono, T | lithium's diversity effect in blood is compositional — consistent with the whole-blood decomposition |
| an effect that **survives** the precision adjustment | lithium changes that lineage's transcriptome (cell-intrinsic, or a sub-type shift inside the lineage) |
| an effect that **disappears** with the precision adjustment | an artefact of bMIND shrinkage, not biology |
| a non-zero indirect effect through the cell mix | the estimates still carry composition — a caveat on that lineage |
| anything in NK or B | uninformative: their detection floors (4.07, 26.11 log2FC) are too high |

Primary confirmatory tests: Shannon total effect (A1) and Aitchison PERMANOVA (B1) for **gran, mono, T**,
Holm-corrected within each test, lithium. Bipolar vs control is **exploratory** throughout.

## 7. Context, so the numbers make sense

The same analysis on **whole blood** (`MEDIATION_MASTER.md`, section 12) found that lithium lowers
transcriptome Shannon diversity through the cell mix (indirect effect −0.034, permutation p = 0.020),
95 % of it from the blood reweighting toward granulocytes, whose transcriptome is the narrowest
(~1,400 effective genes vs ~3,600 for lymphocytes). Within-lineage changes were not detectable there.
This run tests that last statement directly, one lineage at a time.
