# de_v2 — differential expression, no bands, full multiple-testing panel

The current differential expression analysis for GSE124326. Four comparisons,
every correction method applied to all of them, no verdicts.

`DECISIONS.md` is the reasoning. This file is how to run it and where things
are.

---

## The four comparisons

| | contrast | n | tests |
|:--|:--|--:|--:|
| **WB_LI** | whole blood; lithium users vs non-users within bipolar I | 226 | 12,368 |
| **WB_BPD** | whole blood; bipolar I off lithium vs healthy controls | 308 | 12,368 |
| **CT_LI** | per cell lineage; same contrast as WB_LI | 226 | 61,840 |
| **CT_BPD** | per cell lineage; same contrast as WB_BPD | 308 | 61,840 |

Each fitted twice: `_raw` (unadjusted) and `_ilr` (adjusted for the four
cell-composition balances). Eight fits, 148,416 tests in the raw set.

---

## How to run

```bash
source ../env.sh                               # R 4.6.1 on PATH
cd de_v2
Rscript scripts/00_validate_mtc.R              # ~2 min, gates everything below
Rscript scripts/01_refit_and_mtc.R             # ~5 min
Rscript scripts/02_hierarchical.R              # ~1 min
Rscript scripts/03_permutation_fdr.R           # ~3 h, checkpointed
```

`00` must print `ALL CHECKS PASSED` or nothing after it is trustworthy — it
proves the hand-written corrections in `mtc.R` are correct before they touch
real data. `03` resumes from `results/perm_parts/` if interrupted.

**Inputs** (read, never written): `../de_analysis/data/prep.rds` and
`../de_analysis/results/bmind_profiles.rds`. Nothing from v1's `results/` is
read. See `DECISIONS.md` §2.

---

## Which file do I want?

| question | file |
|:--|:--|
| How many genes, under every correction? | `results/summary_by_analysis.csv` |
| Everything about one gene | `results/tests_<analysis>.csv.gz` |
| Which lineage carries the signal? | `results/celltype_per_lineage.csv` |
| Cell type vs whole blood, fairly | `results/celltype_pooled.csv` |
| What does running four comparisons cost? | `results/across_contrast_correction.csv` |
| Are the corrections actually correct? | `results/mtc_validation.csv` |
| FDR without assuming independence | `results/permutation_fdr.csv` |
| All methods side by side | `results/correction_comparison.csv` |

Each `tests_*.csv.gz` has one row per test with `gene`, `lineage`, `logFC`,
`AveExpr`, `t`, `P.Value`, and then `bonferroni`, `holm`, `hommel`, `BH`, `BY`,
`q_storey_spline`, `q_storey_boot`.

---

## Headline numbers

Discoveries at the 0.05 level under each correction, unadjusted fits:

| | tests | π₀ | BH | Storey | BY | Holm |
|:--|--:|--:|--:|--:|--:|--:|
| WB_LI | 12,368 | 0.487 | 1,426 | 2,568 | 303 | 125 |
| WB_BPD | 12,368 | 0.644 | 4 | 5 | 0 | 0 |
| CT_LI | 61,840 | 0.773 | 105 | 120 | 33 | 22 |
| CT_BPD | 61,840 | 0.803 | 0 | 0 | 0 | 0 |

Correcting across all four comparisons (Benjamini–Bogomolov) selects **two of
four** families — both lithium ones. Neither bipolar family is selected, so
under that correction they report nothing.

Per lineage, lithium contrast: gran 109, mono 20, T 28, NK 4, B 0. The B-cell
family is never selected.

**Read the zeros carefully.** Every lineage's detection floor exceeds the
smallest effect anyone called meaningful — B cells would need a 72-million-fold
change. A cell-type zero here is not evidence of absence, and no correction
method can change that. `DECISIONS.md` §11.

---

## Relationship to the other analysis folders

| folder | status | what it is |
|:--|:--|:--|
| `de_v2/` | **current** | This. No bands, full correction panel. |
| `../de_analysis/` | superseded for results | v1, the pre-registered band-based analysis. `config.yaml` is the timestamped pre-registration and `METHODS.md` is still the reference for the gene filter, lineage choice, bMIND defects and detection floors. |
| `../de_analysis/_superseded/` | do not cite | Output of `bmind_de()`, abandoned after two verified defects. |
| `../_archive_exploratory_DE/` | archived | Pre-registration exploratory work. Not read by anything current. |

To find any file in the project: `python ../index/build_index.py --query "SELECT path, description FROM files WHERE path LIKE '%de_v2%'"`.
