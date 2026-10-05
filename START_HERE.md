# Start here

Orientation for a new session. Read this, then `de_v2/DECISIONS.md`.

**The question.** Almost every bipolar patient in a public dataset is medicated,
so a "bipolar signature" in blood may just be a lithium signature. This project
separates them in GEO **GSE124326** (474 subjects, 12,368 genes after filtering).

---

## Read in this order

| | |
|:--|:--|
| `de_v2/DECISIONS.md` | 16 sections. Every analysis decision and why. **Start here.** |
| `de_v2/README.md` | What the current analysis is, how to run it, which file answers what. |
| `de_analysis/METHODS.md` | 21 sections. v1's log, still the reference for the gene filter, lineage choice, the `bmind_de` defects and the detection floors. |
| `de_analysis/config.yaml` | The pre-registration. Historical record — v2 deliberately drops its bands. |
| `DIRECTORY.txt` | Hand-written tour of every folder. |

---

## Navigate with the index, not by guessing

```bash
python index/build_index.py                                  # rebuild (<1s)
python index/build_index.py --query "SELECT path, description FROM files WHERE path LIKE '%sweep%'"
```

Views: `v_current_outputs`, `v_code`, `v_do_not_use`, `v_biggest`, `v_undescribed`.
If `DIRECTORY.txt` and `INDEX.db` disagree, **the index is right** — it is generated.

---

## Which folder is live

| folder | status |
|:--|:--|
| **`de_v2/`** | **current.** No verdict bands, seven correction methods, filter sweeps. |
| `de_analysis/` | superseded for results. Keep: it holds the pre-registration and the detection floors. |
| `de_analysis/_superseded/` | **do not cite.** Output of `bmind_de()`, abandoned after two verified defects. |
| `_archive_exploratory_DE/` | archived. Read by nothing current. Not for citation. |

---

## Five facts that prevent mistakes

1. **No cell-type zero is evidence of absence.** Every lineage's detection floor
   exceeds the 0.5 log2FC meaningfulness threshold — B cells would need a
   72,000,000-fold change. This applies to our own results. `DECISIONS.md` §11.
2. **DEG counts are not comparable across filters.** The denominator *and* the
   multiple-testing burden both move. Use the per-gene stability tables in
   `de_v2/results/sweep/`, not counts. §15.
3. **π₀ is unusable at cell-type level.** Its permutation null is centred near 1
   but so dispersed (5th percentile 0.073) that a low value means nothing. §19
   of METHODS.md. Use it only to power the Storey adjustment.
4. **Don't pool all four comparisons into one BH.** It *raises* the bipolar
   counts (4→16, 0→22) because lithium's signals lift the threshold for
   everyone. Use Benjamini–Bogomolov. §8.
5. **"Off lithium" is not "unmedicated."** The deposit records lithium and
   nothing else. Every bipolar claim inherits this.

---

## Headline results

| | genes at BH 0.05 | |
|:--|--:|:--|
| Lithium, whole blood | **1,426** | robust core of 1,092 under every filter |
| Bipolar, whole blood | **4** | zero survive Rubin pooling, TREAT, or all filters |
| Lithium, cell type | **66** | 79 granulocyte findings hold under every threshold |
| Bipolar, cell type | **0** | but see fact 1 |

92.8% of the lithium genes move in the direction granulocyte fraction predicts;
adjusting for composition shrinks effects by a median 40.8%. The effect that
remains is in granulocytes.

---

## Getting the data

`data/analysis_matrices/` has everything needed to reproduce without re-running:
filtered counts, log2(TPM+1), TMM factors, CIBERSORTx fractions (22 types and 5
lineages), ILR coordinates, all 20 tobacco imputations, covariates, contrast
membership. Files over 14 MB are split into `.partN.csv`.

Two things are **not** in the repo and must be rebuilt:

```bash
python data/make_mixture_full.py          # the 102 MB CIBERSORTx matrix (exceeds GitHub's limit)
Rscript de_analysis/scripts/06_bmind_profiles.R   # bMIND cell-type estimates, ~3.1 h
```

Journal PDFs are excluded too — not ours to redistribute.

---

## Open items

- **Prediction arm (Arm 1) not started.** Design settled: nested repeated CV,
  everything cross-sample in-fold, classic CIBERSORT (not our B-mode fractions,
  which are transductive), permutation null, DeLong with a 0.05 AUC equivalence
  margin. Needs `glmnet`, `pROC` and a tree learner — **none installed**, and
  CLAUDE.md says ask before adding a dependency.
- **CIBERSORTx high-resolution** run in progress; see
  `data/cohort_474/CIBERSORTX_HIRES.md` for settings and the gene-subset files.
- **Correction needed:** CLAUDE.md §1 says "five compositional coordinates".
  Five lineages give **four** ILR balances. Fix before it reaches the paper.

---

## Hard rules

Everything lives in `Implementation/`. Nothing outside it is read or written
without explicit permission. Ask before adding a dependency. Verify numbers by
computing them — do not carry a figure forward from a paper or from notes.
