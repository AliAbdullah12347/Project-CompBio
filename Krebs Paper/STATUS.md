# Status — saved 2026-09-27

Stopped mid-run for battery. Everything below is **on disk and complete**; nothing is
half-written. Pick up at "Next steps".

## Toolchain (installed once, machine-wide — do not reinstall)

- R 4.6.1 at `C:\Users\hp\R\R-4.6.1`, on the user PATH permanently.
- Library `C:\Users\hp\R\win-library\4.6`: limma 3.68.5, edgeR 4.10.5, WGCNA 1.74,
  statmod 1.5.2, MASS 7.3.65, dynamicTreeCut 1.63.1, fastcluster 1.3.0. Bioconductor 3.23.
- Shared by every project under `Implementation/` via `../env.sh` (bash) or `../env.ps1`.
- Chocolatey needs admin and installed nothing; R was installed per-user instead. Don't retry choco.

## Pipeline state — all scripts written and run clean

| Script | Status | Output |
|:--|:--|:--|
| `00_prepare_data.py` | done, 17/17 checks | `data/processed/` |
| `01_extract_ground_truth.py` | done | `data/reference/` (7 CSVs) |
| `02_differential_expression.R` | done | `DE_BD_all444.csv`, `DE_lithium_*.csv` |
| `02b_model_sweep.R` | done | console diagnostic |
| `02c_presheet_resolution.R` | done | `presheet_model_sweep.csv` (81 specs) |
| `03_wgcna.R` | done | `wgcna_modules.csv`, `wgcna_module_trait.csv` |
| `04_celltype.R` | done | `celltype_stepaic.csv`, `celltype_reference_comparison.csv` |
| `05_validate.py` | done | `validation_scorecard.csv` — 11 passed, 0 failed |

Total footprint 120 MB. `./download_data.sh` refetches all inputs; `./run_all.sh` runs everything.

## Reproduction scorecard

| Quantity | Paper | Ours |
|:--|:--|:--|
| Samples post-QC | 444 | 444 |
| Gene filter | "12,344" | **12,353 — matches their DE table exactly** |
| BD DEGs | 6 | 6; their six are our top 10/12,353; logFC r = 0.9998 |
| Lithium + cell types | 233 | 238; logFC r = 0.989 |
| Neutrophil → lithium (stepAIC, cases) | β 0.63, p 0.024 | **β 0.6297, p 0.0234** |
| Cell fractions vs reference (controls) | — | 38.79 / 36.65 / 24.56 / 0.00 % |
| WGCNA modules | 27 | 25 |

## The pre-sheet mystery: SOLVED

The deposited `Li_DEGs_PreCellTypeCorrection` sheet is **not** the lithium main effect. It is
the contrast **bipolar-on-lithium versus controls** (the case and lithium coefficients summed).

Best specification, from `results/presheet_model_sweep.csv` (81 specs tested):

```
contrast   : case + lithium  (BP on-Li vs control), all 444 samples
formula    : ~ case + age + sex + tobacco + group + rin + plate + seqpc1-3 + lithium
r(logFC)   = 0.9993     sign concordance 99.0%
n_sig      = 2,966  vs their 3,031     overlap 2,919  (Jaccard 0.948)
```

A 3-level diagnosis coding (`dx3`) scores marginally higher still (r = 0.99935, 2,903 overlap).

**Consequence:** the paper's text figure of 976 lithium DEGs does not correspond to any sheet
in its own deposit, and "the 976 lithium DEGs" should never be cited without saying which list
is meant. The deposited sheet conflates diagnosis with medication — exactly the confound the
paper set out to disentangle.

Supporting evidence, already computed: our documented lithium model flags **20/20 Anand** and
**130/134 Breen** externally-validated lithium-response genes as significant, versus 16/20 and
91/134 in their deposited sheet, despite having a third as many DEGs.

## WGCNA notes

25 modules, not 27. Correspondence to published modules is good but not exact:

| Ours | Published | Jaccard | r(lithium) | p |
|:--|:--|--:|--:|--:|
| turquoise (2551) | M1 (brown, 2092) | 0.34 | 0.170 | 3.2e-4 |
| green (462) | M26 (brown4, 484) | 0.29 | −0.186 | 7.9e-5 |
| pink (376) | M7 (lightsteelblue1, 700) | 0.43 | −0.176 | 1.9e-4 |
| orange (39) | M9 (plum2, 55) | 0.62 | 0.175 | 2.2e-4 |

Five modules are lithium-associated at Bonferroni, matching the paper's count of five, and the
turquoise/M1 DEG overlap is 662 (paper: 431 for M1) at p = 1.7e-229. Module count and sizes
differ because dynamic tree cut is sensitive to the exact residualisation and to WGCNA version.

## Next steps

1. **Run the adversarial verification pass** — it was the only thing still pending when the
   run was stopped. The three build agents finished; their verifiers did not run. Re-launch the
   workflow (script saved at the path in the session transcript) or spot-check by hand:
   re-run each script from a clean shell and independently recompute β = 0.6297 and r = 0.9993.
2. **Reconcile the WGCNA 25 vs 27.** Worth trying `deepSplit` 1–4 and checking whether
   residualising before vs after voom changes the cut.
3. **Fold findings into the pre-proposal / CLAUDE.md** — specifically the 12,353 correction
   (there is no annotation-level gap) and the fraction-sum tolerance (1e-4, not 1e-6).

## Corrections owed to our own project notes

- `CLAUDE.md` says the gene filter gives 12,353 vs "Krebs report 12,344 … a ~9-gene gap is
  expected and is annotation-level". **Wrong.** Their DE table has 12,353; 12,344 is their
  WGCNA count. No annotation gap exists.
- `CLAUDE.md` says the deposited fractions sum "to 1.000000". They miss by up to 5.8e-5;
  a 1e-6 check wrongly flags 106 of 480 samples.
- Krebs' lithium contrast is **all 240 cases (152 vs 88)**, and keeps `ascertainment group`,
  which is constant within cases. Our BP1-only contrast (74 vs 152) is a different, smaller
  comparison and must drop that covariate or the design is singular.
