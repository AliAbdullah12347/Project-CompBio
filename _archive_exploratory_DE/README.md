# experimentation/

Exploratory differential expression work on GSE124326. Self-contained: nothing
here is read by the pipeline in `../data/`, and nothing here writes outside this
folder. Delete the whole directory and the project is unaffected.

**This is exploratory.** None of it is pre-specified, several analyses were
chosen after seeing earlier results, and the thresholds were deliberately varied
rather than fixed. Results are hypotheses and robustness checks, not
confirmatory findings. Anything promoted to the paper has to be re-declared in
`config.yaml` first and re-run on that basis.

## Layout

```
experimentation/
├── scripts/      every analysis, numbered in the order it was run
├── data/         prep.rds -- the shared objects all scripts start from
├── runs/         one folder per analysis: summary.csv, NOTES.md, artefacts
└── FINDINGS.md   what was found, with the numbers
```

## Foundation

`scripts/00_prep.R` builds `data/prep.rds` once, and every other script starts
from it, so a difference between two analyses is a difference in the analysis
and never a difference in how the data was made. It contains:

| Object | What |
|:--|:--|
| `counts` | 57,773 × 474 raw counts, ENSGR rows dropped, versions stripped |
| `meta` | 474 rows: diagnosis, lithium, covariates, QC flags |
| `TOB` | all 20 completed tobacco imputations |
| `frac` | 22 CIBERSORTx fractions (B-mode run) |
| `lineage`, `lineage_z` | fractions aggregated to 5 lineages, zero-replaced |
| `ILR` | 4 isometric log-ratio balances on a verified orthonormal SBP |
| `CONTRASTS` | the three contrast definitions |

`scripts/common.R` is the single DE engine. Every experiment calls `de_fit()`
and changes exactly one argument, so results are comparable by construction.

### Choices made in the foundation, and why

**Cell composition enters at lineage level, not 22-type.** Fine-grained LM22
output on this dataset has a between-method ICC of 0.389, and the Aitchison
distance between two deconvolution methods is 1.56× the spread between subjects
— method choice moves the composition more than biology does. Aggregating to
five lineages raises median ICC to 0.693 and the myeloid/lymphoid balance to
0.907.

**Zero replacement is multiplicative at 65% of the detection limit**
(Martín-Fernández 2003), because it preserves the ratios among non-zero parts
and every balance is a ratio.

**No depth/RIN screen.** Krebs et al. state twice that they removed no sample
for quality, adjusting instead. The flags are in `meta` so the screened subset
is one line away.

**Assessment group is dropped automatically where it is constant.** All cases
are group A, so it is inestimable in any cases-only design. The engine detects
and reports this rather than producing a silently rank-deficient fit.

## Contrasts

| Name | n | Definition |
|:--|--:|:--|
| `lithium` | 226 | bipolar I only, 74 non-users vs 152 users |
| `casecon` | 474 | all samples, BD vs control |
| `lithium_krebs` | 240 | all cases, diagnosis retained — the published design |

Plus four built in `90_signature_overlap.R` that the shared list does not cover,
each needing its own exposure variable:

| Name | Comparison |
|:--|:--|
| `A_off_vs_ctrl` | BP1 **not** on lithium (74) vs control (234) |
| `B_on_vs_ctrl` | BP1 **on** lithium (152) vs control (234) |
| `C_lithium` | lithium within BP1 |
| `D_allBD_vs_ctrl` | all BD vs control |

## Scripts

Run in this order; each is independent once `00_prep.R` has run.

| Script | Question |
|:--|:--|
| `00_prep.R` | build shared objects |
| `common.R` | the DE engine (sourced, not run) |
| `90_signature_overlap.R` | is the bipolar signature separable from lithium? |
| `91_power_matched.R` | is that separation just a sample-size artefact? |
| `92_is_it_just_composition.R` | closed-form aggregate mediation test |
| `93_annotate.R` | map Ensembl IDs to symbols, interpret gene lists |
| `94_hidden_celltypes.R` | is the "direct effect" composition LM22 cannot see? |
| `95_eos_abundance_vs_induction.R` | cell abundance or within-cell induction? |
| `96_krebs_replication.R` | replicate the source paper; test its gene list off lithium |

| `98_interaction_surrogate.R` | exposure x mediator interaction gate; hidden structure |
| `99_robustness.R` | covariate LOO, sample set, tobacco MI, p-value histograms |

Agent-written scripts (`filter_sweep.R`, `threshold_sweep.R`, `norm_sweep.R`,
`method_sweep.R`, `cellcomp.R`, `ref_krebs.R` and their follow-ups) cover the
threshold, method, normalisation and spike-in axes; each documents itself in
`runs/<name>/NOTES.md`.

A session usage limit killed the 14 parallel agents at their reporting step,
after their analyses had completed. Six left full output (`baseline`,
`cellcomp`, `filter_sweep`, `method_sweep`, `norm_sweep`, `threshold_sweep`);
the rest were rewritten directly as `98_` and `99_`. No adversarial verification
pass ran, so agent-produced numbers in `FINDINGS.md` are single-sourced and
marked where they are.

## Running it

```bash
source ../env.sh          # puts R 4.6.1 on PATH
Rscript scripts/00_prep.R
Rscript scripts/90_signature_overlap.R
```

R 4.6.1 with limma 3.68.5, edgeR 4.10.5, statmod, WGCNA, ggplot2. `sva`,
`DESeq2` and `qvalue` are **not** installed and were not installed — surrogate
variables are built by hand from residual PCA instead, and π₀ is estimated with
a hand-rolled Storey estimator in `common.R`.

## Caveats that apply to everything here

**"Off lithium" is not "unmedicated."** The dataset records lithium use and
nothing else — no antipsychotics, anticonvulsants, antidepressants, dose or
duration. The 74 bipolar I subjects not on lithium may well be on other
psychotropics. Every comparison involving that group carries this limitation.

**The mediator is estimated from the outcome matrix.** Cell fractions come from
the same expression data that is the DE outcome, so adjusting one for the other
is partly self-referential. `94` quantifies one consequence directly.

**Multiple thresholds were tried on purpose.** That is the point of the folder,
and it also means no p-value here is corrected for the number of analyses run.
