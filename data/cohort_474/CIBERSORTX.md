# Running CIBERSORTx on this cohort

Written 28 September 2026. Keep this alongside the results — the settings are
not recoverable from the output files, and using the wrong run in the wrong arm
invalidates it silently.

## The file to upload

```
Implementation\data\cohort_474\mixture_TPM_hgnc_full_474.txt
```

55,765 gene symbols × 474 samples, 107 MB.
SHA-256 in `mixture_full_474.sha256`.

A smaller alternative, `mixture_TPM_hgnc_474.txt` (20,628 symbols, 49 MB),
exists for the case where the web portal rejects 107 MB. For a run **without**
batch correction the two give identical answers, because the support vector
regression only ever touches the intersection with LM22's 547 signature genes.
For a run **with** B-mode batch correction, prefer the full file — B-mode
corrects across the whole mixture matrix, not just the signature genes.

### Verified, 12 of 12

| Check | Result |
|:--|:--|
| Gene rows | 55,765 |
| Sample columns | 474, identical to `sample_list.txt` and the metadata, in order |
| Format | tab-delimited, `GeneSymbol` header |
| Identifiers | HGNC symbols, zero Ensembl IDs, zero duplicates |
| Scale | linear (non-log) TPM |
| Values | no negatives, no NaN, no blanks |
| Normalisation | every column sums to 1e6, max deviation 1.4 × 10⁻⁴ |
| LM22 markers | 29 of 30 present |

The residual 1.4 × 10⁻⁴ is four-significant-figure rounding in the text file,
not an error in the TPM. Before writing, every column sums to exactly 1e6.

**TRBC1 is the single marker miss and is not a pipeline loss** — it is absent
from GENCODE v19 altogether, so it was never in the count matrix. TRBC2, TRAC,
TRDC and TRGC1 are all present.

### Why 55,765 and not 57,820

```
57,820  rows in the deposited count matrix
   -47  ENSGR pseudoautosomal duplicates (they double-count 1,518 reads)
57,773
-2,008  rows sharing an HGNC symbol with another row, summed
55,765
```

No gene was filtered out. CIBERSORT keys on gene symbols and cannot accept
duplicates, so rows mapping to the same symbol must be combined; TPM is
additive within a gene, so they are summed. That step is forced by the tool.

### Genes that are zero in every sample are kept

1,976 symbols are zero across all 474 samples and stay in the file.

They are not biologically absent. **HLA-DPA1 is zero in all 480 deposited
samples and HLA-DPB1 in 474 of 480**, because HTSeq union mode discards reads
that map ambiguously and the MHC class II region is polymorphic enough to
swallow nearly all of them. HLA-DRA survives in only 278 of 480.

This matters twice over. LM22 identifies B cells, dendritic cells and monocytes
partly through class II genes, which is a concrete reason the signature misfits
this dataset — and monocytes here come out at roughly three times their
clinical reference range. And Krebs et al. kept these rows, so removing them
would make this run differ from the one that produced the fractions already
deposited in the GEO metadata, forfeiting a free correctness check.

---

# Step by step

CIBERSORTx runs at <https://cibersortx.stanford.edu> (free academic account) or
as a Docker container. It cannot be run inside a Claude chat.

## Step 1 — upload the mixture

Menu → **Upload Files** → upload `mixture_TPM_hgnc_full_474.txt` as file type
**Mixture**. It must finish uploading before it appears in the dropdown in
step 4.

If the portal rejects it on size, upload `mixture_TPM_hgnc_474.txt` instead, or
switch to the Docker version — the run page itself recommends Docker for large
jobs, and it has no upload limit.

## Step 2 — Analysis Module

Select **`2. Impute Cell Fractions`**.

Module 1 builds a signature matrix from scratch; we are using LM22, which
already exists. Module 3 imputes cell-type-specific expression, which is a
different question from cell proportions.

## Step 3 — Analysis Mode

Select **Custom**.

This is what opens the "Configure Custom Input" panel. The alternative modes
run Stanford's bundled example data.

## Step 4 — the two required files

| Field | Choice |
|:--|:--|
| **Signature matrix file** | `LM22 (22 immune cell types)` |
| **Mixture file** | `mixture_TPM_hgnc_full_474.txt` |

LM22 is the built-in 547-gene, 22-cell-type leukocyte signature. Krebs et al.
and Boltz et al. both used it on this tissue, so choosing it keeps all three
analyses comparable. Do not build a custom signature matrix — a different
signature would make the deposited fractions incomparable to ours.

## Step 5 — batch correction: THIS IS THE SETTING THAT MATTERS

**Run two separate jobs.** They differ only here, and the difference decides
which analysis may use the output.

### Job A — for Arms 2 and 3 and the deconvolution reliability study

| Field | Choice |
|:--|:--|
| Enable batch correction | **CHECKED** |
| Batch correction mode | **B-mode** |
| Optional source GEP file | **`LM22 Source GEP`** |

*Why batch correction:* LM22 was built from microarray profiles of sorted cells;
this data is RNA-seq. Batch correction reduces that platform mismatch.

*Why B-mode and not S-mode:* S-mode is for signature matrices derived from
single-cell RNA-seq. B-mode is for signatures derived from bulk sorted
reference profiles, which is what LM22 is.

*Why the source GEP:* B-mode corrects the mixture against the profiles the
signature was derived from. `LM22 Source GEP` is that file. The field is
labelled optional, but leaving it blank with B-mode selected is a missed step,
not a neutral default.

### Job B — for Arm 1, the prediction arm, only

| Field | Choice |
|:--|:--|
| Enable batch correction | **UNCHECKED** |
| Batch correction mode | n/a |
| Optional source GEP file | leave blank |

*Why:* batch correction is estimated across the entire mixture matrix, so every
sample's fractions end up depending on every other sample's. In cross-validation
that means a held-out sample influenced the features of the training samples.
It is data leakage even though no outcome label is involved, and it would
invalidate the prediction result without producing any visible symptom.

Job B also reproduces Krebs et al.'s method, so its output should closely match
the CIBERSORT fractions already deposited in the GEO metadata. Treat that as a
correctness check: large disagreement means something in the TPM conversion is
wrong.

## Step 6 — Disable quantile normalization

**CHECKED.**

The wording inverts easily: ticking this box turns quantile normalisation
**off**, which is what you want. Stanford's own linked paper recommends
disabling it for RNA-seq, and Krebs et al. disabled it. It is also a
cross-sample operation, so leaving it on would carry the same leakage problem
as batch correction into Job B.

## Step 7 — Run in absolute mode

**UNCHECKED.**

Relative mode returns fractions that sum to 1 per sample. That is what the
compositional analysis needs — the ILR coordinates are built on a closed
composition — and it is the form the deposited Krebs fractions are in, so the
comparison stays direct. Absolute scores can be renormalised back to relative,
but that is an extra conversion with nothing gained.

## Step 8 — Permutations

**100.**

Matches Krebs et al. exactly. Enough resolution to flag samples with p > 0.05,
where the signature does not fit and the fractions should not be trusted. The
run page warns that server load is elevated with a 12-hour execution cap; 1,000
permutations across 474 samples is a real timeout risk for no analytic gain.

## Step 9 — after the run

Download the fractions and save them next to this file under names that record
which job produced them:

```
fractions_bmode_474.csv      <- Job A, batch-corrected
fractions_nobatch_474.csv    <- Job B, no batch correction
```

**Which job produced a file cannot be recovered from its contents.** Using the
batch-corrected set in Arm 1 would invalidate the prediction result silently.

---

# How this compares to the two papers

| | Krebs et al. 2020 | Boltz et al. 2024 | This project |
|:--|:--|:--|:--|
| Tool | classic CIBERSORT | CIBERSORTx | both, deliberately |
| Signature | LM22 | LM22 | LM22 |
| Input | TPM, linear, unfiltered | TPM | TPM, unfiltered, closed to 1e6 |
| Quantile normalisation | disabled | not stated | disabled |
| Batch correction | n/a (classic has none) | **enabled** | off for Arm 1, B-mode for Arms 2–3 |
| Permutations | 100 | not stated | 100 |
| Cell types kept | all 22 | 8 with proportion > 0.02 | 9, non-zero in >50% of samples |
| Cohort | 480 (this dataset) | 1,730 (larger, other data) | 474 (this dataset) |
| Ground truth | none available | CBC for n = 143 | none available |

**Our cell-type rule and Boltz's nearly agree.** Applying Boltz's
"proportion > 0.02" rule to this dataset keeps exactly 8 types; ours keeps 9.
The single difference is activated dendritic cells at mean proportion 0.0049,
which fails their threshold. That is the type already flagged in the project
config as most exposed to LM22 misfit, so the disagreement falls exactly where
expected. Reporting both is cheap and makes the choice defensible.

**Boltz had ground truth and we do not.** They validated against clinic complete
blood counts for 143 subjects: Pearson R² of 0.85 for lymphocytes, 0.76 for
neutrophils, 0.48 for monocytes. Two things follow. Their result is citable
evidence that CIBERSORTx works on whole blood at the lineage level. But
monocytes are the weakest of the three, and monocytes in this dataset sit at
roughly three times their clinical upper bound — the worst-validated cell type
is also the worst-behaved one here. Treat monocyte-containing balances with
corresponding caution.

**Batch correction is the whole reason two jobs exist.** Boltz used it; Krebs
could not. Neither ran a cross-validated prediction model, so neither had to
care about leakage. We do, which turns a free improvement for them into a
disqualifying setting for us.

---

# Prompt for a fresh Claude chat

Paste this and attach the results CSV once a job finishes.

```
I have CIBERSORTx output for a whole-blood RNA-seq cohort: 474 samples,
LM22 signature, relative mode, quantile normalisation disabled, 100
permutations. Attached is the results CSV.

Context:
- Source data is GEO GSE124326 (Krebs et al. 2020), bipolar disorder
  case-control with lithium use recorded. 234 controls, 226 bipolar I,
  14 bipolar II.
- The same dataset already has CIBERSORT fractions in its GEO metadata,
  produced by the original authors with the same method, so my
  no-batch-correction run should closely reproduce them.
- Four LM22 types are exactly zero in all samples (T follicular helper,
  T regulatory, macrophages M1, mast activated). These are structural
  zeros from signature-matrix misfit, not sampling zeros.

Please:
1. Check the fractions sum to 1 per sample using a 1e-4 tolerance, not
   1e-6 (the deposited values miss by up to 5.8e-5 from rounding).
2. Report which cell types are non-zero in more than 50% of samples.
3. Compare against clinical whole-blood reference intervals,
   renormalising to the four comparable classes (neutrophils,
   lymphocytes, monocytes, eosinophils). Expect neutrophils below range
   and monocytes about three times over; that is known LM22 misfit, not
   an error in my run.
4. Report the per-sample P value, correlation and RMSE, and flag any
   sample with P > 0.05 — those are samples where the signature does
   not fit and the fractions should not be trusted.

Do not impute, transform or filter anything. Report what is there.
```
