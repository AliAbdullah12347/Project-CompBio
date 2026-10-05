# CIBERSORTx high-resolution mode — cell-type expression without single-cell data

Companion to `CIBERSORTX.md`, which covers the fractions run we already did.
This covers **Module 3, high-resolution mode**: per-sample, per-gene,
per-cell-type expression estimates — the same kind of output bMIND gives.

Everything marked **verified** below is quoted from Newman et al. 2019
(*Nat Biotechnol* 37:773–782, `Papers- Annotated Bibliography/CIBERSORTx.pdf`).

---

## 1. Yes, this works without single-cell data

Boltz et al. give as their reason for choosing bMIND that CIBERSORTx and
related methods "require a single-cell RNA-seq reference dataset, ideally in a
matched subset of individuals, for all the cell types of interest across all
genes of interest."

**That is not what the method requires.** Newman et al. — Boltz's own reference
11 — state the input to high-resolution purification is a matrix of bulk
profiles plus a signature matrix, and demonstrate it with signatures that
contain no single-cell data at all:

> "To enumerate FL immune proportions, we applied **LM22, a microarray-derived
> signature matrix** for distinguishing 22 human hematopoietic cell subsets in
> bulk tissues"

LM22 is then used for high-resolution purification of B cells (their Fig. 5a is
labelled "high-resolution purification of B cells"). The lung-cancer analysis
uses a signature built from **bulk FACS-sorted RNA-seq**, not single cells, and
is applied to 1,132 TCGA profiles. On bulk versus single-cell signatures they
report:

> "While significant differences in performance between reference signatures
> derived from bulk populations and those derived from scRNA-seq data were
> **not observed**"

Be fair about the caveat: Boltz's sentence says "ideally in a matched subset of
individuals", which admits a charitable reading — that *good* cell-type
expression wants a matched reference they did not have, rather than that one is
strictly required. But the verb is "require", and Newman's own discussion says
bulk-derived signatures performed no worse. Their second reason for bMIND
(computational efficiency, higher average R² against snRNA-seq) does not depend
on this and stands on its own.

---

## 2. What the paper says we need — checked against what we have

| requirement (verified) | our position |
|:--|:--|
| Bulk mixture GEPs in **non-log linear space** — "expression data in M′ and B are represented in non-log linear space" | ✅ `mixture_TPM_hgnc_full_474.txt` is TPM, columns summing to 1e6 |
| A signature matrix | ✅ LM22, built into the portal, the same one our fractions used |
| System overdetermined, **k > c** (samples > cell types) | ✅ 474 > 22 |
| "The largest gains were achieved when analyzing **at least four–fivefold more mixture samples than cell types**" | ✅ 474 / 22 = **21.5×**, far past the recommendation |
| Gene symbols as row IDs | ✅ HGNC symbols |

On cohort size we are in unusually good shape. Newman's own accuracy curve
(Fig. 3e) is built from 302 follicular lymphoma samples; we have 474.

**Two assumptions you are taking on** (verified, their Methods):

> "First, we assume that **each gene can be analyzed independently**. Although
> ignoring gene–gene covariance relationships will probably impact the resolving
> power for some genes…"

> "Second, for a given gene, we assume that **at least some evidence of
> cell-type-specific differential expression is detectable in bulk tissue
> samples**, even if it is statistically insignificant."

The second matters for us: a cell-type effect that is *perfectly* cancelled in
bulk is invisible to this method by construction. That is the same masking
hypothesis we tested and did not support (`de_analysis/METHODS.md` §19), so it
is consistent, but it is an assumption, not a result.

---

## 3. The file to upload

**`Implementation/data/cohort_474/mixture_TPM_hgnc_full_474.txt`** — 107 MB,
55,765 gene symbols × 474 samples, tab-delimited, first column header
`GeneSymbol`.

Use **this exact file and no other**. It is the one that produced
`fractions_bmode_474_CIBERSORTx.csv`. Module 3 re-derives fractions internally
as part of imputation, so any difference in the mixture file or the settings
below produces cell-type profiles that are not consistent with the fractions
already in the analysis.

Do not substitute `mixture_TPM_hgnc_474.txt` (the smaller fallback) or either
file in `cibersortx/data/` (those are the 480- and 431-sample versions).

---

## 4. Portal steps

Settings that must match the fractions run are marked **→ must match**. The
exact field labels in Module 3 differ slightly from Module 2; where I give a
label, confirm it on screen rather than clicking by position.

### Step 1 — sign in
`cibersortx.stanford.edu`. Same account as the fractions run.

### Step 2 — module
Select **`3. Impute Cell Expression`** (not `2. Impute Cell Fractions`, which is
what we ran before).

### Step 3 — analysis mode
Select **Custom**. This opens the configuration panel, as in Module 2.

### Step 4 — the core choice: **High Resolution**, not Group Mode

| option | what it gives | use it? |
|:--|:--|:--|
| **Group mode** | ONE representative profile per cell type for the whole cohort | ✗ no per-sample values, so no differential expression |
| **High Resolution** | per-gene, per-cell-type, **per-sample** | ✅ **this one** |

This is the setting that makes the output comparable to bMIND. Group mode would
give 22 columns total and would be useless for our contrasts.

### Step 5 — files

| field | value |
|:--|:--|
| Mixture file | upload `mixture_TPM_hgnc_full_474.txt` |
| Signature matrix file | **`LM22 (22 immune cell types)`** → must match |

Do not build a custom signature matrix. A different signature makes the result
incomparable to our fractions and to Krebs et al.

### Step 6 — batch correction → must match

| field | value |
|:--|:--|
| Enable batch correction | **CHECKED** |
| Batch correction mode | **B-mode** |
| Source GEP file | **`LM22 Source GEP`** |

B-mode because LM22 is derived from bulk sorted populations, not single cells
(S-mode is the single-cell case). The source GEP field is labelled optional;
leaving it blank with B-mode selected is a missed step, not a valid choice.

### Step 7 — quantile normalisation → must match

**Disable quantile normalisation** (tick the box labelled *Disable quantile
normalization*; the wording inverts easily — ticking it turns QN **off**).

Off is correct for RNA-seq, and Krebs et al. disabled it.

### Step 8 — cell types to impute

If the form asks which of the 22 subsets to impute, select **all of them**, then
discard afterwards rather than choosing now. Deciding which cell types are
interesting before seeing the output is a selection we would have to account
for later.

**Four LM22 types are zero in all 480 samples** — T cells follicular helper,
T cells regulatory (Tregs), Macrophages M1, Mast cells activated. A cell type
with zero abundance everywhere carries no information for the imputation, so
expect those columns to come back empty or all-NA. That is correct behaviour,
not a failed run.

### Step 9 — permutations

Set to **0 / default** for this module if the field appears. Permutations price
the *fractions* significance test; they do nothing for expression imputation and
only add runtime. (Our fractions run used 100, which is already recorded.)

### Step 10 — run

Expect this to take **substantially longer than the fractions job**. High
resolution solves a separate constrained problem per gene per cell type:
roughly 55,765 × 22 systems across 474 samples.

---

## 5. If the web portal will not take it

474 samples × 55,765 genes is a large high-resolution job and the web portal may
time out or refuse the upload. The supported route for jobs this size is the
**CIBERSORTx Docker container**, which runs the identical algorithm locally.

Request a token from the website (Menu → Download), then run roughly:

```bash
docker run -v "$PWD/in":/src/data -v "$PWD/out":/src/outdir \
  cibersortx/hires \
  --username YOUR_EMAIL --token YOUR_TOKEN \
  --mixture mixture_TPM_hgnc_full_474.txt \
  --sigmatrix LM22.txt \
  --classes LM22_source_GEP.txt \
  --QN FALSE
```

Confirm the image name, the flag spellings and which files the container expects
against the site's own download page before relying on this — I have verified
the *method* from the paper, not the current container interface. `--QN FALSE`
is the one flag I would check most carefully, since it is the setting whose
wording inverts.

**Fallback if it is still too slow:** restrict the mixture to the 12,368 genes
that survive our filter **plus every LM22 signature gene** (the signature genes
must stay or the fraction step breaks). Download `LM22.txt` from the portal to
get that gene list, then subset `mixture_TPM_hgnc_full_474.txt` to the union.
Record that you did this — it changes the gene universe, and our percentages are
only comparable over a common universe (`de_analysis/METHODS.md` §3).

---

## 6. What comes back, and the thing to watch for

High-resolution mode returns **one file per cell type**, genes × 474 samples.

**Expect a lot of NA, and expect it to be uneven.** CIBERSORTx applies an
adaptive noise filter that removes genes it cannot estimate reliably:

> "we developed an **adaptive noise filter to eliminate unreliably estimated
> genes** for each cell type"

> "the fraction of recovered genes after adaptive filtration was **proportional
> to both the number of evaluated samples and the proportion of each cell
> subset**"

That second sentence is the one to take seriously. Recovery scales with cell-type
abundance, so neutrophils will come back nearly complete and the rare subsets
will come back mostly empty — the same abundance-driven limit we measured for
bMIND, where the B-cell detection floor was 26.1 log2FC
(`de_analysis/METHODS.md` §17).

**This is a point in CIBERSORTx's favour, not against it.** bMIND returns a
shrunken number for every gene in every lineage whether or not the data support
one — B-cell estimates had ~100× less variance than granulocyte estimates, which
is shrinkage, not biology. CIBERSORTx instead declines to answer. An explicit NA
is more honest than a confidently shrunk value, and it means the "how much of
this is real" question is answered by the method rather than left to us.

### Comparing it to bMIND

Our bMIND results are at 5 aggregated lineages; CIBERSORTx gives 22 LM22 types.
To compare, aggregate the CIBERSORTx output using the lineage map in
`de_analysis/config.yaml` — but note that **expression profiles cannot be summed
the way fractions were**. The bulk contribution of a lineage is the
abundance-weighted average of its members:

> x_lineage = Σ(p_type × x_type) / Σ(p_type)

using each sample's own LM22 fractions as weights. Summing unweighted would
treat a 0.1% subset as equal to a 38% one.

Before running any contrast on the result, check per cell type: how many genes
survived the filter, and how that number tracks abundance. If the rare lineages
come back with only a few hundred genes, the comparison with bMIND is a
comparison of what each method *refuses* to estimate, which is itself the
interesting finding.

---

## 7. Caveats that carry over

- **Circularity is unchanged.** The fractions are estimated from the same
  expression matrix the imputation then decomposes. For the 547 LM22 marker
  genes this is close to tautological. `CLAUDE.md` §4.2 requires excluding
  marker genes from any mediation outcome set, and that applies here too.
- **B-mode is transductive.** It corrects across the whole cohort at once, so
  this output is fine for differential expression (no train/test split) and
  **disqualified from Arm 1's cross-validated prediction**, exactly as recorded
  for the fractions.
- **22 types is finer than we trust.** We aggregated to 5 lineages because
  between-method ICC at 22-type resolution was 0.389 versus 0.693 at lineage
  level (`de_analysis/METHODS.md` §5). Nothing about using CIBERSORTx for
  expression changes that; treat the 22 per-type profiles as inputs to the
  lineage aggregation, not as 22 independent results.
- **A second method is a comparison, not a confirmation.** Both methods take the
  same fractions from the same matrix. Agreement between them is weaker evidence
  than it looks, because the shared input is the part most likely to be wrong.

---

## 8. The web portal caps the gene subset at 1,000 genes

Discovered on the form, not in the paper. The Gene Subset file is **required**
for high-resolution mode, and the portal states:

> "For analysis on whole-transcriptomes, a gene subset file with a list of genes
> of interest is required (**1000 genes or less**), to reduce the runtime of
> CIBERSORTx, since this type of job is computationally intensive. For a full
> transcriptome version, please navigate to the Downloads page to request
> access."

So the whole-transcriptome run is **Docker-only**. Three subset files are
prepared in this folder:

| file | genes | purpose |
|:--|--:|:--|
| `gene_subset_pilot500_474.txt` | 500 | Dry run. 50 genes per expression decile, so its CV distribution resembles the full transcriptome's rather than one corner of it. |
| `gene_subset_LItop1000_474.txt` | 1,000 | A real analysis that fits the cap: the 1,000 strongest whole-blood lithium genes (all BH q < 0.05, range 7.0e-16 to 3.1e-02). Asks which lineage carries the known effect. |
| `gene_subset_12k_474.txt` | 12,363 | The full analysis gene set. **Docker only.** |

### Do NOT split the 12,363 into 13 jobs of 1,000

This looks like an obvious workaround and it is not valid. The paper says the
filter is

> "an adaptive noise filter based on the **transcriptome-wide distribution of
> coefficient of variation (c.v.) values for each cell type**"

The keep/drop threshold is therefore derived from the spread of CV values across
**all genes present in that run**. Thirteen batches of 1,000 would each compute
their own threshold from their own 1,000 genes, so whether a gene survives would
depend on which batch it landed in. The per-gene imputation is independent — the
paper assumes exactly that — but the filtering on top of it is not.

This is the same class of defect as the bMIND chunking problem recorded in
`de_analysis/METHODS.md` §13, where posterior means were clipped to the range of
each 400-gene chunk rather than the whole matrix. Worth stating in the write-up:
two different deconvolution tools, the same failure mode, both only visible by
reading what the method does across genes rather than within one.

### The selection caveat on the lithium subset

`gene_subset_LItop1000_474.txt` is chosen using the whole-blood lithium result.
That makes it valid for **localisation** — these genes move in bulk, which
lineage are they moving in? — and invalid for **discovery**. Do not compute an
FDR over those 1,000 as though they were a fresh sample; they were selected on a
related statistic from the same subjects. Report them as "where the known effect
sits", which is the question `de_analysis/scripts/12_celltype_structure.R`
already answers for bMIND, and which this run would answer by an independent
method.

### What the pilot will and will not tell you

Will: that the settings are right, what the output files look like, the per-gene
runtime, and roughly how recovery scales with cell-type abundance.

Will not: the final NA rates. Because the threshold is transcriptome-wide, a
500-gene run derives it from 500 genes. The decile stratification makes that
closer to the full-transcriptome answer than a random or top-expressed sample
would be, but it is not the same number.

---

## 9. Run failure: `$ operator is invalid for atomic vectors`

Seen on 2 October 2026. Full message:

```
Error in mat$cv : $ operator is invalid for atomic vectors
Calls: CIBERSORTxHiRes -> CIBERSORTxGEP -> runCIBERSORTxGEP
In addition: Warning message:
In mclapply(1:no_cores, res, mc.cores = no_cores, mc.set.seed = FALSE, :
  all scheduled cores encountered errors in user code
Execution halted
```

### Read it bottom-up

The reported error is **not** the failure. `mclapply` catches errors inside its
workers instead of crashing, and returns them as ordinary objects. So:

1. `all scheduled cores encountered errors in user code` — **every** parallel
   worker failed. This is the real event, and its message is swallowed.
2. The results were combined into a plain vector of error objects rather than
   the expected table.
3. `mat$cv` then asked that vector for a `cv` field. `$` does not work on an
   atomic vector, so R reports that instead.

`cv` is almost certainly the coefficient of variation used by the adaptive noise
filter — the paper describes it as working from "the transcriptome-wide
distribution of coefficient of variation (c.v.) values for each cell type". So
execution reached the filtering stage and found nothing to filter.

### Ruled out: degenerate input rows

Checked every gene in all three subset files against
`mixture_TPM_hgnc_full_474.txt`:

| file | found in mixture | all-zero rows | zero-variance rows | expressed in <5% of samples |
|:--|--:|--:|--:|--:|
| pilot500 | 500 / 500 | 0 | 0 | 0 |
| 12k | 12,363 / 12,363 | 0 | 0 | 0 |

A gene that is constant or absent would give an undefined coefficient of
variation and is the obvious way to kill the filter. There are none, so that is
not the cause.

### Prime suspect: the gene subset file format

This was flagged as unverified when the files were written. Header-bearing
variants now exist beside the originals:

```
gene_subset_pilot500_474_hdr.txt
gene_subset_LItop1000_474_hdr.txt
gene_subset_12k_474_hdr.txt
```

Each is the same list with a single line `GeneSymbol` prepended, matching the
first-column header of the mixture file. Try the `_hdr` version of whichever
subset failed. If the parser expects a header and does not get one, the first
gene is consumed as the column name, which can leave the downstream structure
malformed in exactly the way this traceback suggests.

### If running under Docker: surface the real error

The hidden worker message is the thing worth having. Force serial execution so
`mclapply` cannot swallow it — check the container's own options first rather
than guessing flag names:

```bash
docker run --rm cibersortx/hires --help
```

Then re-run with one core, and with the 500-gene pilot rather than the full
set. Note also that Docker Desktop on Windows gets only a slice of host RAM
(often 2-4 GB under WSL2 unless `.wslconfig` says otherwise), and this machine
has 7.7 GB total. Every worker dying at once is the classic out-of-memory
signature, and this project already hit exactly that with bMIND at six workers
(`de_analysis/METHODS.md` §20).

### If running on the web portal

Memory and core count are not yours to control, so the input files are the only
lever. Resubmit the `_hdr` variant, and use `gene_subset_pilot500_474_hdr.txt`
first — a 500-gene job fails fast and costs nothing to retry.

### RESOLVED — the cause was Windows line endings, not the header

The "prime suspect" above was wrong. The header made no difference because the
header was never the problem.

Every gene subset file had been written with **CRLF line endings**, because R's
`writeLines` on Windows opens the connection in text mode and silently
translates `\n` to `\r\n`. The mixture file, built by a different script, had
plain LF.

```
subset file :  A  B  H  D  5  \r  \n      <- trailing carriage return
mixture file:  5  S  _  r  R  N  A  \n    <- clean
```

Inside the Linux container every symbol therefore carried an invisible trailing
`\r`, so `"ABHD5\r"` never matched `"ABHD5"`. **Zero genes matched, in every
chunk.** That is why *all* workers failed rather than some: the failure was in
the input every worker shared, not in any particular gene.

All six files were converted to LF and verified byte-for-byte:

| file | CR bytes | lines | symbols matching the mixture |
|:--|--:|--:|--:|
| `gene_subset_pilot500_474.txt` | 0 | 500 | 500 / 500 |
| `gene_subset_LItop1000_474.txt` | 0 | 1,000 | 1,000 / 1,000 |
| `gene_subset_12k_474.txt` | 0 | 12,363 | 12,363 / 12,363 |
| …and the three `_hdr` variants | 0 | +1 each | all matching |

**Use the plain (non-`_hdr`) files.** The header variants exist only because the
format was briefly suspected; there is no evidence a header is wanted, and a
file whose first data row is the word `GeneSymbol` is a worse guess than one
without it.

### The general lesson, worth carrying

Any file written on Windows and read by a Linux container needs its line
endings checked. The failure mode is silent — the data looks right in an editor,
every symbol is spelled correctly, and nothing matches. Verify with bytes, not
with eyes:

```bash
python -c "b=open('FILE','rb').read(); print('CR bytes:', b.count(b'\r'))"
```

Two further notes on how this was diagnosed, both of which cost time:

* The reported R error named the wrong thing entirely. `mat$cv` is three steps
  downstream of the actual failure, because `mclapply` catches worker errors
  instead of propagating them. Always read an `mclapply` traceback bottom-up.
* A first check for carriage returns used `od -c | grep '\r'` and returned a
  false positive on two files that were already fixed. Counting the bytes
  directly settled it. When a diagnostic and a functional test disagree, trust
  the functional test — here, "do the symbols actually match?".
