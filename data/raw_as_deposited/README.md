# raw_as_deposited — the two source files, exactly as GEO serves them

**These files are never edited, filtered or regenerated.** Nothing in this
folder is derived from anything else. If a number anywhere in the project is
questioned, this is where the answer is checked against.

Everything the project needs from the outside world is these two files. There
is no third source, no supplementary table, and no local annotation file.

## What you need

| File | Size | What it is |
|:--|--:|:--|
| `GSE124326_count_matrix.txt.gz` | 16.8 MB | Raw integer HTSeq union-mode counts, **57,820 genes × 480 samples** |
| `GSE124326_series_matrix.txt.gz` | 79.5 KB | Sample metadata, **35 characteristics fields × 480 samples** |
| `CHECKSUMS.sha256` | 193 B | SHA-256 of both files. `build_cohort.py` verifies these before reading. |

## Where they came from

GEO accession **GSE124326** (Krebs et al. 2020, *Psychological Medicine* 50(15),
2575–2586). Use https, not ftp — the ftp endpoint is unreliable.

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124326/suppl/GSE124326_count_matrix.txt.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124326/matrix/GSE124326_series_matrix.txt.gz
```

Verify after downloading:

```bash
sha256sum -c CHECKSUMS.sha256
```

```
84886e8b6a78d99d2705f56daf337b785668bad0670e411a2eb35ae448ce6e0a  GSE124326_count_matrix.txt.gz
21e2d687e00b44f98955ba5545923cc0c89c0f800812f8204f5fa7192b5b9cc0  GSE124326_series_matrix.txt.gz
```

## Four things that will break a naive read

**Columns are sample titles, not GSM accessions.** The count matrix columns
carry a `.counts` suffix and match the series matrix `title` field 480/480
after it is stripped. They currently appear in metadata order, but join on the
key regardless — relying on position would fail silently if GEO re-deposits.

**Every sample title starts with a digit.** R's `make.names` prefixes each with
`X`. Always read with `check.names = FALSE`.

**Characteristics lines are ragged.** A given field does not sit at the same
line index for every sample — `Sex` is present for 478 of 480. Parse each
`field: value` cell by its own key, never by row order.

**47 rows are `ENSGR` pseudoautosomal duplicates**, carrying 1,518 counts in
total. They duplicate `ENSG0*` genes and should be dropped before gene-level
analysis, leaving 57,773. That drop happens downstream, not here.

The 57,820 row count is exactly the GENCODE v19 gene count, which pins the
annotation to **GENCODE v19 / Ensembl 74 / GRCh37 (hg19)**.

## What is not screened here

The deposited `included in final analysis` flag (444 TRUE) screens neither
sequencing depth nor RNA integrity. Nine flag-passing samples fall below 3M
assigned reads and four have RIN below 5. Those screens are the project's own
and are applied downstream.
