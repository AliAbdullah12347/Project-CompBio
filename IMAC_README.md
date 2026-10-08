# IMAC_README — orientation for the unattended network-analysis run

You are running on a dedicated iMac, for several days, with **no human available**.
Read this file completely before doing anything else. It tells you where everything
is, what is deliberately missing, and what to do when you cannot find something.

If you only read one line: **`source("network/00_load.R"); d <- load_project()`**
gives you every matrix this arm needs, already checked against the verified facts.

---

## 1. What this project is

Almost every bipolar patient in a public dataset is medicated, so a "bipolar
signature" in blood may just be a **lithium** signature. Lithium raises neutrophil
proportion, which shifts bulk expression without any cell changing behaviour. This
project separates those two things in GEO **GSE124326**.

Three arms. **Arm 3, networks, is yours.** Arms 1 and 2 are other people's work;
their results are here for you to use, not to redo.

The network question: **do gene co-expression networks differ across
healthy / bipolar-off-lithium / bipolar-on-lithium, and does any difference
survive removing cell composition?**

Arm 3 is **exploratory**. You are not required to produce a significant finding.
A well-documented negative or inconclusive result is a complete deliverable.

---

## 2. The data, in one function

```r
source("network/00_load.R")
d <- load_project()        # ~40 s first call, cached to network/cache/ after
```

| field | what it is |
|:--|:--|
| `d$logtpm` | 12,368 × 474 numeric, log2(TPM + 1), **genes × samples** |
| `d$counts` | 12,368 × 474 integer, filtered raw counts |
| `d$meta` | 474 rows: `age sex rin plate seqpc1-3 tobacco dx lithium group` |
| `d$groups` | factor: `control` (234) / `bp_nolith` (74) / `bp_lith` (152) |
| `d$ilr` | 474 × 4 ILR balances `b1..b4` |
| `d$frac` | 474 × 5 lineage fractions: `gran mono T NK B` |
| `d$tmm` | library sizes and TMM factors |
| `d$tobacco` | 20 multiple imputations of tobacco status |
| `d$sym` | Ensembl id → HGNC symbol, **all 12,368 map** |
| `d$lm22` | 547 LM22 signature gene symbols (for the circularity check) |

Helper, in the same file:

```r
E <- expr_for(d, group = "control")                       # samples x genes, WGCNA orientation
E <- expr_for(d, group = "bp_lith",
              residualise = ~ age + sex + rin + b1 + b2 + b3 + b4)   # composition removed
```

`load_project()` **stops** rather than returning a wrong object if any shape or
group count is off. If it stops, do not work around it — see §7.

---

## 3. Where everything lives

```
Implementation/                 <- repo root, and the only directory you may touch
├── IMAC_README.md              <- this file
├── START_HERE.md               <- general orientation (written for a laptop session)
├── CLAUDE.md  (one level up)   <- project rules and the VERIFIED FACTS table
├── network/                    <- YOUR working directory. Create what you need here.
│   ├── 00_load.R               <- the loader above. Do not rewrite it; extend it.
│   └── cache/                  <- gitignored scratch
├── data/
│   ├── analysis_matrices/      <- the 14 CSVs the loader reads. Canonical.
│   ├── cohort_474/             <- raw counts, metadata, CIBERSORTx fractions + inputs
│   └── LM22.txt                <- the 547-gene signature matrix
├── de_v2/                      <- CURRENT differential expression results
│   ├── DECISIONS.md            <- 16 sections, every DE decision justified
│   ├── results/                <- per-gene tables, sweeps, corrections
│   └── scripts/mtc.R           <- hand-written multiple-testing library, REUSE THIS
├── de_final/results/           <- the headline DE tables (4 contrasts)
├── mediation/results/          <- Arm 2 output: ILR balances, mediated proportions
├── de_analysis/                <- v1 DE. METHODS.md (21 sections) is still the
│                                  reference for the gene filter and detection floors
├── Boltz Paper/                <- recreation of Boltz et al. 2024
├── Krebs Paper/                <- recreation of Krebs et al. 2020
├── index/build_index.py        <- regenerates INDEX.db / INDEX.csv
└── INDEX.csv, INDEX.db         <- searchable description of every file in the repo
```

**Navigate with the index, not by guessing:**

```bash
python index/build_index.py --query "SELECT path, description FROM files WHERE path LIKE '%sweep%'"
```

Views: `v_current_outputs`, `v_code`, `v_do_not_use`, `v_biggest`, `v_undescribed`.
If `DIRECTORY.txt` and `INDEX.db` disagree, **the index is right** — it is generated.

---

## 4. What is deliberately NOT in this repo

Nothing here is an accident. Do not hunt for these; regenerate or skip.

| missing | why | what to do |
|:--|:--|:--|
| `de_analysis/data/prep.rds` | 74 MB binary cache | **Not needed.** `00_load.R` rebuilds the equivalent from tracked CSVs. |
| `*.rds` generally | gitignored: unreadable and regenerable | regenerate, or work from the CSVs |
| `data/cohort_474/mixture_TPM_hgnc_full_474.txt` | 102 MB, over GitHub's hard limit | `python data/make_mixture_full.py`, verify against `mixture_full_474.sha256` |
| bMIND cell-type expression profiles | ~3.1 h to compute, binary | `Rscript de_analysis/scripts/06_bmind_profiles.R`. Only needed for cell-type-specific networks. **Budget a whole session for it.** |
| journal PDFs, publisher supplementary | not ours to redistribute | cite, do not seek |
| GEO / GENCODE / DICE downloads | large public files | `download_data.sh` in each paper folder |

---

## 5. Verified facts your work must reproduce

From `CLAUDE.md` §2. These were checked against the actual files. **If a computation
disagrees with one of these, the computation is wrong until proven otherwise. Do not
silently adjust a number to match — investigate, and write down what you found.**

- 12,368 genes after filtering (>10 counts in ≥90% of samples); 474 subjects.
- Groups: **234 healthy controls, 74 BP1 off lithium, 152 BP1 on lithium.**
- The lithium contrast is **BP1 only**. No BP2 subject takes lithium.
- `assessment group` is **constant within BP1** — including it makes the design
  matrix singular. Drop it for any BP1-only model.
- QC exclusion is **not random with respect to diagnosis** (35 of 240 controls
  dropped, Fisher p = 3.04e-10). Every case/control comparison runs on a filtered
  control group. Flag this wherever one appears.
- Tobacco is **balanced with lithium inside BP1** (35.1% vs 31.6%, p = 0.70). It is
  a confounder for case/control, not for the lithium contrast.
- LM22 misfits whole blood: monocytes come out ~3× over reference range, and four
  LM22 types are **structural zeros** in all 480 samples. That is a known
  preliminary result, not a bug.

---

## 6. Traps that have already cost this project time

1. **No cell-type zero is evidence of absence.** Every lineage's detection floor
   exceeds the 0.5 log2FC meaningfulness threshold; B cells would need a
   72,000,000-fold change. `de_v2/DECISIONS.md` §11.
2. **Modules must be defined ONCE**, on a fixed reference group, then measured in
   every subsample. If you re-detect modules per draw, module 7 in draw 1 is not
   module 7 in draw 2 and nothing can be pooled across draws. This is the single
   most important design point in the arm.
3. **Zsummary has no calibrated null.** Langfelder et al. (2011) give rule-of-thumb
   cutoffs (<2 not preserved, 2–10 weak, >10 strong), not p-values. Do not push it
   through a normal CDF. **medianRank** is relative to one run only and has no null
   at all. Neither can be multiple-testing corrected.
4. **Subsample draws are not independent.** Two draws of 74 from the 152 lithium
   users share ~36 people on average. Fisher's and Stouffer's methods and averaging
   p-values are all invalid here. BP1-off-lithium is **never** subsampled — it
   already is exactly 74.
5. **Ensembl ids in `data/cohort_474/counts_474.tsv.gz` carry version suffixes**
   (`ENSG00000000003.10`); everything else in the repo is unversioned. Strip with
   `strsplit(x, ".", fixed = TRUE)` — a regex needs `"\\."`, and getting this wrong
   silently produced empty intersections for 45 minutes once. `00_load.R` already
   gives you unversioned ids.
6. **Don't pool all comparisons into one BH.** It *raised* the bipolar DE counts
   (4→16, 0→22) because lithium's strong signal lifts the threshold for everyone.
   Use Benjamini–Bogomolov. `de_v2/DECISIONS.md` §8.
7. **"Off lithium" is not "unmedicated."** The deposit records lithium and nothing
   else. Every bipolar claim inherits this limitation.
8. **`qvalue` and `IHW` are not installed and must not be added.** `de_v2/scripts/mtc.R`
   is a hand-written, validated replacement (18 checks in `00_validate_mtc.R`, all
   passing). Reuse it.
9. **WGCNA needs `enableWGCNAThreads()`** or it runs single-threaded with a warning.
   And on small gene sets `pickSoftThreshold` often returns `powerEstimate = NA` —
   scan a wider `powerVector` (e.g. `1:20`) and fall back to the WGCNA defaults by
   sample size if it is still NA.

---

## 7. When you cannot find something, or something fails

Work through this in order. **Never fabricate a number, and never silently skip a
check.** An unattended run that quietly produces wrong output is worse than one that
stops and says why.

1. **Search the index first.**
   `python index/build_index.py --query "SELECT path, description FROM files WHERE description LIKE '%<term>%'"`
2. **Then grep the docs:** `de_v2/DECISIONS.md`, `de_analysis/METHODS.md`,
   `mediation/METHODS.md`, `START_HERE.md`, `CLAUDE.md`.
3. **Then check §4 above** — it may be deliberately absent with a stated rebuild command.
4. **If it is regenerable, regenerate it**, and log how long it took.
5. **If it is genuinely missing**, write the gap into `network/BLOCKED.md` with:
   what you needed, what you tried, what you did instead. Then **move on to a
   different analysis**. Do not stall the run waiting for a human.
6. **If a verified fact in §5 fails to reproduce**, stop that analysis, record the
   discrepancy in `network/BLOCKED.md` in full detail, and start a different one.
   Do not adjust your pipeline to force the expected number.

---

## 8. Hard rules — violating any of these invalidates the work

- **Everything stays inside `Implementation/`.** Nothing outside it is read or
  written. This holds even with bypass permissions on.
- **Ask before adding a dependency.** On this unattended run that means: get the
  full install list approved in the first twenty minutes, then add nothing further.
  If you later need a package that was not approved, record it in
  `network/BLOCKED.md` and find a way without it.
- **Verify rather than assume.** When a number can be computed from the data,
  compute it and print it. Never carry a figure forward from a paper or from these
  notes without checking.
- **Say when something looks wrong.** If a result contradicts §5, or a method cannot
  be made sound, write that down rather than working around it.
- **Comment the why, not the what.**
- **Pre-specify before you look.** Any threshold, margin or model choice gets written
  into the run's config *before* the result is computed. Changing it afterwards, or
  reporting only the variant that worked, is the one thing that would invalidate the
  whole arm.

---

## 9. Environment

R with: `WGCNA` (plus `impute`, `preprocessCore`, `GO.db`, `AnnotationDbi`,
`dynamicTreeCut`, `fastcluster`), `limma`, `edgeR`. Python 3 with `sqlite3`
(stdlib) for the index.

Check what you have before planning:

```r
for (p in c("WGCNA","limma","edgeR","impute","preprocessCore","dynamicTreeCut","fastcluster"))
  cat(sprintf("%-16s %s\n", p, ifelse(requireNamespace(p, quietly=TRUE), "yes", "MISSING")))
```

Anything missing goes in the single approval request at the start of the run.

---

## 10. Repository conventions

- Results tables are tracked in git; `.rds` caches and large downloads are not.
- Every results folder gets a `METHODS.md` recording decisions and their reasons.
- `python index/build_index.py` regenerates the file index — run it after adding
  files so the next session can find them.
- Commit messages: imperative, one line of what and one of why.
