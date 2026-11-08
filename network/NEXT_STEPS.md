# What is left, and exactly how to run it

**Updated:** 2026-11-08. Everything below is runnable from a clone of this repo.
Ordered by value per hour. Items 1 and 2 are the science; 3 is bulk compute.

---

## Setup (once, on any machine)

```bash
git clone https://github.com/AliAbdullah12347/Project-CompBio.git
cd Project-CompBio
```

The data is in the repo (`data/`, ~12,368 genes × 474 samples), so nothing else to fetch.

R packages required:

```r
install.packages(c("WGCNA","igraph","fastcluster"))
BiocManager::install(c("limma","edgeR","GO.db","org.Hs.eg.db","preprocessCore","AnnotationDbi","impute"))
```

Verify the environment before running anything:

```bash
Rscript -e 'source("network/00_load.R"); d <- load_project(cache=FALSE); cat("OK:", dim(d$logtpm), "\n")'
bash network/lint_expr.sh          # must print: OK: no experiment uses d$logtpm as expression input
bash network/preflight_queue.sh    # must print: 0 hard faults
```

**Never use `d$logtpm` directly for expression input** — use `expr_for()` or
`expr_samples()`. The lint enforces it. Bypassing it caused a two-week detour.

---

## 1. Raise `perm-031` from K=20 to K=200  ← highest value

**What it is.** The only *calibrated* test of a lithium-specific structural effect. It
draws K random sets of 74 samples from all 474, measures preservation against the
control reference, and asks whether real bp_lith sits below that null.

**Why it must be re-run.** At K=20 the smallest achievable empirical p is 1/20 = 0.05.
BH across 71 modules then cannot clear FDR except at exactly p=0. The result is
suggestive and uninterpretable: bp_lith gold sits below the null in **19 of 20** draws
(p=0.05, BH 0.36 — not significant), and 4 modules fall below *all* 20 draws.

**How to run.**

```bash
# 1. raise K
sed -i '' 's/^K               <- 20L.*/K               <- 200L        # raised 2026-11-08 for BH resolution/' \
  network/experiments/perm-031.R
grep -n '^K ' network/experiments/perm-031.R        # confirm it reads 200L

# 2. raise the cap: 200 draws x ~800s = ~44h, so it MUST checkpoint across sessions
awk -F'\t' -v OFS='\t' 'NR==1{print;next} $1=="perm-031"{$3="PENDING";$4=0;$5="";$6="";$7=0;$9=36000;$8="K raised 20->200"} {print}' \
  network/RUNLOG/queue.tsv > /tmp/q && mv /tmp/q network/RUNLOG/queue.tsv

# 3. run it (resumable - it skips completed draws)
bash network/run_next.sh
```

**Cost.** ~800 s per draw. 180 new draws ≈ **40 hours**, spread over many sessions.
It checkpoints per draw, so interruption is safe.

**What to look for.** With K=200 the minimum empirical p is 0.005, which *can* clear BH
across 71 modules. Check whether `purple`, `blue`, `pink`, `grey60` — the four below all
20 draws now — stay below the null. If they do at K=200, that is a real,
multiple-testing-survivable lithium-specific structural effect, and it would be the
arm's first positive structural result.

---

## 2. A legitimate null for the BPD module enrichment

**What it is.** `analyses/019` found 27 of 69 modules rank-enriched for BPD-associated
genes, against 21 for lithium — even though whole blood yields only **4** BPD DEGs at
BH<0.05 (versus 1,426 for lithium). A size-matched null gives median 0, 95th pct 1, so
it is not a size artifact.

**Why it is not yet established.** Module genes are co-expressed **by construction**, so
their DE statistics are correlated. A Wilcoxon assuming independence is anti-conservative,
and random size-matched gene sets are *uncorrelated*, so that null is too permissive. It
rules out size, not correlation.

**How to run the honest version.** Permute **sample labels**, recompute DE, repeat the
module enrichment. That needs the `de_final` pipeline, not the network arm:

```bash
ls de_final/scripts/          # find the DE entry point
# For each of ~200 permutations: shuffle the BPD label column, rerun the DE fit,
# then re-score modules with network/analyses/019-bpd-module-connection/run.R
# against the permuted statistics. Count significant modules per permutation.
```

**What to look for.** If the observed 27 exceeds the 95th percentile of that
permutation distribution, the distributed sub-threshold BPD signal is real. If not,
27/69 is a co-expression artifact and must be withdrawn.

---

## 3. Why does a group-specific decomposition yield 28 survivors where the fixed reference yields 1?

**What it is.** Newly opened today by `trait-053`. On the fixed verified 71-module
reference, exactly **1 of 71** lithium associations survives composition adjustment.
On `grp-051`'s group-specific 69-module `control_A` network, **28 of 69** do. Model and
sample set are now ruled out — it is the module set.

**How to run.** Re-score `grp-051`'s partition using `trait-053`'s exact model matrix, so
only the partition differs:

```r
# in R, from the repo root
source("network/00_load.R"); d <- load_project(cache=TRUE)
g  <- read.csv("network/runs/grp-051/genes_control_A.csv")
E  <- expr_for(d, group=NULL)
ME <- WGCNA::moduleEigengenes(E, colors=g$module, excludeGrey=TRUE)$eigengenes
# then apply the M_051b model from
# network/experiments/trait-053-composition-reconcile.R verbatim
```

**What to look for.** If 28 survives under the identical model, the group-specific
decomposition is finding composition-independent lithium signal the fixed reference
misses — interesting and worth pursuing. If it collapses toward 1, `grp-051b`'s count
came from something else in its pipeline and the fixed-reference answer stands.

---

## 4. Drain the queue (bulk compute, low value per hour)

19 of 52 rows DONE, 32 PENDING, 0 BLOCKED.

```bash
bash network/run_next.sh --peek    # see the next row without claiming it
bash network/run_next.sh           # run one row
nohup bash network/chain.sh 12 > /tmp/chain.log 2>&1 &   # run rows back-to-back for 12h
```

Most remaining rows are variants (`corr-*`, `net-*`, `pow-*`, `mod-*`) that vary one
network parameter. They matter for robustness, not for the central question — which
items 1–3 address. **Do not expect them to change the conclusions.**

Read results with `medianRank` or against gold **from the same run**. Never compare
`Zsummary` across arms: different references and different permutation budgets put them
on different scales. See `notes/2026-10-26-how-to-read-the-preservation-rows.md`.

---

## Operational notes that cost real time to learn

- **Never log out of the iMac.** eduroam's 802.1X credential is user-scope, so logging
  out locks the keychain and the machine loses all network — every scheduled session
  then dies. Use `pmset displaysleepnow` or `bash network/park.sh`. An Ethernet cable in
  `en0` retires this permanently.
- **Editing a script while it runs:** write a temp file and `mv` it. Bash reads by byte
  offset from an open fd, so in-place truncation corrupts the running process.
- **Generated files** (`INDEX.md`, `catalog.csv`, `RANKING.md`, `ranking.csv`,
  `results_index.csv`) — on a merge conflict, **regenerate**, never merge:
  `Rscript network/build_index.R && Rscript network/build_catalog.R`
- **Find anything:** `bash network/find.sh <term>`, or start at `network/INDEX.md`.
- **Rotate `~/.arm3_token`** — it passed through a chat transcript on 2026-10-14.
