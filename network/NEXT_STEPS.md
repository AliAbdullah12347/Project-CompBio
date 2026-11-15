# What is left, and exactly how to run it

**Updated:** 2026-11-15. Everything below is runnable from a clone of this repo.

> **The queue is drained: 62 of 62 rows DONE, nothing running, nothing open.**
> Items 1-4 below are all ANSWERED and are kept for their method notes, not
> because work remains.
>
> **The one open methodological question is now closed.** The 11 Nov note asked
> for `cent-040`'s control-split null. `cent-056` ran it: two disjoint halves of
> the control group give **43 of 71** "differing" modules against cent-040's 52.
> The procedure was measuring the global level and module size, not biology, and
> cent-040's 52 must never be cited as a lithium effect.
>
> **The arm's negative, resting on every test it has:**
>
> - **0 of 25 parameter variants** (`synth-048`, p 0.441 to 0.989).
> - **0 of 71 modules** against a group-label permutation null at K=157 (`perm-031`).
> - **0 of 5 cell lineages** on bMIND deconvolved expression (`bmind-047`, p 0.761
>   to 0.890). This closes the "bulk dilutes a cell-type-specific effect" objection.
> - **0 of 12360 genes** change intramodular connectivity at BH 5%, K=500 (`hub-037`).
> - **70 of 71 modules** stable across 27 variants; `skyblue2` fails on a Z-bracket
>   edge at 8.64.
>
> **Three positives survive all of it:**
>
> - the **eosinophil module** `mediumpurple2`, whose lithium association
>   strengthens under composition adjustment and which `eos-057` shows is not
>   outlier-driven (LOO 0/226, bootstrap CI excludes zero, fragility 13 subjects);
> - a **global 11-22% thinning** of the lithium co-expression graph, reproduced in
>   three independent subsamples and localising to no module;
> - **centrality separates diagnosis** (`netmet-060`/`netmet-062`). Cross-diagnosis
>   pairs agree less than within-kind pairs on all eight gene-level metrics, 10 of
>   10 draws, before and after residualising age, sex, tobacco, RIN, plate and the
>   sequencing PCs -- but the gap shrinks by 33-64%, so roughly half the raw effect
>   was demographic and technical. Lithium behaves as a within-kind pair throughout.
>   See `findings/2026-11-14-centrality-separates-diagnosis-and-half-of-it-is-covariates.md`.
>
> **The permutation-null pair (`perm-054` / `perm-059`) is the other late result.**
> Lithium has 1,426 DEGs and **no** module structure (21 of 69 vs null 16.4,
> p = 0.194); BPD has 4 DEGs and **some** (27 of 69 vs 16.5, p = 0.018).
> `analyses/019`'s "21 modules enriched for lithium" is **withdrawn**, and its 27
> for BPD survives against 16.5 rather than the ~0 the old null implied.
>
> **Seven rows reported a wrong number while exiting 0.** All are fixed and re-run;
> the detail is in `network/cluster/README.md`. `perm-032`'s "8 modules significant"
> is an artifact of p=0 -- use `results_permutation_test_addone.csv`, which gives 0.
>
> **A fetch can silently overwrite a committed result.** On 14 Nov a routine
> `collect.py --fetch` replaced `cent-056`'s 20-split null with a 1-split one:
> K 20 to 1, sd to NA, p 0.048 to 0.5, and it still looked like a result. Every
> sliced row was one fetch away from the same thing. `collect.py` now refuses to
> overwrite a git-tracked, unmodified file whose content would change.
>
> **The bulk compute is on the research cluster, not the iMac.** See
> `network/cluster/README.md`. **Do not read an upload count as a success count**
> - the first drain reported 30 rows uploaded and had finished 2.
> Return codes are in `status/*.tsv`.

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

## 3. Why does a group-specific decomposition yield 28 survivors where the fixed reference yields 1?  -- ANSWERED, 2026-11-11

**It does not.** `recon-055` re-scored `grp-051`'s `control_A` partition
through `trait-053`'s model, changing nothing but the partition, and got
**1 of 69** under both adjusted models -- matching the fixed reference's
1 of 71. `grp-051b`'s own `module_trait_joint.csv` agrees: `li_padj` = 29 is
lithium adjusted for diagnosis and technical covariates but NOT composition;
`li_adj_padj` = 1 is the same model with the ILR balances added. The 28 in
this file was the unadjusted column.

The survivors are `mediumpurple2` in the fixed reference and `orangered3` in
`grp-051` -- the same genes, Jaccard 0.846 against a 0.184 background. Two
decompositions that agree on little else each find one module, and it is the
same module.

The partition does matter under the weaker `li + b1..b4` model, where the two
give 26 and 33 on the full sample set, but that gap closes entirely once
diagnosis and the technical covariates enter.

Written up in
`findings/2026-11-11-the-28-vs-1-discrepancy-was-a-misread-column.md`.

---

## 4. Drain the queue (bulk compute, low value per hour)

21 of 52 rows DONE, 29 RUNNING on the cluster, 0 PENDING, 0 BLOCKED.

```bash
export KUBECONFIG=./kubeconfig-user-aabdullah.yaml
kubectl get job.batch.volcano.sh                   # arm3-run-01..07, arm3-perm031
python network/cluster/collect.py                  # what landed
python network/cluster/collect.py --fetch          # download and unpack
```

On a single machine instead:

```bash
bash network/run_next.sh --peek    # see the next row without claiming it
bash network/run_next.sh           # run one row
nohup bash network/chain.sh 12 > /tmp/chain.log 2>&1 &   # run rows back-to-back for 12h
```

Two rows were outstanding for reasons the queue could not show:

- **`ann-044` had no script.** The drain only deals rows that have a
  `network/experiments/<id>.R`, so it was skipped silently at every attempt.
  Written 9 Nov. It is the row that re-annotates the *current* 71-module
  decomposition — `analyses/012` annotated the retracted 33-module one, so
  every finding in it, `darkturquoise`=ribosomal included, describes modules
  that no longer exist.
- **`net-007` was left `RUNNING`** by the iMac on 8 Nov with no output, and the
  drain filtered on `PENDING`, so nothing ever picked it up.

**`synth-048` must be re-run** once these land. It completed on 8 Nov, but from
`6/25` robustness rows — its "71/72 modules stable across ≥80% of variants"
currently describes six variants.

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
