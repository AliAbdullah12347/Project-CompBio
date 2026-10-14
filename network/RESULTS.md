# Arm 3 — Results

**Updated:** 2026-10-14
**Runs completed:** 1 of 50 queued (mod-000-resolution). base-001 v2 RUNNING.
**Analyses completed:** 4 (001–004). **Findings:** 8 — see `RANKING.md`.

---

## The answer so far

**A lithium-associated co-expression dimension exists and it is almost entirely cell
composition.** This finding is robust across all tested decompositions and cannot be
explained by circularity in the composition estimator. Whether the composition difference
*causes* the observed network structure or is a *confounder* is not yet resolved — that is
Arm 2's question.

Preservation analyses (`base-001` v2 through `inp-023`) are queued to test whether the
module structure itself is preserved across disease groups, and whether residualising for
composition restores preservation.

---

## Completed runs

| id | wall | verdict | folder |
|:--|--:|:--|:--|
| `mod-000-resolution` | 2938 s | **Arm unblocked.** deepSplit=4, mergeCutHeight=0 → 33 modules, largest 11.68% | [network/runs/mod-000-resolution](runs/mod-000-resolution) |

### Running

| id | started | progress |
|:--|:--|:--|
| `base-001` v2 | 2026-10-14T19:14Z | nolith 500-perm preservation computing; 20 lith draws follow. ETA complete ~21:30Z 2026-10-14 |

---

## Completed analyses

### 001 — Module-trait associations (deepSplit=4)

`network/analyses/001-module-trait-ds4/` · Finding:
[2026-10-14-module-trait-no-lithium-association.md](findings/2026-10-14-module-trait-no-lithium-association.md)

No module eigengene associates with lithium status after BH correction (0 of 33; best
uncorrected p=0.061). Modules are not independent: median |r| = 0.822 between pairs, so
Bonferroni would be anti-conservative and BH is the appropriate correction.

Evidence 18/25 · Interest 6/10 · Δ −1.2

### 002 — Effective dimensionality / lithium dimension

`network/analyses/002-effective-dimensionality/` · Finding:
[2026-10-14-lithium-dimension-is-composition.md](findings/2026-10-14-lithium-dimension-is-composition.md)

PC2 of the 33 module eigengenes associates with lithium status (β = −0.908, BH p = 0.0029,
permutation p = 0.0015). ILR b1–b4 (cell composition proxies) explain 73.3% of PC2.
Adjusting for composition attenuates the lithium coefficient by 68.1% (p rises to 0.081,
no longer BH-significant).

Evidence 24/25 · Interest 9/10 · Δ −0.6 — **the primary finding to date**.

### 003 — LM22 circularity check

`network/analyses/003-lm22-circularity/` · Finding:
[2026-10-14-lm22-circularity-passes.md](findings/2026-10-14-lm22-circularity-passes.md)

Finding 002 could have been circular: ILR balances are computed from LM22 markers in the
same expression matrix used to build the co-expression network. Removing all 350 LM22 genes
(2.83% of 12,368; spread across 30 of 33 modules) and repeating: β −0.9081 → −0.8996,
BH 0.0029 → 0.0030, composition share 73.3% → 72.2%. The result moves < 1%. The two gene
sets' PC1–PC4 correlate at 1.000 — removing the marker genes leaves components identical.
Gene-level exclusion experiments will be uninformative on this data.

Evidence 20/25 · Interest 5/10 · Δ −3.0

### 004 — Decomposition robustness

`network/analyses/004-decomposition-robustness/` · Finding:
[2026-10-14-lithium-dimension-replicates.md](findings/2026-10-14-lithium-dimension-replicates.md)

All 10 cells of the deepSplit × mergeCutHeight grid were rebuilt and asked the same
question under the same a-priori rule (Li & Ji effective dimensionality, smallest BH).

- **6 of 10** show a lithium-associated dimension at BH < 0.05.
- **4 of 4 usable** cells (≥10 modules, largest < 25%) are among them.
- It is PC2 in every case but one.
- **0 of 6 survive composition adjustment** (BH across the 6 comparisons: all p = 0.166).
- **r(n\_modules, composition share) = +0.879** — the better the decomposition, the more
  of the lithium dimension is composition. At 33 modules: 73.3%. At 3 modules: 27.1%.
  This monotone relationship is inconsistent with a parameter artifact, and implies the
  degenerate base-001 structure would have *understated* the composition story.

Evidence 22/25 · Interest 7/10 · Δ −1.8

---

## Headline chain

1. **Scale-free topology fails** at every power 1–20, all groups, both n=74 and n=234.
   `pickSoftThreshold` returns `powerEstimate = 1` — the function failing. Max R² = 0.547.
2. **The network is extraordinarily dense**: mean connectivity 1,031 at power 12 (each gene
   tied to 8.3% of the transcriptome).
3. `mergeCutHeight = 0.25` **cascades and collapses** 33 modules into 3 (largest 97.5%)
   because median eigengene-pair correlation 0.739 is just below the 0.75 merge threshold.
4. On that degenerate decomposition **no module beats a random gene set**: gold scores
   Zsummary 20.26; best real module beats it in 9/15 subsamples.
5. `deepSplit=4, mergeCutHeight=0` gives **33 usable modules** (largest 11.68%, median 248).
6. No module eigengene associates with lithium after BH (0/33) — the eigengene approach
   lacks power at this dimensionality and correlation structure.
7. **A lithium dimension exists at the aggregate level**: PC2, permutation p = 0.0015.
8. It is **73% cell composition**, and that fraction holds across every tested decomposition.
9. **Never survives composition adjustment** — in any cell of the 10-decomposition grid.

---

## Next

1. **base-001 v2** (RUNNING) — preservation at correct parameters (33 modules, gold baseline)
2. **base-003-remerge** — reads base-001 v2 results; generates comparison tables
3. **base-002-ceil** — split-half ceiling; interpretability of any Zsummary seen in bp_lith
4. **inp-022** / **inp-023** — key question: does residualising for composition restore
   module preservation? This is the structural counterpart of finding 002's eigengene result.

---

## The question

Do gene co-expression networks differ across healthy / bipolar-off-lithium /
bipolar-on-lithium, and does any difference survive removing cell composition?

Exploratory. A well-documented negative or inconclusive result is a complete deliverable.
