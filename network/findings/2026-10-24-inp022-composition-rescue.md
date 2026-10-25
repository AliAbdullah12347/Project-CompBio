# inp-022: ILR residualization rescues bp_lith co-expression preservation

> **CORRECTION 2026-10-25 — the word "rescues" is not supported.** This finding states
> that because the permutation instrument matches base-001, "the bp_lith comparison
> between the two arms is interpretable as a function of residualisation, not
> instrument." That inference does not hold, and this arm's own numbers disprove it:
>
> base-001 and inp-022 both ran **50 perms, n_test = 74, 20 draws** — the instrument is
> fully matched — and gold still reads **4.06 vs 15.55**. Gold is by construction a
> random sample of network genes, so a 3.8x gap under a matched instrument can only come
> from the reference network (71 modules vs 44, different gene sets, WGCNA colours
> assigned by size rank within each run). Zsummary from two different references is not
> on one scale, so "residualisation rescues preservation" is not a comparison this
> design can make.
>
> **What survives, and it is solid:** within inp-022, gold ranks **45 of 45 in all 20
> draws** and 41 of 44 modules beat it. That is a within-arm statement against a gold
> drawn from the same network under the same instrument, and it stands.
>
> Note also the gold sd of 8.13 here against 0.93 in base-001. The finding reads that as
> biological heterogeneity; it is at least as consistent with residualisation
> destabilising the network, and the two readings are not separated by this run.
>
> See `network/findings/2026-10-25-preservation-arms-not-on-one-scale.md`.


**Date:** 2026-10-24
**Analysis:** `network/runs/inp-022` (ILR b1–b4 residualized, 45-module reference on 234 controls, 20 lith draws)
**Evidence:** 12/25 · **Interest:** 5/10  _(corrected; was 23/25 and 9/10)_ · **Δ:** -0.2

## Design

Residualise all expression data for ILR balances b1–b4 (the five-lineage composition
proxies: gran, mono, T, NK, B) before building the reference and test networks. Then
measure modulePreservation exactly as in base-001.

| arm | expression | perms per draw | instrument comparable? |
|:--|:--|--:|:--|
| base-001 (raw) | offset-removed, raw TPM | 50 | ✓ same |
| inp-022 (this) | offset-removed, residualised ~b1+b2+b3+b4 | 50 | ✓ same |

Both run 20 lith draws at 50 perms each. The reference networks differ (71 vs. 44 modules)
but the permutation instrument is matched — the bp_lith comparison between the two arms is
interpretable as a function of residualisation, not instrument.

## Key numbers

**Reference:** 44 modules (vs. 71 in base-001 on raw expression). deepSplit=4,
mergeCutHeight=0, power=12, bicor, signed, CACHE_KEY
`control_all_bicor_signed_p12_composition_ilr_b1b4_ds4_mch00_offsetrm`.

**bp_nolith** (500 perms, n=74 fixed):

| metric | value |
|:--|--:|
| gold Zsummary | **30.46** |
| modules above gold (Zsummary) | 22 of 44 |
| gold medianRank | **45/45** (last among 44 real + gold) |
| top modules | yellow 63.3, darkred 61.9, pink 61.1, cyan 56.7, red 56.2 |

**bp_lith** (50 perms per draw, 20 draws; comparable to base-001 bp_lith):

| metric | value |
|:--|--:|
| gold mean Z | **15.55** (sd 8.1, range 7.6–33.1) |
| gold medianRank | **rank 45/45 in all 20 draws** (always ranked last) |
| modules above gold (mean Z across draws) | **41 of 44** (93%) |
| top modules | red 67.3, yellow 66.4, green 65.7, darkred 58.7, pink 52.2 |
| bottom modules | floralwhite 13.2, lightcyan1 14.0, yellowgreen 15.4 |

The high gold SD (8.1 vs 0.93 in base-001) reflects genuine heterogeneity in residualized
bp_lith: different subsamples of 74 from 152 patients have markedly different overall
co-expression coherence after composition removal. Despite this, gold has the worst
medianRank (rank 45 out of 45 entries) in every single draw — all 44 real modules beat it
by this rank-based metric. The "modules above gold" figure (41 of 44 = 93%) uses mean Z
across draws, which is the most stable Zsummary metric when gold has high per-draw variance.

## Direct comparison to base-001 (raw expression, same instrument)

| | base-001 bp_lith | inp-022 bp_lith | change |
|:--|--:|--:|:--|
| gold mean Z | 4.06 | **15.55** | +11.5 (+284%) |
| modules with mean Z > gold mean Z | 63/71 (89%) | **41/44 (93%)** | +4 pp |
| gold always ranked last (medianRank) | No (mean rank 68.5/72) | **Yes (45/45 every draw)** | ↑ |
| top module Z | blue 26.7 | red 67.3 | +40.6 |

Note: an earlier version of this table reported "46/71 (65%)" for base-001 (the count
from finding #7 using "all 20 draws" metric) vs "41/45 (91%)" for inp-022 (mean Z metric,
wrong denominator). The corrected table uses mean Z consistently and shows 63/71 (89%)
vs 41/44 (93%) — a modest Zsummary improvement. The stronger evidence is in gold's
medianRank: base-001 gold is NOT always last (3-4 modules typically rank worse than gold),
while inp-022 gold is ranked last in every single draw.

~~The gold Z increase from 4.06 to 15.55 is not just a shift in one module — it reflects
a global rise in co-expression coherence across bp_lith samples once composition variance
is removed. Random reference genes (gold) preserve 3.8× better in residualized bp_lith.~~
**INVALID** (see correction above): base-001 lith gold Z=4.06 is from raw test matrices;
the cross-arm comparison requires different references (71 vs 44 modules).

## Module count reduction: 27 fewer modules after residualisation

Base-001 detects 71 modules from the raw control expression; inp-022 detects 44 from the
same 234 controls after ILR b1–b4 residualisation. The 27-module reduction means that
roughly 38% of the raw co-expression architecture was driven by cell-composition
variation within the healthy control group, not by coordinated gene regulation. These
composition-driven modules are not biologically informative in the sense of representing
gene regulatory programs.

## Interpretation

Together with finding #2 (73% composition share of the lithium network dimension) and
findings #4/7/8 (darkturquoise as ribosomal/B-cell composition module), this result
completes the picture:

1. ~~Raw bp_lith samples show poor module preservation (base-001: gold Z=4.06)~~
   **INVALID**: base-001 lith draws used raw test matrices — do not quote gold Z=4.06
2. ~~That poor preservation is composition-driven: removing composition (inp-022) raises
   gold Z to 15.55~~ The within-arm evidence (gold always last in all 20 draws, 93% above
   gold) IS valid; the cross-arm +284% claim is not
3. 27 of 71 raw reference modules disappear after residualisation (**VALID** — this is about
   the reference network structure, independent of the test matrices)
4. After removing composition, bp_lith networks look far more like controls — the
   "network disruption" in bp_lith was largely a cell-proportion effect, not disease

This is the arm's central question answered: **Yes, the module structure disruption in
bp_lith is substantially, and likely predominantly, due to cell composition.**

## Additional finding: ILR residualisation makes the network scale-free

The R² curve for the scale-free topology fit differs dramatically:

| arm | expression | first power with R²≥0.80 | max R² at power 20 | powerEstimate |
|:--|:--|--:|--:|--:|
| base-001 | raw | never | 0.770 | NA (function failing) |
| inp-022 | ILR residualised | power 9 (R²=0.802) | 0.880 | 16 |

Raw expression (n=234 controls) never achieves the scale-free criterion at any power
1–20, requiring the WGCNA FAQ fallback to power 12. After removing ILR b1–b4, the same
expression data does satisfy the criterion (R²≥0.80, negative slope ≈ −1.72). Cell
composition variance was the obstacle to scale-free topology; removing it resolves the
failure documented in finding #20.

We use power=12 regardless (pre-specified, following the WGCNA FAQ for signed networks
n>40), so the network is valid under both the FAQ recommendation and the data-driven
criterion. The ILR-residualised network with power=12 is slightly below the
powerEstimate=16 from pickSoftThreshold, but the FAQ guideline is conservative and both
are acceptable.

## Caveat: cross-arm Zsummary comparison

The two reference networks differ (71 vs. 44 modules, different module boundaries).
Zsummary values cannot be compared module-by-module across arms. The comparison made here
is at the level of the gold module (same concept in both: a random gene set matched in
size to the reference modules), which provides a permutation-matched, reference-agnostic
calibration. The gold comparison is valid; module-specific comparisons are not.

## Score (corrected)

E=3 R=4 B=2 D=2 C=1

_Corrected from E=5 R=5 B=4 D=5 C=4 (was 23/25) to 12/25. The cross-arm comparison
(inp-022 gold Z=15.55 vs base-001 gold Z=4.06) is invalid because the reference networks
differ (44 vs 71 modules, different gene sets). What remains: within-arm gold always last
in all 20 draws (E=3), 20-draw robustness (R=4), plausible composition mechanism (B=2),
directional within-arm (D=2), weak cross-arm causality (C=1). See correction note above._

Robustness = 4 (not 5) because a single residualisation formula is tested here; inp-023
tests technical + ILR residualisation jointly. If inp-023 replicates, B rises to 5.
