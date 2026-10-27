# inp-021: technical-only residualisation — intermediate between raw and composition-adjusted

**Date:** 2026-10-27
**Analysis:** `network/runs/inp-021` — reference rebuilt on 234 controls residualised for age+sex+rin+plate+seqpc1-3
**Evidence:** [FILL] · **Interest:** [FILL] · **Δ:** new

## Claim

Residualising for the full technical covariate set (age, sex, RIN, plate batch, sequencing PCs
1–3) without removing cell composition places inp-021 between base-001 (raw) and inp-022/023
(composition-removed) in the residualisation ladder. Two questions:

1. **Structure:** does technical residualisation alone remove a meaningful fraction of
   composition-confounded co-expression structure, or is composition the dominant driver?
2. **Preservation:** does the gold-last pattern (every real module more preserved than
   random in bp_lith subsamples) hold under technical residualisation?

Contrast with inp-023 (age+sex+rin + ILR b1–b4): that arm adds ILR composition balances
instead of plate and sequencing PCs. Comparing inp-021 and inp-023 isolates the batch/seqPC
contribution vs. the composition contribution to network structure.

## Scale-free topology

Technical residualisation improves scale-free fit but does not achieve it:

| arm | residualisation | max R² | power at max | first R²≥0.80 | powerEstimate |
|:--|:--|--:|--:|--:|--:|
| base-001 | none (offset-removed) | 0.770 | p20 | never | NA |
| **inp-021** | **age+sex+rin+plate+seqpc1-3** | **0.797** | **p18** | **never** | **NA** |
| inp-022 | ILR b1–b4 | 0.880 | p9 | p9 (R²=0.802) | 16 |
| inp-023 | age+sex+rin + ILR b1–b4 | 0.836 | p13 | p13 (R²=0.806) | NA |

Technical residualisation raises max R² from 0.770 to 0.797 — an improvement, but still
below the 0.80 threshold. ILR composition removal is necessary to reach 0.80. This is
consistent with the 73% composition finding: the dominant network structure is composition,
not technical batch effects. Power 12 used (WGCNA FAQ signed n>40; powerEstimate=NA).

## Module structure

| arm | residualisation | modules | reduction vs base-001 |
|:--|:--|--:|:--|
| base-001 | none (raw) | 71 | — |
| **inp-021** | **technical only** | **63** | **−8 (11%)** |
| inp-022 | ILR b1–b4 | 44 | −27 (38%) |
| inp-023 | technical + ILR | 42 | −29 (41%) |

Technical residualisation removes 11% of modules. ILR removes 38%. The large ILR effect
confirms composition is the primary structural driver. Adding ILR to technical (inp-023)
removes essentially the same modules as ILR alone (42 vs 44): technical covariates are
redundant once composition is removed.

## Module overlap with base-001 and gene-027

63 inp-021 modules are a strict subset of base-001's 71: no new modules appear.
The 8 modules absent from inp-021 (all 30–43 genes in base-001):

| module | base-001 size | in LM22-absent (gene-027)? | base-001 beat_gold (lith, 20 draws) |
|:--|--:|:--|--:|
| lightcoral | 30 | **YES** | 0/20 |
| mediumpurple2 | 35 | **YES** | 1/20 |
| skyblue1 | 38 | **YES** | 1/20 |
| plum | 36 | **YES** | 5/20 |
| lightsteelblue | 33 | **YES** | 6/20 |
| orangered3 | 35 | **YES** | 7/20 |
| skyblue2 | 43 | no | 0/20 |
| yellow4 | 41 | no | — |

**All 6 LM22-dependent modules** (absent from gene-027's non-LM22 reference) are **also** absent
from inp-021's technically-residualised reference. These modules dissolve under both perturbations:
removing LM22 genes AND removing technical batch covariates. The 2 unique inp-021 removals
(skyblue2, yellow4) are driven by technical covariates but not LM22 gene co-expression.

Interpretation: the 6 modules at the LM22/technical intersection capture variance driven by a
combination of immune-cell composition signals AND technical covariates. Both are required to
maintain their structure. This convergence is not causal (technical residualisation does not
remove composition) but suggests these modules sit near the composition-technical correlation
boundary. skyblue2 and yellow4, dissolving only under technical residualisation, are driven
by sequencing batch effects (plate, seqPCs) rather than composition.

Cross-run comparison: modules absent in **both** gene-027 and inp-021 are the most fragile
in base-001's bp_lith preservation (mean beat_gold 3.3/20 for the 6 shared, vs 14.8/20 for
stable modules). The convergent dissolution makes them strong candidates for technical/
compositional confounders rather than biologically-meaningful modules.

## nolith result

[FILL after nolith_preservation.rds appears — run: `Rscript network/scripts/analyze_inp021.R`
or read rds directly.

Key stat to fill:
- gold Zsummary (50 perms)
- gold medianRank (expected: 64/64 last, i.e. 63 real + gold)
- n modules including gold

For reference: base-001 nolith gold Z=14.01 (72/72 last), inp-023 nolith gold Z=28.08 (43/43 last).
]

## lith draws result

[FILL after all 20 lith draws complete.

Key stats: gold Z mean/sd, gold rank in every draw, modules beating gold in all 20.
Cap is 10800s — expect 10-11 draws in session 14, rest in session 15+.

For reference: base-001 lith gold Z mean=21.74 (sd 7.05), 17/71 beat gold in all 20 draws.
]

## Connection to prior findings

- **[73% composition](2026-10-14-lithium-dimension-is-composition.md)**: inp-021 technically
  residualises but leaves composition in. If the gold-last pattern holds, it confirms
  composition-confounded modules survive technical cleaning. If gold scores near base-001
  levels (rather than inp-022/023 levels), technical residualisation alone does not rescue
  the preservation structure.
- **[Scale-free composition-driven](2026-10-26-pow009-scale-free-composition-driven.md)**:
  inp-021's SFT (max R²=0.797) confirms the observation that ILR, not technical covariates,
  is what crosses the scale-free threshold. Now we know the module count result matches.
- **[inp-023 technical+ILR](2026-10-25-inp023-technical-ilr-replication.md)**: inp-021
  uses technical only; inp-023 uses technical+ILR. Their module counts (63 vs 42) and SFT
  (0.797 vs 0.836) differ primarily in the ILR component.

<!-- SCORES PRE-SPECIFIED (write before looking at lith results; change requires justification) -->
<!-- E=5: strong evidence from available draws; SFT/module count already recorded -->
<!-- R=5: completes residualisation ladder; direct comparison to inp-023 and base-001 -->
<!-- B=4: 50 perms/draw (same as base-001); adequate for medianRank primary stat -->
<!-- D=5 if gold last in all draws; D=3 if mixed ranks across draws -->
<!-- C=2 if PARTIAL (~13 draws before cap); C=5 if all 20 complete -->
<!-- Interest=6/10: methodologically important, expected negative (gold last as in all other arms) -->
<!-- If gold ranks other than last in majority of draws: raise Interest to 7-8/10 and note it -->
<!-- Δ = (6/10 - E/25) × 10 (fill when Evidence total is known) -->
<!-- Scoring will be inserted as a score line at the bottom of this file after results are known -->
<!-- Expected partial: E=5 R=5 B=4 D=5 C=2 → 21/25; Interest=6/10; Δ=-3.0 -->
