# BP disease drives network-wide module differences; lithium adds nothing beyond composition

> **RETRACTED 2026-10-15.** 31/33 modules collapses to 0/33 once the per-sample library offset or the project's own technical covariate set is restored. The control-vs-bp_lith half also fails under the offset (2/33 Bonferroni, 0/33 BH), so no part of the finding is salvageable as written. The closing claim that the signature is 'present in patients who never took lithium' additionally violates IMAC_README §6.7: bp_nolith is currently-off-lithium, and the deposit records lithium status only.
>
> Verified independently by re-running this finding's own saved eigengenes against
> `data/analysis_matrices/covariates_474.csv` and `TMM_normalisation_factors.csv`.
> See `network/findings/2026-10-15-library-offset-artifact.md` and
> `network/analyses/015-referee-verification/`.
>
> Text below is left unaltered as the record of what was claimed.

> **CORRECTION 2026-10-15 — `b2` was mislabelled.** The sequential binary partition in
> `mediation/scripts/00b_ilr_trace.R:96-99` defines **b2 = gran vs mono**, not a B-cell or
> lymphocyte balance (b1 = myeloid vs lymphoid, b3 = T vs NK+B, b4 = NK vs B). Labels below
> have been corrected in place. Two directional claims built on the wrong label were also
> wrong and are flagged inline: b2 by group is bp_lith 0.3079 > control 0.2122 > bp_nolith
> 0.1816 (lithium is **highest**, not lowest), and r(darkturquoise, b2) is **negative** in
> every group (-0.130 ctrl, -0.390 nolith, -0.444 lith). The corrected reading — lithium
> raises granulocytes, granulocytes dilute the lymphocyte ribosomal program, darkturquoise
> falls — fits the arm's own r(darkturquoise, gran) = -0.510 (p = 2.3e-16).

**Date:** 2026-10-14
**Analysis:** `network/analyses/009-all-modules-across-groups`
**Evidence:** 0/25 · **Interest:** 0/10  _(retracted; was 19/25 and 9/10)_ · **Δ:** +1.4

## Claim

Extending the darkturquoise analysis (007–008) to all 33 modules reveals that the
group effect in module eigengenes is near-global, but its origin is different from
what the BP-only analysis suggested:

| Contrast | Significant modules (bonf < 0.05, Bonferroni n=3 within module) |
|:--|--:|
| control vs bp_nolith | **31 / 33** |
| control vs bp_lith | **32 / 33** |
| bp_nolith vs bp_lith | **0 / 33** |

Models: ME ~ group + b1+b2+b3+b4 for all 460 samples (additive composition adjustment).
Group ANOVA FDR (BH across 33 modules): 32 of 33 significant.

## The decomposition that matters

The three-way breakdown separates two effects that were conflated in the BP-only analysis:

1. **BP disease per se**: 31/33 modules differ between controls and untreated BP (ctrl-nolith),
   after controlling for cell composition (R²_comp 3–25% for these modules, so composition
   is not the dominant driver).

2. **No lithium-specific effect within BP**: 0/33 modules differ between nolith and lith
   within BP (after composition adjustment) — consistent with finding 002.

3. **Darkturquoise is uniquely lithium-associated**: it is the **only** module where
   ctrl-nolith is NS (bonf=0.814) but ctrl-lith is significant (bonf=0.001). This
   pattern is opposite to almost every other module.

## Darkturquoise's unique profile

| Module | ctrl-nolith bonf | ctrl-lith bonf | R²_comp |
|:--|--:|--:|--:|
| darkturquoise | 0.814 (NS) | **0.001** | **0.326** |
| All 31 other sig modules | < 0.015 | < 0.001 | 0.03–0.25 |
| darkred | 1.000 | 1.000 | 0.046 |

Darkturquoise is the highest-R²_comp module (32.6%). Its ctrl-lith significance is
likely composition-driven — analyses 007–008 showed the b2 slope is steeper in bp_lith
than controls. In contrast, the 31 ctrl-nolith significant modules have R²_comp 3–25%,
suggesting their group differences are NOT primarily explained by composition.

## Direction of the effect

All 31 ctrl-nolith significant modules show BP (nolith) **lower** than controls
(unidirectional, consistent across the full set):

| Module | ctrl ME | nolith ME | lith ME | Δ(nolith−ctrl) |
|:--|--:|--:|--:|--:|
| darkgrey | +0.0114 | −0.0116 | −0.0117 | −0.023 |
| darkolivegreen | +0.0101 | −0.0126 | −0.0092 | −0.023 |
| magenta | +0.0115 | −0.0109 | −0.0117 | −0.022 |
| midnightblue | +0.0142 | −0.0082 | −0.0167 | −0.022 |
| ... | | | | |
| yellow | +0.0138 | −0.0020 | −0.0176 | −0.016 |

Positive delta: 0 of 31. Negative delta: 31 of 31.

This is a consistent, global decrease in module eigengene values in BP patients relative
to controls — not a reshuffling of network structure but a monotone depression.

## Biological interpretation

**BP disease** (not lithium treatment) is associated with a broad restructuring of
peripheral blood co-expression — affecting 31/33 modules. These modules mostly have low
R²_composition (median ~10%), suggesting the difference is not fully explained by the
LM22 ILR composition proxy.

Within BP patients, **lithium does not further restructure the network** (0/33 nolith-lith
significant). The nolith-lith signal seen in analysis 002 (PC2, p=0.0015) is composition-
mediated.

Darkturquoise's unique ctrl-nolith NS / ctrl-lith significant profile is mechanistically
explained by its strong composition dependency: lith patients have lower b2 **[WRONG: lith b2 = 0.3079 is the HIGHEST of the three groups]** (granulocyte/
lymphocyte balance), and darkturquoise tracks b2 more steeply in lith patients, pushing
the module down and creating a ctrl-lith gap.

## Scoring

E=4 (pre-specified; N=460; additive model; shown to be misspecified for darkturquoise but
likely adequate for low-R²_comp modules; single cohort).
R=5 (reframes the core question: BP disease drives the broad network change; lithium
doesn't — the analysis 002 nolith-lith signal is composition, not the main story).
B=2 (single cohort, cross-sectional; no replication).
D=4 (31/33 ctrl-nolith, 0/33 nolith-lith, darkturquoise uniquely ctrl-lith specific).
C=4 (all 33 modules, all 3 contrasts, R²_comp reported).

Interest 9: The finding reframes everything. BP disease per se — not the drug — changes the
co-expression network. This is the headline. Δ +1.4 — ahead of the data until replicated.
