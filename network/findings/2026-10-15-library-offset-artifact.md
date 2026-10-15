# The reference network is built on an un-removed per-sample library offset

**Date:** 2026-10-15
**Analysis:** `network/analyses/015-referee-verification`
**Evidence:** 24/25 · **Interest:** 10/10 · **Δ:** new

## Claim

`expr_for()` in `network/00_load.R` returns raw `log2(TPM+1)`. The TMM normalisation
factors are loaded into `d$tmm` (line 76) and used in **exactly one** file across every
analysis and experiment (`inp-024.R`, still PENDING). The per-sample mean of
`log2(TPM+1)` is therefore left in the matrix, and it dominates the correlation structure
the entire network is built from.

Measured on the control reference set (1,500 random genes, n=234):

| | mean gene–gene r | % of random pairs with r > 0.5 |
|:--|--:|--:|
| raw `log2(TPM+1)` | **+0.535** | **61.9%** |
| per-sample mean removed | +0.002 | 5.4% |

Sixty-two percent of *randomly chosen* gene pairs correlate above 0.5 because of one
scalar per sample. Remove it and that falls to the null expectation.

That scalar is not neutral with respect to the research question:

| | control | bp_nolith | bp_lith |
|:--|--:|--:|--:|
| mean `log2(TPM+1)` | 3.4637 | 3.2688 | 3.1811 |

Welch one-way p = 1.14e-09, and it correlates r = +0.806 with the TMM normalisation
factor it was never divided by.

## What this invalidates

**Finding #12 (`pc1-captures-bp-disease`) is dead.** PC1 is 77.5% of eigengene variance
with all 33 loadings the same sign — a global scaling axis, not a biological factor.
r(PC1, TMM factor) = −0.828 (R² = 0.685) against r(PC1, disease) = +0.289. Partial the
TMM factor out and the disease correlation collapses from +0.289 to **+0.056**.

**Finding #8 (`bp-disease-drives-network-globally`) is dead.** Re-running analysis 009's
own model on its own saved eigengenes, control vs bp_nolith:

| model | raw p<.05 | Bonferroni×3 | BH across 33 |
|:--|--:|--:|--:|
| 009 as published (composition only) | 31/33 | 31/33 | 31/33 |
| + technical (age sex rin plate seqpc1–3) | 9/33 | 1/33 | **0/33** |
| + library offset | 0/33 | 0/33 | **0/33** |
| + both | 0/33 | 0/33 | **0/33** |

The control-vs-bp_lith half, which the referee panel's editor recommended keeping,
survives technical covariates (32/33) but **also collapses under the library offset**
(2/33 Bonferroni, 0/33 BH). It should not be kept either.

**It mechanically explains four existing findings** that were interpreted as results:
scale-free topology failing at every power (#16 — a uniformly dense graph cannot be
scale-free), no module beating a random gene set (#5), the merge threshold being decisive
(#3), and 33 modules spanning only ~6 effective dimensions with median inter-eigengene
|r| = 0.822 (analysis 002). These are all symptoms of one cause.

## What survives

**The headline finding (002) survives, and this is the important part.** Within BP1
(n=226), PC2 against lithium:

| model | β | p |
|:--|--:|--:|
| marginal | −0.9082 | 0.000473 |
| + composition | −0.2891 | 0.0394 |
| + library offset | −0.6450 | 0.00748 |
| + offset AND composition | −0.2904 | 0.0394 |

PC2 is only 10.6% explained by the TMM factor where PC1 is 68.5%. The lithium dimension
is orthogonal to the artifact axis, which is why it survives. The marginal β = −0.9082,
p = 4.73e-4 reproduces finding 002's reported β = −0.908, p = 4.8e-4 to three digits.

## Either branch damages the retracted findings

If the offset is a pure normalisation artifact, #8 and #12 are measuring library
composition. If it is instead real biological composition that LM22 cannot see — note it
is only weakly related to the ILR balances, and LM22 has no erythrocyte class — then #8's
conclusion *inverts*, because #8's entire argument is that the BP effect is **not**
compositional. There is no reading under which they stand as written.

## Status

Verified independently against the referee panel by re-running every number above from
the repository's own saved CSVs and cached project object. Agreement with the panel is
close to exact (panel: +0.544 / 62.9% / +0.001 / 5.0% / offsets 3.464, 3.269, 3.181).

**This is a rebuild of the reference network, not a revision of the writeup.** The
dendrogram, the 33 modules and every preservation statistic in the queue inherit it.

E=5 R=5 B=5 D=5 C=4
