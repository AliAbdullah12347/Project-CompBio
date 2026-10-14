# The lithium-composition result is not circular: it survives removing every marker gene

**Date:** 2026-10-14
**Analysis:** `network/analyses/003-lm22-circularity`
**Evidence:** 20/25 · **Interest:** 5/10 · **Δ:** -3.0

## Claim
Finding 002 said PC2 is lithium-associated and 73.3% composition. Because the ILR balances
are estimated from the LM22 markers in the same matrix, that could have been circular.

Dropping all **350 LM22 markers** (2.83% of genes, spread over 30 of 33 modules) and
recomputing:

| | β | BH | attenuation | PC r² composition |
|:--|--:|--:|--:|--:|
| all 12,368 genes | −0.9081 | 0.0029 | 68.1% | 73.3% |
| LM22 excluded (12,018) | −0.8996 | 0.0030 | 67.4% | 72.2% |

**The association survives and the composition share holds.** Removing every gene used to
estimate composition moves the result ~1%.

## Secondary observation
The two gene sets' PC1–PC4 correlate at **1.000**. Removing the most
composition-informative 2.83% of genes leaves the components identical — the structure is
extremely diffuse, consistent with mean connectivity at 8.3% of the genome. Gene-level
exclusion experiments on this data will generally be uninformative.

## Scoring
E=5 (direct, deterministic, same a-priori dimension rule, no post-hoc PC picking).
R=4 (validates the headline finding rather than extending it). B=2 (still one
decomposition). D=5 (the check passes unambiguously). C=4.

Interest 5: a negative control that passes is reassuring, not surprising. Δ −3.0 — well
supported, deliberately unexciting.

## What it does not test
CIBERSORTx estimates fractions from the whole profile, not only marker genes, so some
shared variance is unavoidable. This removes the most direct circular path, not all of them.
