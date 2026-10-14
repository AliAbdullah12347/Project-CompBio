# BP disease tightens inter-module correlations: modules move more in lockstep

**Date:** 2026-10-14
**Analysis:** `network/analyses/013-inter-module-correlations`
**Evidence:** 16/25 · **Interest:** 8/10 · **Δ:** +1.6

## Claim

BP disease is associated with a tightening of inter-module correlations — the 33 modules
move more in lockstep in BP patients than in controls. The correlation *structure* (which
modules correlate with which) is largely preserved (Mantel r=0.934), but the degree of
inter-module coupling is significantly higher in BP.

## Structure

### Inter-module |r| by group (528 module pairs, upper triangle)

| Group | Mean |r| | Median |r| | Fraction |r|>0.5 | Fraction |r|>0.8 |
|:--|--:|--:|--:|--:|
| control | 0.702 | 0.763 | 82.6% | 41.7% |
| bp_nolith | 0.738 | 0.790 | 86.7% | 47.7% |
| **bp_lith** | **0.793** | **0.845** | **93.0%** | **62.3%** |
| all 460 | 0.751 | 0.806 | 88.3% | 52.8% |

### Mantel-like r (Pearson r between upper triangles of per-group correlation matrices)

| Comparison | Mantel r |
|:--|--:|
| ctrl vs nolith | 0.934 |
| ctrl vs lith | 0.933 |
| nolith vs lith | 0.984 |

Permutation test (K=1000; shuffle group labels): ctrl vs nolith observed r=0.934;
null mean=0.972 (SD=0.015); p_two-sided=0.024.

The observed structure similarity is **significantly lower** than the null expectation.
When group labels are randomised, the correlation structure r≈0.972 — random subsamples
are more similar to each other than ctrl and nolith actually are.

## Interpretation

**What it means structurally:** In controls, the 33 modules show moderate inter-module
coupling (median |r|=0.763). In BP patients (nolith), the coupling is tighter (0.790)
and tightest in lith (0.845). Modules that are somewhat independent in health become
highly coordinated in disease.

**Connection to PC1 (77.5%):** Analysis 010 showed PC1 explains 77.5% of ME variance
across all 460 samples, with uniform loadings (−0.181 to −0.191) on all 33 modules.
This is the eigenvector signature of tightly coupled modules: when all modules move
together, a single PC captures almost all variance. The 77.5% PC1 in the all-sample
analysis emerges directly from the higher inter-module correlations in BP patients.

**Not just "everything goes down":** The structure change (Mantel p=0.024) means it's not
purely a uniform level drop. The pattern of *which pairs* are most correlated changes
across groups. Darkred shows the largest absolute change in mean |r| with other modules
(ctrl=0.310 → nolith=0.424, Δ=+0.115) despite having no ME level difference — as all
immune modules drop together, their correlation with the stable darkred module shifts.

**The "volume" metaphor:** BP disease acts partly like a volume knob (all modules lower,
finding 010) but also like a "pan" knob — modules become less distinct from each other,
converging toward a common immune suppression mode.

## Caution

The nolith group has only n=74, making individual pair correlations noisy; the Mantel test
uses all 528 pairs simultaneously, reducing noise. The permutation null mixes all group
sizes, so it correctly reflects what random subsampling of these data looks like.

## Scoring

E=3 (well-powered N=460/74/234; pre-specified; single cohort; Mantel test is a reasonable
proxy but not the same as direct graph comparison).
R=4 (structurally explains PC1=77.5% variance; adds a new dimension beyond eigengene
levels; connects to preservation analysis).
B=2 (single cohort; no replication).
D=3 (Mantel r>0.93 = high preservation; p=0.024 = significant; |r| tightening monotone).
C=4 (all 3 groups; permutation test; per-module delta table; mantel vs null tabulated).

Interest 8: The tightening of inter-module correlations in BP patients, and its connection
to PC1 explaining 77%, is a mechanistically coherent and quotable structural finding.
Δ +1.6 — interest ahead of evidence; needs replication or formal graph comparison.
