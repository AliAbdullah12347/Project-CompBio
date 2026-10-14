# Darkred is an erythrocyte module: stable because RBC genes are absent from immune composition

**Date:** 2026-10-14
**Analysis:** `network/analyses/011-darkred-annotation`
**Evidence:** 20/25 · **Interest:** 8/10 · **Δ:** 0.0

## Claim

Darkred (164 genes) is the only module with no group effect whatsoever (ctrl-nolith bonf=1.0,
ctrl-lith bonf=1.0, group ANOVA FDR=0.981; analysis 009). This is not random — darkred is
annotated as an **erythrocyte / red blood cell module**, and erythrocyte genes are immune to
the immune cell composition shifts that drive 31/33 other modules.

## GO enrichment (darkred, n=164, universe=12,368)

| Namespace | FDR<0.05 terms | Top term | n_dr / n_all | FDR |
|:--|--:|:--|--:|--:|
| BP | 52 | Erythrocyte development | 9 / 34 | 2.1e-6 |
| CC | 6 | Cortical cytoskeleton | 11 / 64 | 4.7e-7 |
| MF | 1 | Ankyrin binding | 4 / 11 | 0.009 |

Top BP terms: erythrocyte development (FDR=2.1e-6), erythrocyte differentiation,
hemoglobin metabolic process, porphyrin/heme metabolism, myeloid cell differentiation.

## Hub genes (kME in control samples, n=234)

| Rank | Gene | kME | Function |
|--:|:--|--:|:--|
| 1 | ADIPOR1 | 0.960 | Adiponectin receptor |
| 2 | DCAF12 | 0.952 | E3 ubiquitin ligase adaptor |
| 3 | TSPAN5 | 0.951 | Tetraspanin, erythroid |
| 4 | YBX3 | 0.942 | Y-box binding protein, erythroid |
| 5 | FBXO7 | 0.939 | F-box protein, erythropoiesis |
| 6 | **ANK1** | 0.937 | **Ankyrin-1: red blood cell membrane scaffold** |
| 7 | TRIM58 | 0.936 | E3 ligase, erythropoiesis |
| 8 | MPP1 | 0.934 | RBC junctional complex |
| 9 | GUCD1 | 0.934 | Guanylyl cyclase domain |
| 10 | GMPR | 0.924 | Guanosine monophosphate reductase, RBC |

ANK1 (ankyrin-1) is the archetypal red blood cell structural protein. Mutations in ANK1
cause hereditary spherocytosis — its high kME confirms the erythrocyte identity.

## Why darkred is stable across groups

Darkred ME group means: ctrl = +0.00024, nolith = −0.00074, lith = −0.00096 (trivial
differences, SD ≈ 0.045).

Darkred ME vs composition: r_b1=+0.057, r_b2=−0.097, r_b3=−0.024, r_b4=+0.122.

Contrast with darkturquoise: r_b1=−0.304, r_b2=−0.316 (5–6× stronger composition coupling).

**The causal explanation**: peripheral blood RNA-seq includes reticulocytes — young red blood
cells that retain ribosomes and RNA. Reticulocytes express erythrocyte structural and heme
synthesis genes. The LM22 cell type signature matrix (used to compute ILR balances) contains
no erythrocyte cell type — LM22 deconvolves 22 *immune* cell types. Therefore, erythrocyte
gene expression is invisible to the composition proxy, and darkred has near-zero composition
correlations. Because its expression is stable across immune cell mixtures, it doesn't change
between BP and control groups either.

## Completing the network picture

The network is primarily measuring immune cell biology:

| Module | Biological identity | Group effect | Composition |
|:--|:--|:--|:--|
| Darkturquoise (140) | Ribosomal / B-cell translation | ctrl-lith only (composition) | R²=33% |
| Darkred (164) | Erythrocyte / RBC membrane | **None** (FDR=0.981) | R²=4.6% |
| 31 others | Mixed immune | ctrl-nolith + ctrl-lith | R²=3–25% |

The 31 modules that change in BP reflect genuine immune cell biology; the 1 that doesn't
(darkred) is the one that's not measuring immune cells at all.

## Scoring

E=4 (GO enrichment well-powered; hub gene confirms ANK1; clear biological logic; single
cohort).
R=5 (closes the analysis 009 story; mechanistically explains the exception; validates the
immune-biology framing of the full 33-module network).
B=2 (single cohort; no external replication).
D=4 (top GO term FDR=2.1e-6; ANK1 as hub gene; composition correlations 5–6× lower than
darkturquoise; group effect strictly null FDR=0.981).
C=5 (GO across 3 namespaces; kME; composition comparison; group means tabulated; contrast
with darkturquoise).

Interest 8: Closes the exception from finding 009 with a crisp biological mechanism.
Δ +0.4 — evidence and interest near-matched; the erythrocyte story is strong.
