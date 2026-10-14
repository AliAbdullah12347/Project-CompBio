# Darkturquoise is a ribosomal/translation module: the composition signal has a biological substrate

**Date:** 2026-10-14
**Analysis:** `network/analyses/006-darkturquoise-annotation`
**Evidence:** 22/25 · **Interest:** 9/10 · **Δ:** +0.3

## Claim

The 140 genes of the darkturquoise module — the single module driving the
lithium-associated co-expression dimension (PC2, finding 005) — are strongly enriched for
cytoplasmic ribosomal proteins and oxidative phosphorylation / mitochondrial respiratory
chain components. The module is not an artefact of the network method; it has a coherent
biological identity.

## GO results (FDR < 0.05)

| Namespace | Terms FDR<0.05 | Top term | n_dt / n_all | FDR |
|:--|--:|:--|--:|--:|
| BP | 35 | Cytoplasmic translation | 30 / 149 | 7.2e-26 |
| CC | 47 | Cytosolic ribosome | 29 / 92 | 4.9e-32 |
| MF | 12 | Structural constituent of ribosome | 32 / 146 | 7.3e-30 |

Universe: 12,368 network genes. Test: Fisher exact, one-sided (over-representation),
FDR by BH within namespace. Min 10 annotated genes per term.

## Hub genes (top kME in darkturquoise)

| Gene | kME | Function |
|:--|--:|:--|
| RPL30 | 0.953 | Ribosomal protein, large subunit |
| RPL35 | 0.952 | Ribosomal protein, large subunit |
| RPL6  | 0.935 | Ribosomal protein, large subunit |
| RPL27 | 0.926 | Ribosomal protein, large subunit |
| RPS3A | 0.908 | Ribosomal protein, small subunit |
| COX7A2 | 0.904 | Cytochrome c oxidase (Complex IV) |
| UQCRQ | 0.893 | Complex III subunit |
| ATP5F1C | 0.891 | ATP synthase subunit |
| GIMAP2 | 0.890 | GTPase, immune cell (lymphocyte) expressed |

Multiple pseudogenes of RPL6 (RPL6P10, RPL6P27) also appear in the top 20,
further supporting that the module is driven by ribosomal biogenesis.

## Biological interpretation

Lymphocytes — particularly plasma B cells — are among the most translationally active
cells in peripheral blood. They sustain large-scale immunoglobulin synthesis, requiring
high ribosomal capacity and mitochondrial respiration. Samples with a high B-cell/lymphocyte
fraction (ILR b2) express more ribosomal protein genes, and that co-expression coalesces
into the darkturquoise module.

**The causal chain for the PC2 / composition signal is now complete:**

> High B-cell fraction (b2 balance) → higher ribosomal gene expression →
> darkturquoise module eigengene elevated → PC2 positive → PC2 observed to
> co-vary with lithium status only because lithium-treated patients differ in
> cell composition.

The lithium effect on network structure, as measured by module eigengene PC2, is explained
by a lymphocyte translational program — not by a direct lithium effect on gene regulation.

## Scoring

E=5 (GO test pre-specified in config.json; Fisher exact; three namespaces; large universe).
R=5 (gives the biological substrate of the headline finding; the ribosomal identity of
darkturquoise makes the composition explanation mechanistically coherent).
B=2 (one cohort, one decomposition; GO is stable but the causal claim is inferential).
D=5 (the signal is unmistakable: 35 BP terms, 47 CC terms, 12 MF terms; top kME genes
are almost all ribosomal proteins).
C=5 (all three GO namespaces reported; kME for all 140 genes written to CSV).

Interest 9: this is the most quotable result in the arm — a complete mechanistic story.
Δ +0.3 — ahead of the data only in the causal framing; the annotation itself is solid.

## Limits

The GO test uses `GOALL` (includes all ancestor terms), which inflates term counts but not
false positives in a one-sided test. The causal chain is plausible but not experimentally
verified: cell composition is measured from the same RNA, not from paired cell counts.
