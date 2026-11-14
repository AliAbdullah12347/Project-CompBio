# Centrality separates diagnosis — and between a third and two thirds of it is age, tobacco and RNA quality

**Date:** 2026-11-14
**Runs:** `network/runs/netmet-060` (raw) · `netmet-062` (covariate-residualised)
**Aggregation:** `netmet-061` · `netmet-063`
**Scripts:** `network/experiments/netmet-060.R`, `netmet-061.R` (+ `netmet-062.R`, `netmet-063.R` shims)

## What was asked

Seven metric families over bipolar and non-bipolar networks, then the same for
lithium: module eigengene correlation, intramodular connectivity, gene
significance, topological overlap, differential connectivity, network
propagation proximity, and degree / betweenness / maximal-clique centrality.

## The design, and why it is not cent-040 again

`cent-040` reported **52 of 71 modules differing** in degree between control and
lithium. `cent-056` then split the control group in half -- same people, nothing
to find -- and got **43**. Two structural causes, neither statistical:

1. **Unmatched n.** control had 234 samples against bp_nolith's 74. bicor
   variance falls with n, so the 234-sample adjacency is systematically tighter
   and every derived metric inherits it.
2. **No same-kind reference.** "52 modules differ" has no scale without knowing
   what the procedure returns on two samples that cannot differ.

So every group here is cut to **n = 74**, and each of 10 draws builds **four**
networks -- `ctrlA`, `ctrlB` (disjoint controls), `nolith`, `lith` -- giving all
six pairwise comparisons, of which two are within-kind and four cross-diagnosis.

**The headline test is cross-diagnosis pairs against within-kind pairs, paired
by draw.** It is structurally fair in a way the ctrlA-vs-ctrlB null is not: half
the pairs on each side are random-vs-random and half are fixed-vs-random, so the
sampling-noise asymmetry cancels. (`nolith` is the only group not redrawn, being
all 74 of them.) The pre-specified NULL contrast is reported too, but it draws
**both** sides where the real contrasts draw one, so it carries more noise than
either and is conservative for BIPOLAR and lenient for LITHIUM.

## Result 1 — diagnosis separates on every metric, raw

Spearman agreement between the two groups of a pair, mean over 10 draws.

| metric | within-kind | cross-diagnosis | gap | draws | p |
|:--|--:|--:|--:|--:|--:|
| eigenvector | 0.835 | 0.502 | **0.333** | 9/10 | 0.008 |
| degree | 0.867 | 0.730 | 0.137 | 10/10 | 0.006 |
| PPR proximity (BPD seeds) | 0.844 | 0.742 | 0.102 | 10/10 | 0.006 |
| log-MCC | 0.783 | 0.702 | 0.080 | 10/10 | 0.006 |
| kIN | 0.907 | 0.834 | 0.073 | 10/10 | 0.006 |
| PPR proximity (LI seeds) | 0.891 | 0.823 | 0.068 | 10/10 | 0.006 |
| betweenness | 0.707 | 0.645 | 0.062 | 10/10 | 0.006 |
| kME | 0.892 | 0.839 | 0.053 | 10/10 | 0.006 |

All eight, BH < 0.007.

## Result 2 — and between a third and two thirds of it is covariates

The cross-diagnosis pairs are not covariate-matched and the within-kind pairs
nearly are:

| | age gap | tobacco gap | RIN gap |
|:--|--:|--:|--:|
| ctrlA vs ctrlB | 1.2 yr | 4.3 pts | 0.07 |
| **ctrlA vs nolith** | **9.5 yr** | **16.8 pts** | **0.20** |
| nolith vs lith | 5.1 yr | 3.9 pts | 0.05 |

Matching the controls is not available: only **54 of 234** reach the patients'
median age (54.8) and the design needs 148. So `netmet-062` residualises
expression on `age + sex + tobacco + rin + plate + seqpc1-3` over all 474
samples first -- fitted once and subset, never per group, for the same reason
the library offset is. Verified: max |cor(gene, age)| falls 0.445 to 0.000.

| metric | raw gap | residualised gap | removed |
|:--|--:|--:|--:|
| degree | 0.137 | 0.049 | 64% |
| kIN | 0.073 | 0.028 | 62% |
| PPR (LI) | 0.068 | 0.027 | 60% |
| log-MCC | 0.080 | 0.038 | 53% |
| kME | 0.053 | 0.025 | 53% |
| PPR (BPD) | 0.102 | 0.050 | 51% |
| betweenness | 0.062 | 0.032 | 48% |
| eigenvector | 0.333 | 0.222 | 33% |

**Every metric survives** -- all eight still separate at BH < 0.04, 10 of 10
draws except eigenvector at 8 of 10 -- but roughly half the apparent effect was
demographics and RNA quality. Anyone reporting the raw numbers as a bipolar
network signature would be reporting age and tobacco for about half of it.

Residualisation is conservative here, not neutral: age and tobacco both
correlate with diagnosis, so it absorbs part of the diagnosis effect too. The
surviving gap is a lower bound.

## Result 3 — lithium does not restructure the network

`nolith vs lith` is one of the two **within-kind** pairs, and it behaves like
one. It never falls below the within-kind baseline on any metric, raw or
residualised, and against the pre-specified NULL it sits at or **above** it on
all ten quantities. Lithium leaves centrality, connectivity, topological
overlap and propagation proximity where it found them.

That agrees with the rest of the arm -- 0 of 25 parameter variants, 0 of 71
modules, 0 of 12,360 genes, 0 of 5 lineages -- and with the finding that
lithium's transcriptomic effect is large, global and compositional. A mixture
shift moves every gene a little and reorders none of them.

## Result 4 (POST-HOC, exploratory) — lithium-treated patients sit nearer controls

Not among the three pre-specified contrasts. Noticed in the six-pair matrix and
reported as a hypothesis, not a result.

After residualisation, control-vs-lith agreement exceeds control-vs-nolith on
**7 of 8 metrics** (BH 0.019-0.037, 8-9 of 10 draws), most strongly on
eigenvector centrality (0.382 to 0.575). In the raw data the same comparison is
null or slightly negative.

Two reasons to be careful, and one reason it is not nothing:

- It is post-hoc, from one dataset, at modest significance.
- Residualisation removing the dominant demographic axis is exactly the kind of
  operation that can manufacture structure in what is left.
- But the sampling artifact runs **against** it: `nolith` is fixed and `lith` is
  redrawn, so control-vs-lith carries more noise and should score *lower*.

If it holds up it would be the first thing in this arm suggesting lithium moves
patients toward the control configuration rather than simply adding an effect of
its own. It needs an independent cohort, not another test here.

## Hub genes, aggregated over all 10 draws

Single-draw hubs are noise -- `grp-051` found two halves of one control group
share only 29% of their top-200. These are mean ranks across draws.

| metric | control | bipolar off-lithium | on-lithium |
|:--|:--|:--|:--|
| degree | UFM1, THAP12, PCNP, MICU2, SRP9 | **MYO1F, IL6R, OSBPL2, AGO4, NCF2, VSIR** | CBL, KMT2C, MED13L, IGF2R, KDM2A |
| eigenvector | UFM1, THAP12, MICU2, ZBTB6, PCNP | **OSBPL2, MYO1F, IL6R, RASSF2, NCF2, LYN** | ACAP2, TAOK1, KMT2C, MED13L, DOCK2 |
| betweenness | ZBTB25, RPL18, OXLD1, SIVA1, JMJD4 | GVQW3, NDUFS7, CFAP68, ZFTRAF1, S100A9 | RPL18, SPAG16, ABI2, NDUFS7, APH1A |

The off-lithium hub set is myeloid -- `MYO1F`, `IL6R`, `NCF2`, `VSIR`, `LYN`,
`PTPRE`, in modules `white`, `ivory`, `lightpink` -- where controls are
housekeeping and lithium is chromatin and signalling. Descriptive only: these
are each group's own network, and no null has been run on hub identity here.

## Three cautions on the metrics themselves

- **log-MCC is size-biased.** It returns the largest module (`turquoise`,
  `blue`) as the top hubs in every group, because a bigger module has more
  maximal cliques. Rank it within modules, never globally.
- **Eigenvector centrality carries the biggest effect and the most noise.** Two
  control halves agree at only 0.639 after residualisation, against 0.808 for
  degree. It is the leading axis of the graph, so it moves most and is least
  stable -- consistent with the arm's earlier finding that PC1 captures disease.
- **The "modules differing" counts are the cent-040 statistic** and are reported
  for continuity only. They inherit the draw-structure asymmetry and should not
  be quoted without the within-kind baseline beside them.

## What is NOT established

- That the surviving difference is diagnosis. It is everything that
  distinguishes a patient from a control and is not in the covariate set: the
  non-random QC filter on controls (35 of 240 dropped, Fisher p = 3.04e-10),
  other medications, illness state. There is no BMI and no other-medication
  field in this deposit.
- That the modules themselves differ. They do not -- `base-002-ceil` found no
  preservation difference at this n, and `me_cor_spearman` barely moves here
  (0.915 to 0.908). The picture is a network that keeps its modules and rewires
  inside them.
- Anything about effect size in biological units. These are rank correlations
  between centrality vectors.
