# Alpha and beta diversity, mediated through cell mix — research and plan

**Status: research only. No analysis has been run.**
Two read-only checks of the data were made so that the plan below is correct for *this*
dataset (section 7). Nothing was written, fitted or saved by them.

Written 2026-10-18 after the weekly check-in.

---

## 0. The request, verbatim

> "Alpha diversity and beta diversity, and whether changes of alpha and beta diversity are
> mediated through cell mix. Do it in gene expression in mediation.
> Research more into that, understand how this kind of Alpha diversity is done on the gene
> expression level, and also how it is done on higher layer such as pathway analysis."

---

## 1. The task in one paragraph

Borrow two ideas from ecology. **Alpha diversity** is how varied *one* sample is;
**beta diversity** is how different samples are *from each other*. Apply both to blood gene
expression: alpha = how evenly a person's blood spreads its expression across genes, beta =
how different two people's whole expression profiles are. Then ask three things:
(1) does lithium change alpha diversity, (2) does lithium change beta diversity, and
(3) **is either change mediated through the cell mix** — the same exposure → mediator →
outcome structure as the existing mediation arm, but with transcriptome-wide *summaries* as
the outcome instead of 12,018 separate genes. Repeat at the pathway level, where "genes"
are replaced by biological programmes.

## 2. Why this is being asked — the aim

It answers the exact wall the gene-by-gene mediation hit.

| Problem found in the mediation arm | What diversity changes |
|:--|:--|
| 12,018 separate tests; per-gene significance needed \|z\| ≥ 4.60, and the best attainable was 2.92 | **One** test per diversity measure. Significance needs \|z\| ≥ 1.96, which the lithium→cell-mix link (t = 2.0–2.55) can in principle support |
| The count of "mediated genes" was meaningless — the null alone produced up to 7,150 | No gene counts. A single global answer |
| The proportion mediated was indistinguishable from noise (permutation p = 0.49) | The primary evidence becomes a direct test of the indirect effect, not a ratio |
| The bipolar-vs-control contrast found almost nothing gene by gene | A global measure can detect a diffuse shift that no single gene carries |

So the aim is: **ask the mediation question once, about the whole transcriptome, with enough
power to get an answer.**

---

## 3. The concepts

### 3.1 Where the words come from

Whittaker (1960, 1972) split biodiversity into three levels for plant communities:

| Level | Ecology meaning |
|:--|:--|
| **alpha** | diversity *within* one site — how many species, and how evenly abundant |
| **beta** | difference *between* sites — how different their species make-up is |
| **gamma** | diversity of the whole region pooled together |

Microbiome research adopted the same vocabulary, which is almost certainly where it came
from in the check-in: per-sample alpha diversity (Shannon, Simpson, richness) and
between-sample beta diversity (Bray–Curtis, UniFrac, Aitchison), tested with PERMANOVA.

### 3.2 The translation to gene expression

| Ecology | This project |
|:--|:--|
| a species | a **gene** |
| a species' abundance | the gene's **expression** (TPM or reads) |
| one site / community | **one person's blood sample** |
| alpha diversity of a site | how evenly that person's expression is spread across genes |
| beta diversity between sites | how different two people's expression profiles are |
| gamma | the pooled cohort |
| a higher taxonomic rank (genus, family) | a **pathway / gene set** |

### 3.3 Alpha diversity on gene expression — what already exists

**The founding definition.** Martínez & Reyes-Valdés (2008, *PNAS*) defined **transcriptome
diversity** as the Shannon entropy of a sample's expression frequency distribution:

```
H_j = − Σ_i p_ij · log(p_ij)        p_ij = expression of gene i / total expression in sample j
```

High H: expression spread evenly over many genes. Low H: dominated by a few highly expressed
genes. The same paper also defines **gene specificity** (mutual information between sample
and transcript) and **specialization** (expression-weighted mean specificity); a 2010 follow-up
found tumours *less* specialized than normal tissue.

**Why it matters for us specifically.** García-Nieto, Wang & Fraser (2022, *PLoS Comput Biol*)
showed transcriptome diversity is **one of the largest single sources of variation in RNA-seq
data** — it is the strongest factor captured by PEER hidden-confounder factors, and it
associates with both technical variables and biological ones such as cell composition. Two
consequences: lithium could plausibly move it (through composition), and technical
covariates **must** be adjusted for or the outcome will mostly measure the library.

**The modern unifying framework: Hill numbers** (Hill 1973; Jost 2006; Chao et al. 2014).
One formula, an order *q* that controls how much weight dominant genes get:

| q | Hill number | In expression terms |
|:--|:--|:--|
| 0 | richness | number of genes detected — very sensitive to sequencing depth |
| 1 | exp(Shannon H) | **"effective number of equally expressed genes"** |
| 2 | 1 / Σ p² (inverse Simpson) | effective number of *dominant* genes |

Hill numbers are in "effective number of genes", which is far easier to explain than nats.

### 3.4 Beta diversity on gene expression — what already exists

**Step 1, a distance between every pair of people.** Choices:

| Distance | Notes |
|:--|:--|
| **Aitchison** (Euclidean on log-ratio / CLR expression) | RNA-seq is compositional (Quinn et al. 2018). **Same geometry as our ILR mediator** — Euclidean distance between two people's ILR balances *is* the Aitchison distance between their cell mixes. Mediator and outcome measured on one coherent scale |
| Jensen–Shannon | the information-theory partner of Shannon alpha diversity |
| Bray–Curtis | the ecology default; dominated by abundant genes |
| 1 − correlation | common in expression clustering; ignores scale |

**Step 2, test it.** Two separate questions, two separate tests:

| Question | Test | Reference |
|:--|:--|:--|
| Do lithium users' profiles sit in a **different place**? (location) | **PERMANOVA** / distance-based redundancy analysis (dbRDA) | Anderson 2001; Legendre & Anderson 1999; McArdle & Anderson 2001 |
| Are lithium users **more or less varied among themselves**? (dispersion) | **PERMDISP** | Anderson 2006 |

This is not borrowed loosely: distance-based regression was introduced *for gene expression*
by Zapala & Schork (2006, *PNAS*, "MDMR"). Gene-level variance differences have their own
literature — "differential variability" (Ho et al. 2008; Phipson & Oshlack 2014) — and the
ecological idea that stressed systems vary *more* than healthy ones is the "Anna Karenina
principle" (Zaneveld et al. 2017).

**A trap that applies to us.** PERMANOVA assumes equal dispersion. Anderson & Walsh (2013)
showed it is robust to unequal dispersion in **balanced** designs but **not unbalanced** ones.
Ours is 152 vs 74. So PERMDISP is not optional: a significant PERMANOVA cannot be read as a
location shift until dispersion has been checked.

### 3.5 The higher layer: pathways

Genes → pathways is the expression analogue of species → genus or species → functional group.
There are two established routes, and they suit different questions.

| Route | How | Good for |
|:--|:--|:--|
| **Aggregate** | sum each gene's TPM into its pathway (a gene in *k* sets contributes 1/*k* to each, so nothing is double counted), then proportions across pathways | **alpha**: Shannon over pathways = how evenly blood spreads its expression across biological programmes |
| **Score** | single-sample pathway scores: ssGSEA (Barbie et al. 2009) or GSVA (Hänzelmann et al. 2013) | **beta**: a samples × pathways matrix, then Euclidean distance |

Two cautions found in the research:
- **No published precedent applies Shannon entropy to GSVA/ssGSEA scores**, and those scores
  can be negative, so entropy is undefined on them. Alpha must use the aggregation route.
- **Circularity is worse at this level.** Many curated sets (e.g. Hallmark inflammatory
  response, interferon, IL6/JAK/STAT3) are dominated by myeloid-cell genes, so a pathway score
  partly re-measures the cell mix. LM22 signature genes must be removed before aggregation as
  a sensitivity check.

Recommended gene sets: **MSigDB Hallmark** (50 curated, low-redundancy sets; Liberzon et al.
2015) as primary; Reactome top-level as sensitivity.

### 3.6 A decomposition that fits this problem exactly

This is the single most useful idea the research turned up.

Bulk blood is a **mixture** of cell types weighted by their fractions. For any mixture of
distributions, Shannon entropy splits exactly (this is the definition of the generalised
Jensen–Shannon divergence, Lin 1991):

```
H(bulk)  =  Σ_k f_k · H(cell type k)   +   JSD_f(cell types)
 gamma   =  alpha: within-cell-type        beta: between-cell-type
            diversity, weighted by mix     divergence
```

This is Jost's (2007) alpha/beta/gamma partitioning with cell types as the "sites" and the
cell fractions as their weights — and Jost showed **Shannon is the only standard diversity
measure that partitions this way when the weights are unequal**, which cell fractions always
are.

It gives a model-free version of the mediation question for alpha diversity: a change in
bulk diversity splits *exactly* into a part from changed within-cell-type diversity (cells
behaving differently) and a part from the between-cell-type term (the mixture). There is also
a hard mechanistic link: **JSD_f ≤ H(f)**, so the between-cell-type part is capped by the
diversity of the cell mix itself — a shift toward granulocyte dominance lowers the ceiling on
bulk transcriptome diversity.

Caveats: it needs per-sample cell-type expression (bMIND), which is noisy (the DE arm's
detection floors were 0.73 log2FC for granulocytes and 26 for B cells); the weights must be
RNA shares, which CIBERSORTx fractions approximate; and bMIND output must be on the linear
scale. The identity will hold only up to a residual, which must be reported.

---

## 4. How this maps onto our mediation framework

**The orientation matters, and most published methods have it backwards for us.**

The well-known distance-based mediation methods — MedTest (Zhang, Wei & Chen 2018), MODIMA
(Hamidi, Wallace & Alekseyenko 2019), PERMANOVA-med (Yue & Hu 2022) — were built for the
microbiome as the **mediator** (exposure → microbiome → outcome). Here the high-dimensional
data is the **outcome** and the mediator is low-dimensional:

```
lithium  ──►  cell mix (4 ILR balances)  ──►  transcriptome diversity
exposure      mediator                        outcome
```

| | consequence |
|:--|:--|
| **alpha** diversity is one number per person | the arm's validated estimator (`scripts/med.R`) applies **unchanged**, with a single outcome |
| **beta** diversity is a distance matrix | needs a distance-based adaptation (section 5.2). The published tools do not transfer directly |

MODIMA is additionally unsuitable as a primary test because it **cannot adjust for
confounders** (noted by Yue & Hu 2022), and age is imbalanced with lithium here (p = 0.0007).

---

## 5. The exact tests

### Summary table

| # | Question | Statistic | Null hypothesis | Inference |
|:--|:--|:--|:--|:--|
| **A1** | Does lithium change alpha diversity? | total effect TE of lithium on Shannon H | TE = 0 | t-test in `H ~ lithium + C`; bootstrap CI |
| **A2** | **Is that change mediated by cell mix?** | **natural indirect effect NIE** | **NIE = 0** | **percentile bootstrap 95 % CI over subjects (B = 10,000) + permutation test (B = 10,000)** |
| A3 | Does cell mix modify lithium's effect? | joint F on the 4 lithium × balance terms | no interaction | F-test, 4 df |
| **B1** | Does lithium shift whole-transcriptome location? | PERMANOVA pseudo-F for lithium given C | no location difference | 9,999 permutations, residuals permuted under the reduced model (Freedman & Lane 1983) |
| **B2** | Does lithium change heterogeneity? | PERMDISP F on distances to group centroids | equal dispersion | 9,999 permutations |
| **B3** | **Is the location shift mediated by cell mix?** | (i) lithium's pseudo-F given C **and** the ILR balances; (ii) share of lithium's adjusted R² that is shared with cell mix; (iii) joint significance of both paths | no remaining direct effect / no mediation | (i) partial dbRDA permutation test; (ii) bootstrap CI; (iii) both paths significant at α = 0.05 (MacKinnon et al. 2002) |
| B4 | Is a dispersion change mediated? *(only if B2 is significant)* | NIE with each person's distance to the non-user centroid as the outcome | NIE = 0 | as A2 |
| P | All of the above at the pathway level | as A and B, on Hallmark pathways | as above | as above |

### 5.1 Alpha diversity (A1–A3)

**The outcome, defined in advance:**

| Choice | Primary | Sensitivity | Why |
|:--|:--|:--|:--|
| measure | Shannon **H** (nats), reported as **exp(H)** | Hill q = 0 and q = 2 | H is additive (section 3.6) and close to normal; exp(H) is "effective number of genes" |
| scale | **TPM** proportions | read-count proportions | TPM approximates transcript proportions; counts match García-Nieto et al. |
| genes | **protein-coding**, the 12,368 filtered, **minus the 350 LM22 markers** | + markers; all genes | see section 7: structural non-coding RNAs dominate the counts; markers excluded for circularity, matching the arm's primary rule |
| covariates | the arm's set **+ log library size** | without library size | diversity is a known technical variable |

**The model.** Identical to the arm, with a single outcome:

```
M_k = b0_k + b1_k·lithium + b2_k'·C + e_k                         (k = 1..4)
H   = t0 + t1·lithium + Σ t2_k·M_k + Σ t3_k·(lithium·M_k) + t4'·C + e
NIE = Σ (t2_k + t3_k)·b1_k          NDE = t1 + Σ t3_k·b0_k          TE = NDE + NIE
```

Fitted by ordinary least squares (no voom weights: H is a per-person summary, not a count).

**The primary evidence for mediation is the test of NIE, not the proportion mediated.** The
arm already showed PM can sit near 0.4 under pure noise. PM is reported only if TE is
significant, and only as a description.

**Cost:** H is computed once per person from that person's own counts, so resampling only
reshuffles rows. 10,000 bootstrap replicates take seconds.

**Optional mechanistic layer:** the exact decomposition of section 3.6, using bMIND
estimates, to show how much of any change in H is within-cell-type versus mixture.

### 5.2 Beta diversity (B1–B4)

**Distance:** primary **Aitchison** (same gene set as alpha); secondary Jensen–Shannon and
Bray–Curtis.

**A useful simplification.** Aitchison distance is Euclidean, so distance-based redundancy
analysis on it is *exactly* ordinary multivariate regression on the principal components of
the log-ratio expression matrix. Every beta test can therefore be done with linear models,
and bootstrapping becomes resampling rows of that matrix — avoiding the zero distances that
duplicated subjects would otherwise create.

**B3 in more detail**, since it is the actual mediation question:

1. **Direct-effect test** — does lithium still explain transcriptome structure once the cell
   mix is conditioned on? (partial dbRDA, lithium given C + 4 ILR balances).
2. **Shared share** — variation partitioning (Borcard et al. 1992; Peres-Neto et al. 2006):
   `1 − R²adj(lithium | C, M) / R²adj(lithium | C)`, with a bootstrap CI. Read as "the
   fraction of lithium's association with whole-transcriptome structure that it shares with
   the cell mix". This is the multivariate *difference method*: associational, not a formal
   causal decomposition, and it should be described that way.
3. **Joint significance** — mediation declared only if *both* lithium → cell mix (MANOVA of
   the balances on lithium given C) *and* cell mix → transcriptome given lithium (partial
   dbRDA) are significant. This is the structure PERMANOVA-med uses, turned around to fit an
   outcome that is high-dimensional.
4. *(Secondary)* Principal-coordinate axes as outcomes: run `med.R` on each leading axis,
   giving an NIE and proportion per axis. PC1 of whole blood is very likely the
   neutrophil–lymphocyte axis, so expect strong mediation there — informative, but partly
   built in.

### 5.3 Pathway level (P)

1. Map genes to the 50 Hallmark sets.
2. **Alpha:** aggregate TPM by fractional assignment, proportions over pathways, Shannon H.
   Same model and tests as A1–A3.
3. **Beta:** GSVA scores per person, Euclidean distance, same tests as B1–B3.
4. Sensitivity: remove LM22 markers before aggregation and scoring.

### 5.4 Which exposure

| Exposure | Role |
|:--|:--|
| **lithium**, within bipolar I (152 vs 74) | **primary** — matches the arm's pre-specified primary model |
| bipolar off lithium vs controls (74 vs 234) | **exploratory** — a global measure might catch a diffuse signal the per-gene DE missed; carries the QC-filtered-controls caveat and needs assessment group as a covariate |

---

## 6. What we might expect (stated before looking)

- **Alpha:** neutrophils carry little RNA and a narrow, specialised programme, while
  monocytes and lymphocytes are transcriptionally richer. Lithium raises granulocytes relative
  to monocytes (b2, p = 0.012), so bulk diversity plausibly **falls** on lithium, with a large
  share mediated by the mix. This is a hypothesis, not a finding.
- **Beta location:** almost certainly significant, given 1,272 genes with a total effect;
  conditioning on cell mix should shrink it substantially.
- **Beta dispersion:** no prior expectation either way.

---

## 7. Pitfalls specific to this dataset — found by checking, not assumed

Two read-only checks were run (nothing written).

**1. Globin is not a problem.** Whole blood can be ~46 % globin reads without depletion.
Krebs et al. used TruSeq **Ribo-Zero Gold**, which depletes globin; the residual is median
**0.02 %** (max 2.5 %) and does not differ by lithium (Wilcoxon p = 0.55). Excluding globin
genes is a minor sensitivity only.

**2. Structural non-coding RNAs dominate instead — this one matters.** Ribo-Zero libraries
are total RNA, not poly-A selected, so abundant structural RNAs are captured:

| gene | share of all filtered counts |
|:--|--:|
| RN7SL1 (signal-recognition-particle RNA) | 6.86 % |
| RN7SL2 | 4.79 % |
| RN7SK | 3.45 % |
| MALAT1 | 2.33 % |
| RPPH1 | 1.14 % |
| **top 10 genes together** | **20.9 %** |

A diversity measure on all genes would largely measure *how much 7SL RNA each library
happened to capture* — a library-prep artefact. **Restrict to protein-coding genes.** The
GENCODE v19 annotation with gene biotypes is already on disk (two copies), so no download is
needed. Mitochondrially encoded genes should also be a sensitivity exclusion.

**3. Circularity** — the mediator was estimated from the same matrix. Exclude the LM22
markers (primary, as in the arm) and rerun with Krebs' independent deconvolution (as S2 did).

**4. Unbalanced design** → PERMDISP before interpreting PERMANOVA (section 3.4).

**5. The proportion mediated is not evidence of mediation** — established in the arm
(permutation p = 0.49). Evidence comes from the NIE test and joint significance.

---

## 8. What running it would require

**New R packages — none are installed, and each needs your approval** (CLAUDE.md: ask before
adding a dependency):

| package | source | needed for |
|:--|:--|:--|
| `vegan` | CRAN | PERMANOVA, dbRDA, PERMDISP, variation partitioning — the standard implementations |
| `GSVA` | Bioconductor | pathway scores (beta at pathway level) |
| `msigdbr` | CRAN | Hallmark gene sets |
| `energy` | CRAN | *optional*, only to compare against MODIMA |

The alpha analysis needs **no** new packages: Shannon and Hill numbers are one line each, and
the estimator already exists. With Aitchison distance, the beta tests can also be written as
linear models, so `vegan` could in principle be avoided — but it is the reference
implementation, and validating our own against it would need it anyway.

**Data on disk:** counts, gene lengths, the GENCODE v19 GTF, both ILR mediators, LM22. One
open item: the per-sample bMIND estimates needed for section 3.6 — `Boltz Paper/data/
processed/bmind_cts.rds` exists but came from the Boltz reproduction, so which cell types and
scale it holds must be confirmed before use.

**Compute:** small. Alpha: seconds. Beta: minutes. GSVA: minutes.

**Before running:** a short addendum to `PRESPEC.md` fixing the primary outcome, gene set,
distance and tests above — the same discipline as the main arm.

---

## 9. Open questions — RESOLVED 2026-10-18

Ali delegated these decisions ("make these decisions on your own … make sure it is
justifiable"). Each was decided and justified in **`PRESPEC.md` section 8.1**, in writing,
before any diversity analysis ran:

| Question | Decision |
|:--|:--|
| Which meaning of alpha/beta? | **All three**: (A) per-person / between-person = primary; (B) within- vs between-cell-type decomposition = secondary; (C) the cell mix's own diversity = secondary |
| Location, dispersion, or both? | **Both**; dispersion is mediated only if lithium changes dispersion |
| Lithium only, or bipolar vs control? | **Lithium primary; bipolar vs control exploratory** |
| Which pathway database? | **MSigDB Hallmark** primary, **KEGG legacy** sensitivity; aggregation by fractional counting rather than GSVA |

Two refinements made while implementing, both recorded:
- **Hill q = 0 (richness) is computed on all protein-coding genes, not the filtered set**, because
  the filter keeps only genes detected in ≥ 90 % of samples, which makes richness nearly constant
  by construction.
- **Path a of the beta joint-significance test uses MANOVA (Pillai)** exactly as declared here in
  section 5.2. It was briefly coded as an RDA, which weights balances by raw variance and so
  answers a different question; the RDA value is now reported beside the declared one.

Packages installed with approval: `vegan` 2.7.6 and `msigdbr` 26.1.1 plus 10 dependencies,
**22.8 MB**, and a **84.6 MB** MSigDB 2026.1.Hs gene-set cache. `GSVA` and `energy` were not needed.
The Hallmark and KEGG sets used are frozen in `data/genesets_msigdb_2026.1.Hs.rds`.

---

## 10. References

✓ = verified during this research (October 2026). Others are standard references cited from
the literature and should be checked against the source before going in the paper.

**Diversity foundations**
- Whittaker RH (1960). Vegetation of the Siskiyou Mountains, Oregon and California. *Ecological Monographs* 30:279–338.
- Whittaker RH (1972). Evolution and measurement of species diversity. *Taxon* 21:213–251.
- Shannon CE (1948). A mathematical theory of communication. *Bell System Technical Journal* 27:379–423.
- Hill MO (1973). Diversity and evenness: a unifying notation and its consequences. *Ecology* 54:427–432.
- Jost L (2006). Entropy and diversity. *Oikos* 113:363–375.
- ✓ Jost L (2007). Partitioning diversity into independent alpha and beta components. *Ecology* 88:2427–2439.
- Chao A et al. (2014). Rarefaction and extrapolation with Hill numbers. *Ecological Monographs* 84:45–67.
- Lin J (1991). Divergence measures based on the Shannon entropy. *IEEE Transactions on Information Theory* 37:145–151.

**Diversity of transcriptomes**
- ✓ Martínez O, Reyes-Valdés MH (2008). Defining diversity, specialization, and gene specificity in transcriptomes through information theory. *PNAS* 105(28):9709–9714.
- ✓ García-Nieto PE, Wang B, Fraser HB (2022). Transcriptome diversity is a systematic source of variation in RNA-sequencing data. *PLoS Computational Biology* 18(3):e1009939. *(The bioRxiv preprint was titled "…source of bias…".)*
- ✓ Ho JWK et al. (2008). Differential variability analysis of gene expression and its application to human diseases. *Bioinformatics* 24(13):i390–i398.
- Phipson B, Oshlack A (2014). DiffVar: a new method for detecting differential variability. *Genome Biology* 15:465.
- Zaneveld JR, McMinds R, Vega Thurber R (2017). Stress and stability: applying the Anna Karenina principle to animal microbiomes. *Nature Microbiology* 2:17121.

**Distance-based methods**
- Anderson MJ (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology* 26:32–46.
- Anderson MJ (2006). Distance-based tests for homogeneity of multivariate dispersions. *Biometrics* 62:245–253.
- ✓ Anderson MJ, Walsh DCI (2013). PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: what null hypothesis are you testing? *Ecological Monographs* 83:557–574.
- Legendre P, Anderson MJ (1999). Distance-based redundancy analysis. *Ecological Monographs* 69:1–24.
- McArdle BH, Anderson MJ (2001). Fitting multivariate models to community data. *Ecology* 82:290–297.
- ✓ Zapala MA, Schork NJ (2006). Multivariate regression analysis of distance matrices for testing associations between gene expression patterns and related variables. *PNAS* 103(51):19430–19435.
- ✓ Zapala MA, Schork NJ (2012). Statistical properties of multivariate distance matrix regression for high-dimensional data analysis. *Frontiers in Genetics* 3:190.
- Borcard D, Legendre P, Drapeau P (1992). Partialling out the spatial component of ecological variation. *Ecology* 73:1045–1055.
- Peres-Neto PR et al. (2006). Variation partitioning of species data matrices. *Ecology* 87:2614–2625.
- Freedman D, Lane D (1983). A nonstochastic interpretation of reported significance levels. *Journal of Business & Economic Statistics* 1:292–298.
- Quinn TP et al. (2018). Understanding sequencing data as compositions. *Bioinformatics* 34:2870–2878.

**Mediation with distance-based or high-dimensional data**
- ✓ Zhang J, Wei Z, Chen J (2018). A distance-based approach for testing the mediation effect of the human microbiome. *Bioinformatics* 34(11):1875–1883. *(MedTest)*
- ✓ Hamidi B, Wallace K, Alekseyenko AV (2019). MODIMA, a method for multivariate omnibus distance mediation analysis. *Genes* 10(7):524.
- ✓ Yue Y, Hu Y-J (2022). Extension of PERMANOVA to testing the mediation effect of the microbiome. *Genes* (PMC9222534). *(PERMANOVA-med; R package `LDM`, function `permanovaFL`)*
- MacKinnon DP et al. (2002). A comparison of methods to test mediation and other intervening variable effects. *Psychological Methods* 7:83–104.

**Pathway scoring**
- Barbie DA et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature* 462:108–112. *(ssGSEA)*
- ✓ Hänzelmann S, Castelo R, Guinney J (2013). GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics* 14:7.
- Liberzon A et al. (2015). The Molecular Signatures Database hallmark gene set collection. *Cell Systems* 1:417–425.

**Online sources consulted**
- García-Nieto et al. 2022 — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8982896/
- Martínez & Reyes-Valdés 2008 — https://pmc.ncbi.nlm.nih.gov/articles/PMC2443819
- MODIMA — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6679421/
- PERMANOVA-med — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9222534/
- `permanovaFL` documentation — https://search.r-project.org/CRAN/refmans/LDM/html/permanovaFL.html
- MedTest — https://researchwith.njit.edu/en/publications/a-distance-based-approach-for-testing-the-mediation-effect-of-the/
- MDMR (Zapala & Schork 2012) — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3461701/
- Jost 2007 appendix — https://esapubs.org/Archive/ecol/E088/145/appendix-A.htm
- GSVA — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3618321/
- Ho et al. 2008 — https://academic.oup.com/bioinformatics/article/24/13/i390/226572
- Anderson & Walsh (via PRIMER-e) — https://learninghub.primer-e.com/link/1028
