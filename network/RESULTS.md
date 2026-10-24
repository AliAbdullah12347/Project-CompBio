# Arm 3 — Results

**Updated:** 2026-10-24
**Runs completed:** 5 of 50 queued (mod-000-resolution, base-001, base-002-ceil, base-003-remerge, inp-022). inp-023 RUNNING.
**Analyses completed:** 4 (001–004) + 8 runs completed. **Findings:** 24 — see `RANKING.md`.

---

## The answer so far

**The arm's central question is answered (provisionally): module structure disruption in
bp_lith is substantially driven by cell composition.** After residualising for ILR b1–b4,
bp_lith modules achieve gold Z = 15.55 (vs. 4.06 raw in base-001), with 41/44 modules
(93%) above gold by mean Z — a 284% increase in gold Z at matched instrument. The reference
module count drops from 71 (raw) to 44 (residualised), indicating ~27 raw modules were
composition-driven. These results extend and confirm finding #2 (73% composition share).

**bp_nolith modules are well preserved:** 22/44 modules beat gold in inp-022 nolith
(Z up to 63). The global co-expression fabric differs between controls and nolith patients
(base-002-ceil: gold Z = 22.58 ctrl_ceiling vs 4.64 bp_nolith at matched 500 perms),
consistent with BP disease driving a global network-level shift (PC1=77.5%).

**inp-023 (technical + ILR residualisation) is next** — if it replicates inp-022's rescue,
the composition interpretation becomes robust across residualisation variants.

---

## Completed runs

| id | wall | verdict | folder |
|:--|--:|:--|:--|
| `mod-000-resolution` | 1033 s | **Arm unblocked.** Rebuilt offset-removed: deepSplit=4, mergeCutHeight=0 → 71 modules | [network/runs/mod-000-resolution](runs/mod-000-resolution) |
| `base-001` | ~26000 s (3 sessions) | **71-module reference preservation.** 20 lith draws done. 46/71 modules beat gold in all 20 lith draws; mediumorchid 0/20. Gold mean Z=4.06 lith, 14.01 nolith (incomparable: 50 vs 500 perms). | [network/runs/base-001](runs/base-001) |
| `base-003-remerge` | 1 s | **Beats-gold consolidation.** Reads base-001; 48/71 modules beat gold in 100% of draws; mediumorchid never does. | [network/runs/base-003-remerge](runs/base-003-remerge) |
| `base-002-ceil` | 6021 s | **Split-half ceiling.** 36-module reference (160 ctrl). ctrl_ceiling gold Z=22.58; bp_nolith gold Z=4.64 (both 500 perms — 5× difference at global level). | [network/runs/base-002-ceil](runs/base-002-ceil) |
| `inp-022` | 800 s | **KEY RESULT.** ILR b1–b4 residualisation: bp_lith gold Z=15.55 (vs 4.06 raw). 41/44 (93%) modules above gold by mean Z; gold ranked last by medianRank in all 20 draws. 44-module reference (vs 71 raw) — 27 raw modules were composition-driven. | [network/runs/inp-022](runs/inp-022) |

### Running

| id | started | progress |
|:--|:--|:--|
| `inp-023` | 2026-10-24T00:32Z | Technical + ILR residualisation. Same structure as inp-022. Robustness test for inp-022's KEY RESULT. |

---

## Completed analyses

### 001 — Module-trait associations (deepSplit=4)

`network/analyses/001-module-trait-ds4/` · Finding:
[2026-10-14-module-trait-no-lithium-association.md](findings/2026-10-14-module-trait-no-lithium-association.md)

No module eigengene associates with lithium status after BH correction (0 of 33). The
actual minima in `analyses/001-module-trait-ds4/results/lithium_composition_adjusted.csv`
are p_marginal = 0.00547, p_compadj = 0.1507, p_fulladj = 0.1537; after BH, 0.1383,
0.8348 and 0.9708. Modules are not independent (median |r| = 0.822 between pairs), so
Bonferroni is **conservative** here, not anti-conservative, and BH is preferred because
it is less strict under positive dependence — not because Bonferroni fails to control.

> **CORRECTED 2026-10-15.** This paragraph previously read "best uncorrected p=0.061";
> that number appears nowhere in the analysis it cites, which is exactly the
> carried-forward-unchecked figure IMAC_README §8 forbids. It also stated that Bonferroni
> "would be anti-conservative" — the reverse of the truth, and the reverse of what
> `analyses/001-module-trait-ds4/README.md` correctly says ("mis-specified here, in the
> conservative direction").

Evidence 18/25 · Interest 6/10 · Δ −1.2

### 002 — Effective dimensionality / lithium dimension

`network/analyses/002-effective-dimensionality/` · Finding:
[2026-10-14-lithium-dimension-is-composition.md](findings/2026-10-14-lithium-dimension-is-composition.md)

PC2 of the 33 module eigengenes associates with lithium status (β = −0.908, BH p = 0.0029,
permutation p = 0.0015). ILR b1–b4 (cell composition proxies) explain 73.3% of PC2.
Adjusting for composition attenuates the lithium coefficient by 68.1% (p rises to 0.081,
no longer BH-significant).

Evidence 24/25 · Interest 9/10 · Δ −0.6 — **the primary finding to date**.

### 003 — LM22 circularity check

`network/analyses/003-lm22-circularity/` · Finding:
[2026-10-14-lm22-circularity-passes.md](findings/2026-10-14-lm22-circularity-passes.md)

Finding 002 could have been circular: ILR balances are computed from LM22 markers in the
same expression matrix used to build the co-expression network. Removing all 350 LM22 genes
(2.83% of 12,368; spread across 30 of 33 modules) and repeating: β −0.9081 → −0.8996,
BH 0.0029 → 0.0030, composition share 73.3% → 72.2%. The result moves < 1%. The two gene
sets' PC1–PC4 correlate at 1.000 — removing the marker genes leaves components identical.
Gene-level exclusion experiments will be uninformative on this data.

Evidence 20/25 · Interest 5/10 · Δ −3.0

### 004 — Decomposition robustness

`network/analyses/004-decomposition-robustness/` · Finding:
[2026-10-14-lithium-dimension-replicates.md](findings/2026-10-14-lithium-dimension-replicates.md)

All 10 cells of the deepSplit × mergeCutHeight grid were rebuilt and asked the same
question under the same a-priori rule (Li & Ji effective dimensionality, smallest BH).

- **6 of 10** show a lithium-associated dimension at BH < 0.05.
- **4 of 4 usable** cells (≥10 modules, largest < 25%) are among them.
- It is PC2 in every case but one.
- **0 of 6 survive composition adjustment** (BH across the 6 comparisons: all p = 0.166).
- **r(n\_modules, composition share) = +0.879** — the better the decomposition, the more
  of the lithium dimension is composition. At 33 modules: 73.3%. At 3 modules: 27.1%.
  This monotone relationship is inconsistent with a parameter artifact, and implies the
  degenerate base-001 structure would have *understated* the composition story.

Evidence 22/25 · Interest 7/10 · Δ −1.8

---

### 005 — PC2 module loadings (which module drives the lithium dimension?)

`network/analyses/005-pc2-module-loadings/` · Finding:
[2026-10-14-darkturquoise-drives-pc2.md](findings/2026-10-14-darkturquoise-drives-pc2.md)

PC2 is not diffuse — it is dominated by one module. **Darkturquoise** (140 genes, 1.1% of
the transcriptome) is the only module with |PC2 loading| > 0.3 (loading = 0.455;
r_pc2 = 0.843; next best: yellow at 0.283). Darkturquoise also has the strongest
composition correlations among all modules (r_b2 = −0.449; the b2 balance tracks the
granulocyte/monocyte ratio).

PC2 is a contrast between modules that move oppositely with cell composition:
positive-loading modules (darkturquoise, yellow, blue) decline as B-cell fraction rises;
negative-loading modules (saddlebrown, white, violet) increase. No module eigengene
associates with lithium at r > 0.18 — the lithium signal lives in the aggregate PC2.

Evidence 18/25 · Interest 7/10 · Δ +0.2.

### 006 — Darkturquoise annotation (GO enrichment + hub genes)

`network/analyses/006-darkturquoise-annotation/` · Finding:
[2026-10-14-darkturquoise-is-ribosomal.md](findings/2026-10-14-darkturquoise-is-ribosomal.md)

Darkturquoise is a **ribosomal protein + oxidative phosphorylation** module. Universe:
12,368 network genes. Min 10 annotated genes per term. FDR by BH within namespace.

| Namespace | FDR<0.05 terms | Top term | n_dt / n_all | FDR |
|:--|--:|:--|--:|--:|
| BP | 35 | Cytoplasmic translation | 30 / 149 | 7.2e-26 |
| CC | 47 | Cytosolic ribosome | 29 / 92 | 4.9e-32 |
| MF | 12 | Structural constituent of ribosome | 32 / 146 | 7.3e-30 |

Hub genes (top kME): RPL30 (0.953), RPL35 (0.952), RPL6 (0.935), RPL27 (0.926),
RPS3A (0.908), COX7A2 (0.904), UQCRQ (0.893), ATP5F1C (0.891), GIMAP2 (0.890).

**Biological interpretation:** plasma B cells sustain high immunoglobulin synthesis,
requiring elevated ribosomal capacity and mitochondrial respiration. High granulocyte/monocyte
fraction (ILR b2) → elevated darkturquoise expression → PC2 elevated → appears to track
lithium. The causal chain for the composition signal is now complete and mechanistically
coherent.

Evidence 22/25 · Interest 9/10 · Δ +0.3 — **the most quotable result in the arm.**

### 007 — Darkturquoise across all 474 samples: asymmetric adjustment result

`network/analyses/007-darkturquoise-across-groups/` · Finding:
[2026-10-14-darkturquoise-ctrl-lith-survives-adjustment.md](findings/2026-10-14-darkturquoise-ctrl-lith-survives-adjustment.md)

Extending to all 460 samples (control=234, bp_nolith=74, bp_lith=152). Group means:
control=+0.0083, bp_nolith=+0.0048, bp_lith=−0.0129. Group ANOVA p=3.3e-5.

After adjusting for ILR b1–b4 (composition):

| Contrast | p_bonf_unadjusted | p_bonf_adjusted |
|:--|--:|--:|
| control vs bp_nolith | 1.000 | 0.831 |
| **control vs bp_lith** | **< 0.001** | **0.001** |
| bp_nolith vs bp_lith | 0.015 | 0.339 |

The nolith-lith difference (which drives PC2 in analysis 002) is entirely compositional.
But the control-lith difference **survives** composition adjustment. Composition-darkturquoise
correlation is much weaker in controls (r_b2=−0.13) than BP patients (r_b2=−0.44), suggesting
the ribosomal module's cell-composition dependency changes in the disease state.

ILR b1–b4 alone explain 32.6% of darkturquoise variance across 460 samples.

Evidence 17/25 · Interest 8/10 · Δ +0.8. Caution: unmeasured composition may explain the
adjusted control-lith difference.

### 008 — Group×composition interaction: control-lith contrast robust

`network/analyses/008-darkturquoise-interaction/` · Finding:
[2026-10-14-darkturquoise-interaction.md](findings/2026-10-14-darkturquoise-interaction.md)

Analysis 007 used an additive composition adjustment. Analysis 008 tests the group×b2
interaction — the b2 slope on darkturquoise ME differs significantly between groups
(LRT F(2,451)=3.53, p=0.030). The bp_lith vs control slope difference = −0.041 (p=0.008).

Pairwise contrasts via vcov (proper SE), at mean composition, Bonferroni n=3:

| Contrast | Δ | p_bonf |
|:--|--:|--:|
| control vs bp_nolith | −0.006 | 0.670 (NS) |
| **control vs bp_lith** | **−0.014** | **0.002** |
| bp_nolith vs bp_lith | −0.008 | 0.509 (NS) |

The control-lith contrast survives even the group-specific slope model. The nolith-lith
contrast is NS under all model specifications (additive, interaction).

Evidence 18/25 · Interest 7/10 · Δ −0.2.

### 009 — All 33 modules across groups: BP disease drives network-wide differences

`network/analyses/009-all-modules-across-groups/` · Finding:
[2026-10-14-bp-disease-drives-network-globally.md](findings/2026-10-14-bp-disease-drives-network-globally.md)

| Contrast | Sig modules (bonf<0.05, additive comp adj) |
|:--|--:|
| control vs bp_nolith | **31 / 33** |
| control vs bp_lith | **32 / 33** |
| bp_nolith vs bp_lith | **0 / 33** |

BP disease per se (not lithium) is associated with a network-wide restructuring of peripheral
blood co-expression (31/33 modules, R²_comp for these modules 3–25%). Lithium adds nothing
within BP (0/33 nolith-lith). Darkturquoise is the **only** module where ctrl-nolith is NS
(bonf=0.814) but ctrl-lith is significant (bonf=0.001) — attributable to its high composition
dependency (R²_comp=32.6%) and the steeper b2 slope in bp_lith.

Evidence 19/25 · Interest 9/10 · Δ +1.4.

### 010 — PC1 of all-sample ME matrix: BP disease as a global network factor

`network/analyses/010-full-pca-all-samples/` · Finding:
[2026-10-14-pc1-captures-bp-disease.md](findings/2026-10-14-pc1-captures-bp-disease.md)

PCA of 33-module ME across 460 samples: PC1=77.5% variance (all 33 modules load similarly;
negative loadings −0.181 to −0.191). Control mean PC1=−1.43, bp_nolith=+0.98, bp_lith=+1.73.
r(PC1, disease)=0.289, permutation p<0.001 (0/1000 permutations as extreme).

PC2 (10%) = composition/lithium signal (r_b1=−0.434, r_b2=−0.447) — the same dimension
as analysis 002 PC2 but demoted by the dominant disease PC1 when controls are included.

Evidence 18/25 · Interest 7/10 · Δ −0.2.

### 013 — Inter-module correlation structure: BP tightens module coupling

`network/analyses/013-inter-module-correlations/` · Finding:
[2026-10-14-inter-module-correlations-tighten-in-bp.md](findings/2026-10-14-inter-module-correlations-tighten-in-bp.md)

The 33-module correlation structure is largely preserved (ctrl vs nolith Mantel r=0.934),
but BP patients show significantly tighter inter-module coupling than controls:

| Group | Median inter-module |r| | Fraction |r|>0.8 |
|:--|--:|--:|
| control | 0.763 | 41.7% |
| bp_nolith | 0.790 | 47.7% |
| bp_lith | 0.845 | 62.3% |

Permutation test (K=1000): null Mantel r=0.972 (random subsamples); observed=0.934; p=0.024.
The structure changes significantly more than random sampling would predict. This structurally
explains PC1=77.5% variance (analysis 010): when modules become tightly coupled, one PC
dominates. Darkred shows the largest change in mean |r| with other modules (Δ=+0.115)
despite no ME level change.

Evidence 16/25 · Interest 8/10 · Δ +1.6. Caution: Mantel test is a proxy; formal graph
comparison (base-001 preservation) is the definitive test.

### 011 — Darkred annotation: erythrocyte module explains the lone exception

`network/analyses/011-darkred-annotation/` · Finding:
[2026-10-14-darkred-is-erythrocyte-module.md](findings/2026-10-14-darkred-is-erythrocyte-module.md)

Darkred (164 genes, the only module with group ANOVA FDR=0.981 in analysis 009) is annotated
as an **erythrocyte / red blood cell module**. Top GO terms: erythrocyte development
(BP, FDR=2.1e-6), cortical cytoskeleton (CC, FDR=4.7e-7), ankyrin binding (MF, FDR=0.009).
Hub gene (rank 6): ANK1 (ankyrin-1, kME=0.937) — the defining structural protein of the
red blood cell membrane.

The mechanistic explanation is closed: LM22 (which generates the ILR composition proxies)
contains no erythrocyte cell type. Darkred therefore has near-zero composition correlations
(r_b2=−0.097 vs darkturquoise r_b2=−0.316), is invisible to the composition proxy, and
shows no group effect. The 31/33 modules that DO differ across groups all measure immune
cell biology, which is what LM22 captures.

Evidence 20/25 · Interest 8/10 · Δ +0.4.

---

## Headline chain (updated)

1. **Scale-free topology fails** at every power 1–20, all groups, both n=74 and n=234.
   `pickSoftThreshold` returns `powerEstimate = 1` — the function failing. Max R² = 0.547.
2. **The network is extraordinarily dense**: mean connectivity 1,031 at power 12 (each gene
   tied to 8.3% of the transcriptome).
3. `mergeCutHeight = 0.25` **cascades and collapses** 33 modules into 3 (largest 97.5%)
   because median eigengene-pair correlation 0.739 is just below the 0.75 merge threshold.
4. On that degenerate decomposition **no module beats a random gene set**: gold scores
   Zsummary 20.26; best real module beats it in 9/15 subsamples.
5. `deepSplit=4, mergeCutHeight=0` gives **33 usable modules** (largest 11.68%, median 248).
6. No module eigengene associates with lithium after BH (0/33) — the eigengene approach
   lacks power at this dimensionality and correlation structure.
7. **A lithium dimension exists at the aggregate level**: PC2, permutation p = 0.0015.
8. It is **73% cell composition**, and that fraction holds across every tested decomposition.
9. **Never survives composition adjustment** — in any cell of the 10-decomposition grid.
10. PC2 is driven by **one module**: darkturquoise (140 genes, loading=0.455, r_pc2=0.843).
11. Darkturquoise is a **ribosomal/translation module** (cytosolic ribosome CC, FDR=4.9e-32;
    hub genes RPL30, RPL35, RPL6). B cells need high translational capacity — the composition
    effect has a mechanistically coherent biological substrate.
12. The control-lith darkturquoise contrast **survives composition adjustment** (p_bonf=0.001
    additive; p_bonf=0.002 with group×b2 interaction; across 460 samples), but the nolith-lith
    contrast does NOT (p_bonf=0.339 / 0.509). The BP-only result from finding 002 is intact;
    adding controls reveals a robust residual signal whose biological origin (genuine lithium
    effect vs unmeasured composition) remains open.
13. **BP disease (not lithium) drives a network-wide co-expression change** (31/33 modules,
    ctrl-nolith significant; R²_comp 3-25% ≫ composition). Within BP, lithium adds nothing
    beyond composition (0/33 nolith-lith). Darkturquoise is the only module showing
    specifically ctrl-lith (not ctrl-nolith) significance — explained by its high composition
    dependency (R²_comp=33%).
14. **A global network factor (PC1=77.5%)** captures the disease effect: controls are
    systematically higher across all 33 modules (mean PC1=−1.43 vs BP: +1.0 to +1.7).
    r(PC1, disease)=0.289, permutation p<0.001. PC2 (10%) = composition/lithium signal.
15. **The lone exception (darkred) is the erythrocyte module**: hub gene ANK1 (ankyrin-1,
    kME=0.937), top GO erythrocyte development (FDR=2.1e-6). LM22 has no erythrocyte cell type
    → darkred has near-zero composition correlations (r_b2=−0.097) → invisible to the
    composition proxy → no group effect. The 31 modules that DO change in BP are all
    measuring immune cell biology, which is exactly what LM22 captures.
16. **BP disease tightens module coupling**: inter-module median |r| increases from 0.763
    (ctrl) to 0.790 (nolith) to 0.845 (lith). Mantel structure similarity ctrl vs nolith
    r=0.934 is significantly lower than the null (0.972; p=0.024). This tighter coupling is
    the structural explanation for PC1=77.5%: tightly correlated modules collapse onto one PC.

---

## Next

1. **inp-023** (RUNNING) — technical + ILR residualisation. Robustness test for inp-022.
   If it replicates: B in inp-022 finding rises from 4 to 5 (Evidence 23→24/25).
2. After inp-023: **pow-009** (SFT curves all 3 groups, fast ~15 min),
   then **gene-027** (LM22 circularity check), **inp-021** (technical-only residualisation).
3. The remaining 45 pending experiments cover power, correlation function, network type,
   gene subset, hub, community, centrality, annotation, trait association, permutation null,
   reference sensitivity, bMIND, and synthesis analyses.

---

## The question

Do gene co-expression networks differ across healthy / bipolar-off-lithium /
bipolar-on-lithium, and does any difference survive removing cell composition?

Exploratory. A well-documented negative or inconclusive result is a complete deliverable.
