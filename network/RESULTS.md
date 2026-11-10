# Arm 3 — Results

**Updated:** 2026-10-31T23:00Z
**Runs completed:** 13 of 51 queued (mod-000-resolution, inp-022, base-001, base-003-remerge, base-002-ceil, inp-023, pow-009, gene-027, inp-021, perm-031, **trait-033/034/035 DONE 2026-10-31T21:42Z**).
**Running:** none.
**Done 2026-11-01:** pres-030 — density>>connectivity preservation split is the same in bp_nolith and bp_lith (not lithium-specific); see findings/2026-10-31-pres030-density-connectivity-components.md. **Pending (next):** samp-043.
**Analyses completed:** 4 (001–004) + runs done. **Findings:** 34 — see `RANKING.md`.

---

## The answer so far

**The arm's central question is answered by base-002-ceil (DONE 2026-10-26T07:46Z).**

base-002-ceil scored one reference network (67 modules on 160 controls) against all three
test groups under matching conditions (same reference, same n=74, verified offset-removed
expression). Key result:

| arm | test group | gold Zsummary | gold rank | real modules > gold per draw |
|:--|:--|--:|:--|:--|
| base-002-ceil | ctrl_ceiling (74 disjoint ctrl) | 35.78 | 68/68 last | 15/67 (22%) |
| base-002-ceil | bp_nolith (n=74, 500 perms) | 23.16 | 68/68 last | 47/67 (70%) |
| base-002-ceil | bp_lith (20 draws, 50 perms) | 22.22 (sd 7.34) | 68/68 last, all 20 draws | 47.0/67 mean per draw |

**bp_lith and bp_nolith are nearly indistinguishable** (module Z correlation r=0.93;
mean module Z: nolith=31.05, lith=32.28). Both are below the within-controls ceiling
(gold Z 35.78). The ceiling can detect the between-population signal when it is present
(35.78 vs 23.16 in nolith) — but there is no meaningful nolith-vs-lith difference.

This is a **clean negative on the matched instrument.** The lithium-composition result
(73% of the lithium network dimension is cell composition; finding #2, E=24/25) is
therefore the final picture for this arm's central question.

**Summary table — all completed runs, within-arm gold comparison:**

| arm | reference | gold Z (nolith) | gold Z (lith, mean±sd) | gold rank |
|:--|:--|--:|--:|:--|
| base-001 | 71-mod raw | 14.01 | 21.74 ± 7.05 | always last |
| inp-022 | 44-mod ILR | 30.46 | 15.55 ± 8.13 | always last |
| base-002-ceil | 67-mod raw (n=160) | 23.16 | 22.22 ± 7.34 | always last |
| inp-023 | 42-mod tech+ILR | 28.08 | 16.98 ± 10.79 | always last |
| inp-021 | 63-mod technical | 31.13 | **DONE** — 20/20 lith draws; gold Z mean=32.60 sd=3.81 range=[26.71,45.00] | nolith: last; lith: last in all 20 |

Cross-arm Zsummary comparison remains invalid (reference size determines baseline).

**inp-023 DONE 2026-10-26T12:05Z:** Gold ranks last in ALL 20 draws (43/43); 40/42 modules
beat gold by mean Z (95%); 16/42 beat gold in every single draw. Replicates inp-022 pattern
under stricter (technical + ILR) residualisation. E=24/25. See
`findings/2026-10-25-inp023-technical-ilr-replication.md`.

**NOTE 2026-10-25 — raw matrix bug found (base-001 and base-003-remerge now re-run and
DONE):** base-001 lith draws, all of base-002-ceil, and base-003-remerge used raw
expression (`t(d$logtpm)`) for test matrices, bypassing `expr_for()`. Scripts fixed in
commit 1177d8c. All corrected. Do not quote the contaminated "base-001 lith gold Z=4.06"
figure — the verified gold Z is 21.74. See
`findings/2026-10-25-raw-test-matrices-bypassed-normalisation.md`.

---

## Completed runs

| id | wall | verdict | folder |
|:--|--:|:--|:--|
| `mod-000-resolution` | 1033 s | **Arm unblocked.** Rebuilt offset-removed: deepSplit=4, mergeCutHeight=0 → 71 modules | [network/runs/mod-000-resolution](runs/mod-000-resolution) |
| `inp-022` | 800 s | **KEY RESULT.** ILR b1–b4 residualisation: bp_lith gold Z=15.55. 41/44 (93%) modules above gold by mean Z; gold ranked last in all 20 draws. 44-module reference — 27 raw modules were composition-driven. E=13/25 (R raised 4→5 after inp-023 replication). | [network/runs/inp-022](runs/inp-022) |
| `base-001` | 3603 s | **VERIFIED.** 71-module reference, corrected test matrices. nolith gold Z=14.01 (72/72 last, 65/71 above); lith gold Z=21.74 (sd 7.05, 72/72 last in all 20 draws, 17/71 modules beat gold in all 20 draws). Gold's high Zsummary but last medianRank: size effect confirmed (cor 0.879). E=23/25. | [network/runs/base-001](runs/base-001) |
| `base-003-remerge` | 1 s | Read from base-001 v2 results; verified 20/20 draws complete. Reports consistent with base-001 log. No new finding. | [network/runs/base-003-remerge](runs/base-003-remerge) |
| `base-002-ceil` | 5357 s (attempt 4) | **DONE 2026-10-26T07:46Z.** Split-half ceiling, corrected matrices, all groups. ctrl_ceiling gold Z=35.78, nolith gold Z=23.16, lith mean gold Z=22.22 (sd 7.34). KEY NEGATIVE: lith≈nolith (r=0.93). E=24/25. | [network/runs/base-002-ceil](runs/base-002-ceil) |
| `inp-023` | 10311 s | **DONE 2026-10-26T12:05Z.** Technical+ILR residualisation: 42-module reference. nolith gold Z=28.08 (43/43 last, 25/42 above); lith gold mean Z=16.98 (sd 10.79, 43/43 last in all 20 draws, 40/42 beat gold by mean Z, 16/42 in ALL draws). Replicates inp-022. E=24/25. | [network/runs/inp-023](runs/inp-023) |
| `pow-009` | 86 s | **DONE 2026-10-26T14:12Z.** SFT curves on offset-removed expression, all 3 groups. Scale-free still fails (control max R²=0.770, bp_nolith=0.755, bp_lith=0.701; all below 0.80 threshold). Improvement from raw (0.66/0.41/0.54) confirms offset was part of the problem; residual failure is composition-driven (inp-022 SFT fits at power 9). E=22/25. | [network/runs/pow-009](runs/pow-009) |
| `gene-027` | 14093 s | **DONE 2026-10-26T23:15Z.** LM22 circularity check: 65-module reference on 12018 non-LM22 genes. nolith gold Z=12.23 (66/66 last); lith mean gold Z=21.41 (sd 6.67, range 9.00–33.02), 66/66 last in ALL 20 draws. 27/65 modules beat gold in all 20 draws. **Circularity check NEGATIVE** — same pattern as base-001 with LM22 genes excluded. E=21/25. | [network/runs/gene-027](runs/gene-027) |
| `inp-021` | 826 s (draw 20, 2nd attempt) | **DONE 2026-10-31T20:00Z.** Technical-only residualisation (age+sex+rin+plate+seqpc1-3): 63-module reference. nolith gold Z=31.13 (64/64 last); lith gold Z mean=32.60 (sd 3.81, range [26.71,45.00]), 64/64 last in ALL 20 draws. 10/63 modules beat gold in all 20 draws; 24/63 never beat gold. SFT max R²=0.797 (below 0.80; ILR required). 6 LM22-dependent modules absent. E=24/25. | [network/runs/inp-021](runs/inp-021) |
| `perm-031` | 4587 s (6 remaining draws) | **DONE 2026-10-31T21:25Z.** Group-label permutation null K=20: 74 random samples from all 474 vs base-001 bp_lith Zsummary. Gold null mean=32.02 (sd 5.21, range [17.94,39.02]); observed=21.74 < null in 19/20 draws. 44/72 modules less preserved in bp_lith than null mean (61%). **4 BH-significant (emp_p=0): blue, grey60, pink, purple.** Draw 18 null gold Z=17.94 < obs, so gold not BH-sig (padj=0.36). E=20/25. | [network/runs/perm-031](runs/perm-031) |
| `trait-033` | 13 s | **DONE 2026-10-31T21:31Z.** Eigengene ~ lithium + dx + age + sex + rin + plate, all 474, 71 modules. **29/71 lithium-significant** (BH <5%): 13 elevated, 16 reduced. Also: 46 RIN-sig, 52 plate-sig, 6 sex-sig, 4 dx-sig (blue, grey60, midnightblue, tan), 0 age-sig. Convergent with perm-031: blue, grey60, midnightblue overlap. E=24/25. See trait-034/035 for composition overlap. | [network/runs/trait-033](runs/trait-033) |
| `trait-034` | 15 s | **DONE 2026-10-31T21:34Z.** Eigengene ~ mono + T + NK + B (gran dropped), all 474, 71 modules. **211 BH-sig pairs** (T:56/71, NK:52/71, mono:54/71, B:49/71; gran:0). ALL 71 modules sig for any fraction. All 29 lithium-sig (trait-033) are fraction-sig: 0 lithium-only modules. E=22/25. | [network/runs/trait-034](runs/trait-034) |
| `trait-035` | 12 s | **DONE 2026-10-31T21:42Z.** Eigengene ~ ILR b1–b4, all 474, 71 modules. **217 BH-sig pairs** (b1:59, b2:50, b3:58, b4:50). ALL 71 modules ILR-significant. **CRITICAL: all 29 lithium-sig modules are ILR-sig → 0 lithium-specific modules.** Confirms 73% composition finding at module level. E=24/25; Interest=7/10. | [network/runs/trait-035](runs/trait-035) |

### Previously reset (now re-run and DONE)

`base-001` and `base-003-remerge` were reset 2026-10-25 (raw matrix bug in lith draws,
commit 1177d8c). Both re-ran successfully 2026-10-26 with corrected `expr_for()` test matrices.
Contaminated archived results remain in `_raw_test_contaminated_20261025T193947Z/`; do not use.

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

> **RETRACTED 2026-10-15.** Eigengenes were computed on the raw (library-offset-contaminated)
> expression matrix. The 31/33 result collapses to 0/33 once the library offset is removed.
> See finding `2026-10-14-bp-disease-drives-network-globally.md` for the full retraction note.
> Evidence now **0/25**. The 0/33 nolith-lith result (finding #001) is not affected.

`network/analyses/009-all-modules-across-groups/` · Finding:
[2026-10-14-bp-disease-drives-network-globally.md](findings/2026-10-14-bp-disease-drives-network-globally.md)

### 010 — PC1 of all-sample ME matrix: BP disease as a global network factor

> **RETRACTED 2026-10-15.** Same contamination as analysis 009 — eigengenes on raw expression.
> PC1=77.5% disease result is invalid. Evidence now **0/25**.
> See `2026-10-14-pc1-captures-bp-disease.md` for the full retraction note.

`network/analyses/010-full-pca-all-samples/` · Finding:
[2026-10-14-pc1-captures-bp-disease.md](findings/2026-10-14-pc1-captures-bp-disease.md)

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
13. ~~**BP disease (not lithium) drives a network-wide co-expression change** (31/33 modules,
    ctrl-nolith significant)~~ **RETRACTED 2026-10-15.** 31/33 collapses to 0/33 on
    offset-removed expression. The only valid sub-claim: 0/33 nolith-lith (finding #001).
14. ~~**A global network factor (PC1=77.5%)** captures the disease effect~~ **RETRACTED
    2026-10-15.** PC1 result was on raw-expression eigengenes. Evidence now 0/25.
15. **The lone exception (darkred) is the erythrocyte module**: hub gene ANK1 (ankyrin-1,
    kME=0.937), top GO erythrocyte development (FDR=2.1e-6). LM22 has no erythrocyte cell type
    → darkred has near-zero composition correlations (r_b2=−0.097) → invisible to the
    composition proxy → no group effect.
16. **BP disease tightens module coupling**: inter-module median |r| increases from 0.763
    (ctrl) to 0.790 (nolith) to 0.845 (lith). Mantel structure similarity ctrl vs nolith
    r=0.934 is significantly lower than the null (0.972; p=0.024). This tighter coupling is
    the structural explanation for PC1=77.5%: tightly correlated modules collapse onto one PC.

---

## Session 10 additions (2026-10-26T19:35Z)

### Infrastructure fix: cap audit

10 experiments had 5400s cap but run 500 nolith perms on the full ~12000-gene data set.
Measured nolith timing: ~6480s (from gene-027 first attempt). All 10 raised to 10800s:
pres-030, mod-013–mod-020, gene-029. Documented in METHODS.md 2026-10-26T19:20Z.

### gene-027 design adjustment: NOLITH_PERMS 500→50

A parallel session identified that 500 nolith perms at 65 modules takes ~6480s (exceeds old
7200s cap). The nolith perms were reduced to 50 and cap raised to 18000s. At 50 perms:
nolith takes ~720s (same as one lith draw). Circularity check still valid — we compare medianRank
(perm-independent) and relative Z patterns, not absolute Z values. Documented in gene-027.R
and METHODS.md 2026-10-26.

## Session 12 additions (2026-10-26T20:20Z)

### gene-027 preliminary results (RUNNING, 4/20 lith draws done as of 20:17Z)

- nolith: gold Z=12.23 (66/66 last) — matches base-001 pattern. Circularity check NEGATIVE.
- lith draws 1-4: gold last in ALL 4 (Z=30.16, 26.42, 21.89, 33.02; mean=27.87). Higher mean
  than base-001 (21.74) because fewer modules (65 vs 71) → less competition for gold rank.
- All 20 draws expected ~23:29Z, within 18000s cap. Finding template pre-written.
- Finding `network/findings/2026-10-26-gene027-lm22-circularity-preservation.md` —
  fill [FILL] sections once all 20 draws complete, then rebuild RANKING.md.

### Documentation cleanup

- Added retraction notes to analyses 009 and 010 in this document (above). Both were run on
  raw (library-offset-contaminated) eigengenes. 31/33 result → 0/33 on clean data. Headline
  chain items #13 and #14 struck.

## Next

1. **pres-030** — density vs connectivity decomposition of module preservation (~3h cap); queue position #1
2. **inp-024** — residualisation ladder next step (technical+ILR+interaction, if queued)
3. 38 more pending experiments cover gene subsets, module detection sensitivity, correlation
   function, network type, power, annotation, community, centrality, hub stability, bMIND.

**trait-033/034/035 DONE summary (2026-10-31T21:42Z):** CRITICAL FINDING: All 29 lithium-significant
modules (trait-033 BH 5%) are ILR-significant (trait-035); 0 lithium-only modules. ALL 71 modules
associate with ILR composition in both ILR (trait-035) and raw fraction (trait-034) approaches.
Confirmed at the module eigengene level: the 73% composition finding from analysis 002 extends to show
that every apparent lithium module association is captured by ILR. The joint model (lithium + ILR
simultaneously) is the decisive test; pending.

---

## The question

Do gene co-expression networks differ across healthy / bipolar-off-lithium /
bipolar-on-lithium, and does any difference survive removing cell composition?

Exploratory. A well-documented negative or inconclusive result is a complete deliverable.

---

---

---

---

---

---

## Cross-run synthesis (2026-11-10T07:54:55)

Generated by synth-048. Stability threshold: 80% of completed variants.

**Progress**: 46/46 experiments with results

| Experiment | Status |
|---|---|
| base-001 | complete |
| pres-030 | complete |
| inp-021 | complete |
| inp-022 | complete |
| inp-023 | complete |
| inp-024 | complete |
| gene-025 | complete |
| gene-026 | complete |
| gene-027 | complete |
| gene-028 | complete |
| corr-005 | complete |
| corr-006 | complete |
| net-007 | complete |
| net-008 | complete |
| pow-010 | complete |
| pow-011 | complete |
| ref-002 | complete |
| ref-004 | complete |
| pow-012 | complete |
| gene-029 | complete |
| mod-013 | complete |
| mod-014 | complete |
| mod-015 | complete |
| mod-016 | complete |
| mod-017 | complete |
| mod-018 | complete |
| mod-019 | complete |
| mod-020 | complete |
| perm-031 | complete |
| perm-032 | complete |
| ann-041 | complete |
| ann-042 | complete |
| trait-033 | complete |
| trait-034 | complete |
| trait-035 | complete |
| hub-036 | complete |
| hub-037 | complete |
| cent-040 | complete |
| comm-038 | complete |
| comm-039 | complete |
| samp-043 | complete |
| filt-044 | complete |
| tob-045 | complete |
| bmind-046 | complete |
| bmind-047 | complete |
| synth-048 | complete |

### Primary result (base-001)

bp_lith preservation (mean across 20 draws):

| Module | Zsummary mean (±SD) | medianRank mean | Interpretation |
|---|---|---|---|
| antiquewhite4 | 18.32 ± 2.73 | — | strongly preserved (Z≥10) |
| bisque4 | 23.15 ± 4.40 | — | strongly preserved (Z≥10) |
| black | 30.53 ± 3.24 | — | strongly preserved (Z≥10) |
| blue | 44.93 ± 3.79 | — | strongly preserved (Z≥10) |
| brown | 42.86 ± 5.33 | — | strongly preserved (Z≥10) |
| brown4 | 32.40 ± 4.47 | — | strongly preserved (Z≥10) |
| coral1 | 24.07 ± 2.25 | — | strongly preserved (Z≥10) |
| coral2 | 11.02 ± 1.90 | — | strongly preserved (Z≥10) |
| cyan | 79.93 ± 11.62 | — | strongly preserved (Z≥10) |
| darkgreen | 17.27 ± 2.56 | — | strongly preserved (Z≥10) |
| darkgrey | 42.30 ± 6.27 | — | strongly preserved (Z≥10) |
| darkmagenta | 38.82 ± 4.38 | — | strongly preserved (Z≥10) |
| darkolivegreen | 27.32 ± 4.50 | — | strongly preserved (Z≥10) |
| darkorange | 24.34 ± 3.32 | — | strongly preserved (Z≥10) |
| darkorange2 | 38.36 ± 3.62 | — | strongly preserved (Z≥10) |
| darkred | 50.53 ± 5.06 | — | strongly preserved (Z≥10) |
| darkseagreen4 | 20.71 ± 3.53 | — | strongly preserved (Z≥10) |
| darkslateblue | 21.29 ± 5.16 | — | strongly preserved (Z≥10) |
| darkturquoise | 27.35 ± 4.87 | — | strongly preserved (Z≥10) |
| floralwhite | 13.94 ± 1.68 | — | strongly preserved (Z≥10) |
| gold | 21.74 ± 7.05 | — | strongly preserved (Z≥10) |
| green | 39.47 ± 4.79 | — | strongly preserved (Z≥10) |
| greenyellow | 47.42 ± 4.10 | — | strongly preserved (Z≥10) |
| grey60 | 24.41 ± 2.87 | — | strongly preserved (Z≥10) |
| honeydew1 | 19.22 ± 3.74 | — | strongly preserved (Z≥10) |
| ivory | 35.05 ± 3.12 | — | strongly preserved (Z≥10) |
| lavenderblush3 | 17.95 ± 2.80 | — | strongly preserved (Z≥10) |
| lightcoral | 13.96 ± 1.96 | — | strongly preserved (Z≥10) |
| lightcyan | 50.65 ± 7.85 | — | strongly preserved (Z≥10) |
| lightcyan1 | 12.69 ± 1.57 | — | strongly preserved (Z≥10) |
| lightgreen | 55.23 ± 5.63 | — | strongly preserved (Z≥10) |
| lightpink4 | 34.22 ± 4.23 | — | strongly preserved (Z≥10) |
| lightsteelblue | 17.66 ± 2.03 | — | strongly preserved (Z≥10) |
| lightsteelblue1 | 22.91 ± 4.49 | — | strongly preserved (Z≥10) |
| lightyellow | 32.90 ± 4.10 | — | strongly preserved (Z≥10) |
| magenta | 35.77 ± 4.48 | — | strongly preserved (Z≥10) |
| maroon | 13.66 ± 1.47 | — | strongly preserved (Z≥10) |
| mediumorchid | 11.30 ± 1.70 | — | strongly preserved (Z≥10) |
| mediumpurple2 | 13.08 ± 1.50 | — | strongly preserved (Z≥10) |
| mediumpurple3 | 33.00 ± 5.06 | — | strongly preserved (Z≥10) |
| midnightblue | 53.64 ± 7.52 | — | strongly preserved (Z≥10) |
| navajowhite2 | 28.53 ± 3.17 | — | strongly preserved (Z≥10) |
| orange | 40.93 ± 5.56 | — | strongly preserved (Z≥10) |
| orangered3 | 16.54 ± 2.16 | — | strongly preserved (Z≥10) |
| orangered4 | 36.86 ± 4.45 | — | strongly preserved (Z≥10) |
| paleturquoise | 45.96 ± 5.04 | — | strongly preserved (Z≥10) |
| palevioletred3 | 23.32 ± 3.19 | — | strongly preserved (Z≥10) |
| pink | 22.04 ± 2.44 | — | strongly preserved (Z≥10) |
| plum | 15.21 ± 2.83 | — | strongly preserved (Z≥10) |
| plum1 | 26.32 ± 4.96 | — | strongly preserved (Z≥10) |
| plum2 | 30.59 ± 2.51 | — | strongly preserved (Z≥10) |
| purple | 49.80 ± 6.67 | — | strongly preserved (Z≥10) |
| red | 49.63 ± 6.84 | — | strongly preserved (Z≥10) |
| royalblue | 27.12 ± 2.94 | — | strongly preserved (Z≥10) |
| saddlebrown | 19.59 ± 3.13 | — | strongly preserved (Z≥10) |
| salmon | 23.39 ± 4.09 | — | strongly preserved (Z≥10) |
| salmon4 | 25.55 ± 3.71 | — | strongly preserved (Z≥10) |
| sienna3 | 27.02 ± 5.00 | — | strongly preserved (Z≥10) |
| skyblue | 29.98 ± 4.45 | — | strongly preserved (Z≥10) |
| skyblue1 | 14.55 ± 2.32 | — | strongly preserved (Z≥10) |
| skyblue2 | 8.64 ± 1.12 | — | weakly preserved (2≤Z<10) |
| skyblue3 | 26.44 ± 4.43 | — | strongly preserved (Z≥10) |
| steelblue | 34.91 ± 3.86 | — | strongly preserved (Z≥10) |
| tan | 42.19 ± 4.36 | — | strongly preserved (Z≥10) |
| thistle1 | 22.63 ± 3.41 | — | strongly preserved (Z≥10) |
| thistle2 | 33.17 ± 3.73 | — | strongly preserved (Z≥10) |
| turquoise | 60.97 ± 5.50 | — | strongly preserved (Z≥10) |
| violet | 33.29 ± 4.12 | — | strongly preserved (Z≥10) |
| white | 61.14 ± 6.30 | — | strongly preserved (Z≥10) |
| yellow | 54.12 ± 5.23 | — | strongly preserved (Z≥10) |
| yellow4 | 16.10 ± 2.00 | — | strongly preserved (Z≥10) |
| yellowgreen | 29.60 ± 6.06 | — | strongly preserved (Z≥10) |

bp_nolith preservation (Zsummary):

- antiquewhite4: Z=25.57 (strongly preserved (Z≥10))
- bisque4: Z=21.41 (strongly preserved (Z≥10))
- black: Z=34.76 (strongly preserved (Z≥10))
- blue: Z=48.09 (strongly preserved (Z≥10))
- brown: Z=38.59 (strongly preserved (Z≥10))
- brown4: Z=39.30 (strongly preserved (Z≥10))
- coral1: Z=20.80 (strongly preserved (Z≥10))
- coral2: Z=15.03 (strongly preserved (Z≥10))
- cyan: Z=75.31 (strongly preserved (Z≥10))
- darkgreen: Z=14.35 (strongly preserved (Z≥10))
- darkgrey: Z=32.76 (strongly preserved (Z≥10))
- darkmagenta: Z=36.74 (strongly preserved (Z≥10))
- darkolivegreen: Z=33.77 (strongly preserved (Z≥10))
- darkorange: Z=24.69 (strongly preserved (Z≥10))
- darkorange2: Z=35.91 (strongly preserved (Z≥10))
- darkred: Z=48.54 (strongly preserved (Z≥10))
- darkseagreen4: Z=23.07 (strongly preserved (Z≥10))
- darkslateblue: Z=22.76 (strongly preserved (Z≥10))
- darkturquoise: Z=29.05 (strongly preserved (Z≥10))
- floralwhite: Z=15.90 (strongly preserved (Z≥10))
- gold: Z=14.01 (strongly preserved (Z≥10))
- green: Z=41.04 (strongly preserved (Z≥10))
- greenyellow: Z=38.12 (strongly preserved (Z≥10))
- grey60: Z=26.33 (strongly preserved (Z≥10))
- honeydew1: Z=16.66 (strongly preserved (Z≥10))
- ivory: Z=35.39 (strongly preserved (Z≥10))
- lavenderblush3: Z=20.72 (strongly preserved (Z≥10))
- lightcoral: Z=13.26 (strongly preserved (Z≥10))
- lightcyan: Z=53.40 (strongly preserved (Z≥10))
- lightcyan1: Z=15.86 (strongly preserved (Z≥10))
- lightgreen: Z=44.54 (strongly preserved (Z≥10))
- lightpink4: Z=34.31 (strongly preserved (Z≥10))
- lightsteelblue: Z=18.45 (strongly preserved (Z≥10))
- lightsteelblue1: Z=22.06 (strongly preserved (Z≥10))
- lightyellow: Z=41.76 (strongly preserved (Z≥10))
- magenta: Z=41.09 (strongly preserved (Z≥10))
- maroon: Z=12.64 (strongly preserved (Z≥10))
- mediumorchid: Z=8.34 (weakly preserved (2≤Z<10))
- mediumpurple2: Z=10.61 (strongly preserved (Z≥10))
- mediumpurple3: Z=36.83 (strongly preserved (Z≥10))
- midnightblue: Z=59.78 (strongly preserved (Z≥10))
- navajowhite2: Z=28.16 (strongly preserved (Z≥10))
- orange: Z=31.03 (strongly preserved (Z≥10))
- orangered3: Z=17.36 (strongly preserved (Z≥10))
- orangered4: Z=34.50 (strongly preserved (Z≥10))
- paleturquoise: Z=44.08 (strongly preserved (Z≥10))
- palevioletred3: Z=26.56 (strongly preserved (Z≥10))
- pink: Z=23.61 (strongly preserved (Z≥10))
- plum: Z=17.76 (strongly preserved (Z≥10))
- plum1: Z=32.96 (strongly preserved (Z≥10))
- plum2: Z=33.81 (strongly preserved (Z≥10))
- purple: Z=58.60 (strongly preserved (Z≥10))
- red: Z=48.90 (strongly preserved (Z≥10))
- royalblue: Z=31.98 (strongly preserved (Z≥10))
- saddlebrown: Z=14.97 (strongly preserved (Z≥10))
- salmon: Z=23.02 (strongly preserved (Z≥10))
- salmon4: Z=31.21 (strongly preserved (Z≥10))
- sienna3: Z=27.73 (strongly preserved (Z≥10))
- skyblue: Z=31.66 (strongly preserved (Z≥10))
- skyblue1: Z=16.26 (strongly preserved (Z≥10))
- skyblue2: Z=10.53 (strongly preserved (Z≥10))
- skyblue3: Z=22.23 (strongly preserved (Z≥10))
- steelblue: Z=38.93 (strongly preserved (Z≥10))
- tan: Z=44.56 (strongly preserved (Z≥10))
- thistle1: Z=30.16 (strongly preserved (Z≥10))
- thistle2: Z=35.06 (strongly preserved (Z≥10))
- turquoise: Z=53.64 (strongly preserved (Z≥10))
- violet: Z=32.51 (strongly preserved (Z≥10))
- white: Z=57.59 (strongly preserved (Z≥10))
- yellow: Z=58.63 (strongly preserved (Z≥10))
- yellow4: Z=13.05 (strongly preserved (Z≥10))
- yellowgreen: Z=30.97 (strongly preserved (Z≥10))

### Stability across sensitivity variants

Stability threshold: ≥80% of 27 completed variants agree on Z-bracket.

| Module | Base-001 Z (mean) | Z-bracket | Variants | Agree | Stable |
|---|---|---|---|---|---|
| antiquewhite4 | 18.32 | Z≥10 | 14 | 14 | YES |
| bisque4 | 23.15 | Z≥10 | 17 | 17 | YES |
| black | 30.53 | Z≥10 | 27 | 27 | YES |
| blue | 44.93 | Z≥10 | 27 | 27 | YES |
| brown | 42.86 | Z≥10 | 27 | 27 | YES |
| brown4 | 32.40 | Z≥10 | 17 | 17 | YES |
| coral1 | 24.07 | Z≥10 | 14 | 14 | YES |
| coral2 | 11.02 | Z≥10 | 13 | 13 | YES |
| cyan | 79.93 | Z≥10 | 26 | 26 | YES |
| darkgreen | 17.27 | Z≥10 | 24 | 24 | YES |
| darkgrey | 42.30 | Z≥10 | 23 | 23 | YES |
| darkmagenta | 38.82 | Z≥10 | 21 | 21 | YES |
| darkolivegreen | 27.32 | Z≥10 | 21 | 21 | YES |
| darkorange | 24.34 | Z≥10 | 22 | 22 | YES |
| darkorange2 | 38.36 | Z≥10 | 17 | 17 | YES |
| darkred | 50.53 | Z≥10 | 24 | 24 | YES |
| darkseagreen4 | 20.71 | Z≥10 | 15 | 15 | YES |
| darkslateblue | 21.29 | Z≥10 | 17 | 17 | YES |
| darkturquoise | 27.35 | Z≥10 | 24 | 24 | YES |
| floralwhite | 13.94 | Z≥10 | 18 | 18 | YES |
| gold *(pseudo-module, not counted)* | 21.74 | Z≥10 | 27 | 27 | YES |
| green | 39.47 | Z≥10 | 27 | 27 | YES |
| greenyellow | 47.42 | Z≥10 | 26 | 26 | YES |
| grey60 | 24.41 | Z≥10 | 26 | 26 | YES |
| honeydew1 | 19.22 | Z≥10 | 15 | 15 | YES |
| ivory | 35.05 | Z≥10 | 18 | 18 | YES |
| lavenderblush3 | 17.95 | Z≥10 | 15 | 14 | YES |
| lightcoral | 13.96 | Z≥10 | 6 | 6 | YES |
| lightcyan | 50.65 | Z≥10 | 26 | 25 | YES |
| lightcyan1 | 12.69 | Z≥10 | 19 | 19 | YES |
| lightgreen | 55.23 | Z≥10 | 26 | 26 | YES |
| lightpink4 | 34.22 | Z≥10 | 15 | 15 | YES |
| lightsteelblue | 17.66 | Z≥10 | 6 | 6 | YES |
| lightsteelblue1 | 22.91 | Z≥10 | 20 | 20 | YES |
| lightyellow | 32.90 | Z≥10 | 26 | 25 | YES |
| magenta | 35.77 | Z≥10 | 27 | 26 | YES |
| maroon | 13.66 | Z≥10 | 15 | 15 | YES |
| mediumorchid | 11.30 | Z≥10 | 13 | 13 | YES |
| mediumpurple2 | 13.08 | Z≥10 | 6 | 6 | YES |
| mediumpurple3 | 33.00 | Z≥10 | 20 | 20 | YES |
| midnightblue | 53.64 | Z≥10 | 26 | 25 | YES |
| navajowhite2 | 28.53 | Z≥10 | 15 | 15 | YES |
| orange | 40.93 | Z≥10 | 22 | 22 | YES |
| orangered3 | 16.54 | Z≥10 | 6 | 6 | YES |
| orangered4 | 36.86 | Z≥10 | 20 | 20 | YES |
| paleturquoise | 45.96 | Z≥10 | 21 | 21 | YES |
| palevioletred3 | 23.32 | Z≥10 | 17 | 17 | YES |
| pink | 22.04 | Z≥10 | 27 | 27 | YES |
| plum | 15.21 | Z≥10 | 7 | 7 | YES |
| plum1 | 26.32 | Z≥10 | 21 | 19 | YES |
| plum2 | 30.59 | Z≥10 | 17 | 17 | YES |
| purple | 49.80 | Z≥10 | 26 | 26 | YES |
| red | 49.63 | Z≥10 | 27 | 27 | YES |
| royalblue | 27.12 | Z≥10 | 24 | 24 | YES |
| saddlebrown | 19.59 | Z≥10 | 22 | 22 | YES |
| salmon | 23.39 | Z≥10 | 26 | 26 | YES |
| salmon4 | 25.55 | Z≥10 | 17 | 17 | YES |
| sienna3 | 27.02 | Z≥10 | 21 | 21 | YES |
| skyblue | 29.98 | Z≥10 | 22 | 22 | YES |
| skyblue1 | 14.55 | Z≥10 | 8 | 8 | YES |
| skyblue2 | 8.64 | 2≤Z<10 | 10 | 1 | NO |
| skyblue3 | 26.44 | Z≥10 | 21 | 21 | YES |
| steelblue | 34.91 | Z≥10 | 22 | 22 | YES |
| tan | 42.19 | Z≥10 | 26 | 26 | YES |
| thistle1 | 22.63 | Z≥10 | 17 | 17 | YES |
| thistle2 | 33.17 | Z≥10 | 17 | 17 | YES |
| turquoise | 60.97 | Z≥10 | 27 | 27 | YES |
| violet | 33.29 | Z≥10 | 21 | 21 | YES |
| white | 61.14 | Z≥10 | 22 | 22 | YES |
| yellow | 54.12 | Z≥10 | 27 | 27 | YES |
| yellow4 | 16.10 | Z≥10 | 10 | 10 | YES |
| yellowgreen | 29.60 | Z≥10 | 21 | 21 | YES |

The 1 module(s) that fail: skyblue2 (base-001 Z = 8.64, 1 of 10 variants agree). The brackets are Z<2, 2<=Z<10 and Z>=10, so a module whose Z sits near 2 or near 10 can disagree with every variant while moving hardly at all. Check the distance to the boundary before reading any of these as unstable.

### Is the lithium network less preserved? (every variant)

Paired Wilcoxon on medianRank, bp_lith (mean over draws) against bp_nolith, over the modules each variant shares. **0 of 25 variants** reach p<0.05; the range is 0.441 to 0.989.

base-001 itself gives p = 0.444. base-003-remerge is base-001 copied into its own run directory and is excluded, or the primary run would be counted as agreeing with itself.

| Variant | Modules | lith medianRank | nolith medianRank | p |
|---|---:|---:|---:|---:|
| inp-021 | 64 | 30.7 | 31.5 | 0.441 |
| net-007 | 62 | 31.4 | 30.5 | 0.525 |
| pow-010 | 54 | 27.3 | 27.5 | 0.559 |
| net-008 | 68 | 36.3 | 35.0 | 0.589 |
| mod-019 | 39 | 19.6 | 19.0 | 0.593 |
| mod-013 | 25 | 11.6 | 13.0 | 0.662 |
| mod-014 | 54 | 27.6 | 26.5 | 0.668 |
| gene-026 | 31 | 15.3 | 15.0 | 0.688 |
| corr-005 | 73 | 34.6 | 36.0 | 0.707 |
| inp-023 | 43 | 19.9 | 22.0 | 0.709 |
| ref-004 | 82 | 42.3 | 41.0 | 0.726 |
| pow-011 | 64 | 32.5 | 31.0 | 0.753 |
| inp-022 | 45 | 22.0 | 22.0 | 0.778 |
| gene-025 | 20 | 10.6 | 10.5 | 0.827 |
| gene-027 | 66 | 33.4 | 32.5 | 0.863 |
| mod-020 | 20 | 10.6 | 10.5 | 0.863 |
| mod-016 | 66 | 34.2 | 32.5 | 0.883 |
| corr-006 | 64 | 31.9 | 31.5 | 0.894 |
| mod-015 | 72 | 37.6 | 36.0 | 0.927 |
| mod-018 | 24 | 12.6 | 13.0 | 0.933 |
| gene-029 | 60 | 29.1 | 30.5 | 0.937 |
| inp-024 | 67 | 35.5 | 34.0 | 0.947 |
| gene-028 | 10 | 6.0 | 6.0 | 0.961 |
| mod-017 | 42 | 19.6 | 22.5 | 0.978 |
| pow-012 | 72 | 36.3 | 36.0 | 0.989 |

Rows absent from the table ran a different comparison and have no results_nolith.csv to pair against: ref-002, and base-002-ceil, which is the matched-size ceiling rather than a parameter variant.

### Annotation (ann-041, ann-042)

ann-041 DEG overlap: 19 significant module×DEG-list pairs (BH<5%)
  - CT_LI: darkred,lightyellow,mediumpurple2,yellow
  - WB_LI: brown4,darkred,greenyellow,grey60,ivory,lightpink4,lightyellow,mediumpurple2,midnightblue,paleturquoise,palevioletred3,tan,thistle2,white,yellow

ann-042 LM22/CT-DE overlap: 14 significant pairs (BH<5%)

### Filter robustness (filt-044)

- antiquewhite4: 13% survive all 31 filters (worst filter: 13%)
- bisque4: 41% survive all 31 filters (worst filter: 41%)
- black: 13% survive all 31 filters (worst filter: 13%)
- blue: 4% survive all 31 filters (worst filter: 4%)
- brown: 6% survive all 31 filters (worst filter: 6%)
- brown4: 41% survive all 31 filters (worst filter: 41%)
- coral1: 52% survive all 31 filters (worst filter: 52%)
- coral2: 4% survive all 31 filters (worst filter: 4%)
- cyan: 8% survive all 31 filters (worst filter: 8%)
- darkgreen: 1% survive all 31 filters (worst filter: 1%)
- darkgrey: 53% survive all 31 filters (worst filter: 53%)
- darkmagenta: 46% survive all 31 filters (worst filter: 46%)
- darkolivegreen: 13% survive all 31 filters (worst filter: 13%)
- darkorange: 9% survive all 31 filters (worst filter: 9%)
- darkorange2: 65% survive all 31 filters (worst filter: 65%)
- darkred: 30% survive all 31 filters (worst filter: 30%)
- darkseagreen4: 32% survive all 31 filters (worst filter: 32%)
- darkslateblue: 0% survive all 31 filters (worst filter: 0%)
- darkturquoise: 0% survive all 31 filters (worst filter: 0%)
- floralwhite: 13% survive all 31 filters (worst filter: 13%)
- green: 2% survive all 31 filters (worst filter: 2%)
- greenyellow: 12% survive all 31 filters (worst filter: 12%)
- grey60: 2% survive all 31 filters (worst filter: 2%)
- honeydew1: 10% survive all 31 filters (worst filter: 10%)
- ivory: 46% survive all 31 filters (worst filter: 46%)
- lavenderblush3: 7% survive all 31 filters (worst filter: 7%)
- lightcoral: 7% survive all 31 filters (worst filter: 7%)
- lightcyan: 0% survive all 31 filters (worst filter: 0%)
- lightcyan1: 1% survive all 31 filters (worst filter: 1%)
- lightgreen: 24% survive all 31 filters (worst filter: 24%)
- lightpink4: 61% survive all 31 filters (worst filter: 61%)
- lightsteelblue: 27% survive all 31 filters (worst filter: 27%)
- lightsteelblue1: 8% survive all 31 filters (worst filter: 8%)
- lightyellow: 23% survive all 31 filters (worst filter: 23%)
- magenta: 0% survive all 31 filters (worst filter: 0%)
- maroon: 19% survive all 31 filters (worst filter: 19%)
- mediumorchid: 2% survive all 31 filters (worst filter: 2%)
- mediumpurple2: 0% survive all 31 filters (worst filter: 0%)
- mediumpurple3: 49% survive all 31 filters (worst filter: 49%)
- midnightblue: 26% survive all 31 filters (worst filter: 26%)
- navajowhite2: 33% survive all 31 filters (worst filter: 33%)
- orange: 26% survive all 31 filters (worst filter: 26%)
- orangered3: 17% survive all 31 filters (worst filter: 17%)
- orangered4: 31% survive all 31 filters (worst filter: 31%)
- paleturquoise: 56% survive all 31 filters (worst filter: 56%)
- palevioletred3: 35% survive all 31 filters (worst filter: 35%)
- pink: 0% survive all 31 filters (worst filter: 0%)
- plum: 25% survive all 31 filters (worst filter: 25%)
- plum1: 25% survive all 31 filters (worst filter: 25%)
- plum2: 36% survive all 31 filters (worst filter: 36%)
- purple: 30% survive all 31 filters (worst filter: 30%)
- red: 31% survive all 31 filters (worst filter: 31%)
- royalblue: 32% survive all 31 filters (worst filter: 32%)
- saddlebrown: 10% survive all 31 filters (worst filter: 10%)
- salmon: 4% survive all 31 filters (worst filter: 4%)
- salmon4: 4% survive all 31 filters (worst filter: 4%)
- sienna3: 7% survive all 31 filters (worst filter: 7%)
- skyblue: 20% survive all 31 filters (worst filter: 20%)
- skyblue1: 8% survive all 31 filters (worst filter: 8%)
- skyblue2: 2% survive all 31 filters (worst filter: 2%)
- skyblue3: 5% survive all 31 filters (worst filter: 5%)
- steelblue: 44% survive all 31 filters (worst filter: 44%)
- tan: 14% survive all 31 filters (worst filter: 14%)
- thistle1: 4% survive all 31 filters (worst filter: 4%)
- thistle2: 57% survive all 31 filters (worst filter: 57%)
- turquoise: 2% survive all 31 filters (worst filter: 2%)
- violet: 21% survive all 31 filters (worst filter: 21%)
- white: 72% survive all 31 filters (worst filter: 72%)
- yellow: 44% survive all 31 filters (worst filter: 44%)
- yellow4: 17% survive all 31 filters (worst filter: 17%)
- yellowgreen: 26% survive all 31 filters (worst filter: 26%)

### Tobacco sensitivity (tob-045)

Zsummary range across 5 tobacco imputations (sensitivity, 1 draw each):

- black: Z mean=53.91, range=28.29 (min=41.73, max=70.03)
- blue: Z mean=60.63, range=45.07 (min=45.79, max=90.87)
- brown: Z mean=53.64, range=22.22 (min=38.06, max=60.27)
- cyan: Z mean=38.34, range=45.38 (min=19.45, max=64.83)
- darkgreen: Z mean=28.12, range=22.26 (min=20.33, max=42.59)
- darkgrey: Z mean=36.95, range=24.49 (min=21.74, max=46.23)
- darkmagenta: Z mean=25.62, range=26.17 (min=14.06, max=40.24)
- darkolivegreen: Z mean=27.93, range=25.34 (min=13.93, max=39.26)
- darkorange: Z mean=26.37, range=48.46 (min=10.61, max=59.07)
- darkred: Z mean=41.27, range=39.04 (min=20.11, max=59.15)
- darkturquoise: Z mean=25.55, range=21.86 (min=16.30, max=38.17)
- gold: Z mean=13.92, range=8.42 (min=9.97, max=18.39)
- green: Z mean=48.94, range=13.92 (min=42.05, max=55.97)
- greenyellow: Z mean=47.91, range=34.74 (min=32.04, max=66.78)
- grey60: Z mean=43.24, range=41.34 (min=23.86, max=65.20)
- ivory: Z mean=14.74, range=10.82 (min=9.08, max=19.90)
- lightcyan: Z mean=44.86, range=46.72 (min=27.57, max=74.29)
- lightcyan1: Z mean=10.95, range=8.46 (min=7.94, max=16.40)
- lightgreen: Z mean=28.28, range=29.81 (min=18.95, max=48.76)
- lightsteelblue1: Z mean=14.61, range=11.97 (min=8.98, max=20.95)
- lightyellow: Z mean=28.69, range=21.88 (min=22.09, max=43.98)
- magenta: Z mean=46.03, range=21.74 (min=32.39, max=54.12)
- mediumpurple3: Z mean=14.98, range=11.47 (min=10.97, max=22.44)
- midnightblue: Z mean=29.12, range=30.76 (min=17.98, max=48.75)
- orange: Z mean=25.35, range=22.96 (min=15.47, max=38.42)
- orangered4: Z mean=21.24, range=29.77 (min=10.04, max=39.81)
- paleturquoise: Z mean=30.57, range=33.02 (min=19.60, max=52.62)
- pink: Z mean=31.05, range=14.89 (min=24.11, max=39.00)
- plum1: Z mean=20.85, range=17.43 (min=11.67, max=29.11)
- purple: Z mean=42.91, range=29.60 (min=26.08, max=55.68)
- red: Z mean=54.98, range=59.73 (min=22.64, max=82.38)
- royalblue: Z mean=21.93, range=9.37 (min=17.66, max=27.02)
- saddlebrown: Z mean=20.75, range=12.00 (min=13.48, max=25.48)
- salmon: Z mean=45.92, range=41.81 (min=17.48, max=59.29)
- sienna3: Z mean=27.14, range=27.31 (min=17.63, max=44.94)
- skyblue: Z mean=29.42, range=32.87 (min=15.77, max=48.63)
- skyblue3: Z mean=22.29, range=16.72 (min=16.55, max=33.28)
- steelblue: Z mean=21.03, range=12.06 (min=14.58, max=26.64)
- tan: Z mean=47.28, range=27.65 (min=32.24, max=59.89)
- turquoise: Z mean=49.95, range=36.65 (min=35.26, max=71.91)
- violet: Z mean=26.13, range=14.26 (min=19.48, max=33.75)
- white: Z mean=26.25, range=28.85 (min=16.45, max=45.31)
- yellow: Z mean=60.92, range=41.91 (min=37.99, max=79.90)
- yellowgreen: Z mean=20.74, range=22.67 (min=12.32, max=34.99)
- darkorange2: Z mean=15.48, range=3.59 (min=13.53, max=17.12)
- floralwhite: Z mean=14.72, range=6.10 (min=11.23, max=17.33)
- bisque4: Z mean=12.05, range=1.04 (min=11.53, max=12.57)
- brown4: Z mean=16.53, range=4.07 (min=14.49, max=18.56)
- darkslateblue: Z mean=14.86, range=0.00 (min=14.86, max=14.86)

### Per-lineage networks (bmind-047)

- B: 18 modules; Z range [5.38, 33.44]
- gran: 14 modules; Z range [4.07, 27.44]
- mono: 20 modules; Z range [1.64, 34.05]
- NK: 16 modules; Z range [3.76, 28.65]

### Conclusions

1. **Reference network**: 71 modules on whole-blood controls (bicor, signed, power 12, deepSplit 4, no merge), excluding WGCNA's gold pseudo-module.

2. **The lithium network is not less preserved.** 0 of 25 parameter variants reach p<0.05 on the primary contrast (range 0.441 to 0.989), and the primary run itself gives p = 0.444. Resolution, soft power, correlation function, network type, gene set and reference group were all varied; none of them changes the answer.

3. **No module-specific lithium effect.** 0 of 71 modules clear BH 5% against the group-label permutation null at K=157 draws; the smallest unadjusted p is 0.038. The four modules that led at K=20 were a floor artifact -- 1/(K+1) was 0.05 there, so nothing could clear BH whatever the biology.

4. **Zsummary is descriptive.** No p-value is derived from Zsummary itself: draws share the reference network and are not independent. Every p above comes from a permutation or a rank test, not from the Z scale.

5. **The exception lives in findings/, not here.** One module survives adjustment for cell composition -- see findings/2026-11-09-eosinophil-module-survives-composition.md. This row reads preservation results only and cannot see it.

