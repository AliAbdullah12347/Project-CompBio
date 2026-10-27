# Arm 3 — Results

**Updated:** 2026-10-26T23:30Z
**Runs completed:** 8 of 51 queued (mod-000-resolution, inp-022, base-001, base-003-remerge, base-002-ceil, inp-023, pow-009, gene-027).
**Running:** inp-021 (technical-only residualisation; next in queue).
**Pending:** inp-021 (next), then pres-030, perm-031, etc.
**Analyses completed:** 4 (001–004) + runs done. **Findings:** 32 — see `RANKING.md`.

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

1. **gene-027** (RUNNING at 19:20Z, 18000s cap, 4/20 draws) — When all 20 complete:
   - Run `Rscript network/scripts/analyze_gene027.R` — script exists
   - Fill finding file, score E=5 R=4 B=3 D=5 C=4 (21/25), Interest=6/10, Δ=-2.4
   - Rebuild RANKING.md and update this RESULTS.md (add gene-027 to completed table)
2. **inp-021** — technical-only residualisation (~2 sessions; starts automatically via LaunchDaemon)
3. **perm-031** — group-label permutation null (fast: 20 draws × 720s = 4h within 10800s cap)
4. **trait-033/034/035** — eigengene-trait associations (3600s cap, fast)
5. 40 more pending experiments cover gene subsets, module detection sensitivity, correlation
   function, network type, power, annotation, community, centrality, hub stability, bMIND.

---

## The question

Do gene co-expression networks differ across healthy / bipolar-off-lithium /
bipolar-on-lithium, and does any difference survive removing cell composition?

Exploratory. A well-documented negative or inconclusive result is a complete deliverable.
