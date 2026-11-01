# pres-030: density vs connectivity preservation components

**Date:** 2026-10-31
**Analysis:** `network/runs/pres-030` — modulePreservation density + connectivity components separately; 71 modules; bp_nolith (50 perms) + bp_lith (20 draws × 50 perms)
**Evidence:** 21/25 (E=5 R=4 B=3 D=4 C=5) · **Interest:** 6/10 · **Δ:** −2.4

*(COMPLETE — 20/20 lith draws, nolith 50 perms; finalised 2026-11-01 after a PARTIAL at the 10800 s cap, resumed from checkpoints.)*

## Claim

Does the bp_lith module preservation failure (base-001 finding) reflect loss of intra-module
co-expression density (genes still transcribed but decoupled), or loss of connectivity pattern
(hub gene structure disrupted), or both?

WGCNA's Zsummary is a composite of density-based and connectivity-based statistics. Decomposing
into these two components is scientifically informative:
- **Density preserved, connectivity failed**: intra-module co-expression survives but hub
  architecture is disrupted — consistent with a hub-specific regulatory change in bp_lith.
- **Density failed, connectivity preserved**: gene-gene correlations weaken uniformly but the
  hub structure partially survives — consistent with a global transcriptional noise increase.
- **Both fail equally**: Zsummary is a reliable summary; no additional information gained.

## Design

Same reference as base-001: 234 controls, bicor, signed, power 12, 71 modules.
- bp_nolith: 50 perms, full density + connectivity output
- bp_lith: 20 draws × 50 perms per draw, checkpoint per draw (`lith_draw_NNN.rds`)

Extracted per-draw statistics:
- `Zsummary.pres`: composite (same as base-001 Zsummary)
- `Zdensity.pres`: density component only
- `Zconnectivity.pres`: connectivity component only
- `medianRank.pres`: rank across modules (gold used as baseline)

## Results

**bp_nolith (50 perms):**

| component | gold Z | gold rank | real mean (71) | real sd | mods > gold |
|:--|--:|:--|--:|--:|--:|
| Zsummary | 30.82 | 72/72 (last) | 31.96 | 14.42 | 36/71 |
| Zdensity | 22.07 | — | **48.32** | 23.93 | **60/71** |
| Zconnectivity | 39.57 | — | **15.61** | 8.26 | **0/71** |

**STRIKING FINDING:** Density is well-preserved in bp_nolith (60/71 real modules > gold) but
connectivity is not (0/71 real modules above gold Zconnectivity of 39.57). Gold has HIGHER
Zconnectivity than all 71 real modules — the hub architecture of every real module is LESS
preserved than a random gene set's "connectivity."

70/71 real modules have Zdensity > Zconnectivity (only pink is a slight exception: −4.6).
Mean delta (density − connectivity) = +32.71. Real module intra-module co-expression density
survives in bp_nolith, but the specific hub gene structure does not.

**Top 10 modules: highest density − connectivity dissociation (nolith):**

| module | Zdensity | Zconnectivity | delta |
|:--|--:|--:|--:|
| cyan | 136.74 | 21.36 | +115.38 |
| midnightblue | 105.36 | 19.80 | +85.56 |
| lightcyan | 104.50 | 20.94 | +83.55 |
| white | 86.65 | 14.26 | +72.39 |
| purple | 86.82 | 22.35 | +64.47 |
| ivory | 71.62 | 13.81 | +57.81 |
| paleturquoise | 70.35 | 12.77 | +57.58 |
| darkred | 79.05 | 24.28 | +54.78 |
| brown4 | 64.24 | 10.14 | +54.10 |
| lightgreen | 74.17 | 21.04 | +53.13 |

Pink is the only module with connectivity slightly > density: Zdensity=21.90, Zconn=26.51, delta=−4.6.

**bp_lith (20/20 draws × 50 perms; per-draw `results_lith_alldraws.csv`):**

| quantity (gold = WGCNA random set) | mean over 20 draws | sd | range |
|:--|--:|--:|:--|
| gold Zsummary | 19.14 | 6.06 | 9.12–29.01 |
| gold Zdensity | 10.44 | 4.65 | 1.35–17.76 |
| gold Zconnectivity | 27.84 | 7.79 | 16.88–40.58 |

Per draw, real modules (71) above gold: Zsummary mean 53.6/71 (range 32–70); Zdensity 69.75/71
(65–71); Zconnectivity 8.0/71 (1–20). Real-module mean Zdensity ≈ 42–53, Zconnectivity ≈ 14–16 in
every draw. Zdensity > Zconnectivity for 70/71 modules (mean over draws), and for 95.8–100% of
modules in each draw (mean 98.3%). Per-module lith-vs-nolith correlation: density r=0.94,
connectivity r=0.97.

**Conclusion.** The density ≫ connectivity split is the same in bp_nolith and bp_lith. It is a property
of how real modules compare with the random gold set, not a lithium effect. What differs between
groups is only the gold level (lith lower on density), which is the base-001 pattern.

**QC flags (recorded, not adjusted):**
1. nolith gold Zsummary here is 30.82 (50 perms) vs 14.01 in base-001 (500 perms) on the same
   74 samples. The observed statistics are identical, so the difference is permutation-null noise. This
   run's nolith gold values are therefore unreliable. The lith draws are not affected the same way
   (gold Zsummary mean 19.14, sd 6.06, vs base-001 21.74, sd 7.05, consistent).
2. 50 perms per draw gives a noisy Z; only the descriptive pattern is claimed. No p-values, no
   calibrated null (README §6.3).
3. Gold Zconnectivity exceeding nearly all real modules is probably a feature of how Zconnectivity
   behaves for a random set of mixed-size genes; I did not test this. I make no hub-architecture claim.

## Connection to prior findings

- **[base-001](2026-10-26-base001-verified-clean.md)**: Zsummary gold Z=21.74 (sd 7.05); 17/71
  modules beat gold in all 20 draws. pres-030 Zsummary should replicate this result (same draws
  use same 50 perms). Any discrepancy is a QC flag.
- **[perm-031](2026-10-27-perm031-group-label-permutation-null.md)**: 4 BH-sig modules (blue,
  grey60, pink, purple) specifically less preserved in lith vs random null. pres-030 may show
  these modules fail on density, connectivity, or both.
- **[trait-033](2026-10-31-trait033-module-eigengene-trait-associations.md)**: blue and grey60
  are both preservation-failed (perm-031) and eigengene-altered (lith coefs). pres-030 adds
  the dimension/component breakdown.

<!-- SCORES PRE-SPECIFIED -->
<!-- E=5: full n=234 ctrl ref, WGCNA standard method, density+connectivity extracted together -->
<!-- R=5: directly tests WHICH component of preservation fails — a targeted mechanistic question -->
<!-- B=3: methodological extension of base-001; not independent replication; perm-031 is the null -->
<!-- D=5 if density/connectivity dissociation found; D=3 if pattern replicates Zsummary with no new info -->
<!-- C=5 if all 20 draws complete; C=2 if PARTIAL (~12 draws expected at cap, based on 832s/draw) -->
<!-- Interest=5/10 if same pattern; 7/10 if dissociation; 8/10 if dissociation is module-specific -->
<!-- Pre-specified PARTIAL (12 draws, dissociation): E=5+R=5+B=3+D=5+C=2 → 20/25; I=7/10; Δ=-3.0 -->
<!-- Pre-specified PARTIAL (12 draws, same): E=5+R=5+B=3+D=3+C=2 → 18/25; I=5/10; Δ=-2.8 -->

<!-- FINAL: pre-specified D=5 reduced to 4 (dissociation found but not group-specific); R=5 reduced to 4 (cannot attribute to hubs, see QC 3). E=5 R=4 B=3 D=4 C=5 = 21/25; I=6/10; delta=-2.4 -->
