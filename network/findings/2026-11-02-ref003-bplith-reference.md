# ref-003: bp_lith (n=152) as reference — preservation in bp_nolith and control

**Date:** 2026-11-02
**Analysis:** `network/runs/ref-003` — 77 modules (+gold) defined once on all 152 bp_lith; tested in bp_nolith (n=74, one run, 50 perms) and control (20 draws of 74 from 234, 50 perms each). Power 12, bicor, signed.
**Evidence:** 18/25 (E=4 R=3 B=4 D=3 C=4) · **Interest:** 4/10 · **Δ:** −0.8

## Claim
Do modules defined in lithium-treated patients preserve in off-lithium patients and in healthy controls?

## Results (Zsummary is descriptive only — no p-values, no calibrated null)
- **bp_nolith:** module Zsummary median 31.1 (range 11.0–69.6); no module below 10. Gold Z=39.3, medianRank 78 (last of 78); 28/77 modules exceed gold Z.
- **control (20 draws):** module mean Zsummary median 27.2 (range 9.2–62.3); 2 modules mean <10. Gold Z mean 19.0 (sd 4.3, range 10.5–25.9), medianRank 77.75 (last); 57/77 modules exceed gold mean.
- Module Zsummary ordering agrees between the two test groups (Spearman 0.85 across 77 modules).
- Scale-free R² at power 12 is 0.71 (curve in `sft_r2_curve.csv`); not asserted as a fit.

## Interpretation and limits
- Descriptively, a lithium-patient reference network is preserved in both test groups, and gold ranks last in both, as in base-001. No evidence here that lithium patients have a distinctive network that others lack.
- **Not a lithium test:** Z is not comparable across test groups with different n/draw structure (nolith single run vs control mean of subsampled draws), and cell composition is not removed (finding #2).
- **Discrepancy recorded:** the script header and config note say 500 perms for bp_nolith, but the code reduced it to 50 (`NOLITH_PERMS`, METHODS.md 2026-10-26) and results show nperms=50. Nolith Z values therefore rest on 50 perms and have no across-draw spread.
- Draws revised 100 -> 20 for budget (breadth over precision), as in base-001. Control draws share samples and are not independent.
- The queue spec's ref-002 (BLOCKED) is not covered here.

<!-- SCORES: E=4 R=3 B=4 D=3 C=4 -->
