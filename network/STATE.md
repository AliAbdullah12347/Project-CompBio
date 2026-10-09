# Arm 3 — Run State

**Updated:** 2026-10-09T~24:30Z (session 2)
**Status:** base-001 RUNNING (PID 52678). Scripts written for next 6 rows. Awaiting
first preservation result before populating RESULTS.md numbers.

## Environment (all verified, not assumed)
- R 4.2.2, 4 cores, 16 GB RAM, ~189 GB free.
- All ten packages installed and loading: WGCNA 1.74, igraph 2.3.4, impute 1.72.3,
  preprocessCore 1.60.2, GO.db 3.16.0, AnnotationDbi 1.60.2, dynamicTreeCut 1.63.1,
  fastcluster 1.3.0, limma 3.54.2, edgeR 3.40.2. `enableWGCNAThreads()` → 3 workers.
- `Rscript network/00_load.R` prints the required line exactly:
  `12368 genes x 474 samples | control 234, bp_nolith 74, bp_lith 152`
- Headless CLI authenticates (`claude -p` → CLI_OK) using a one-year token in
  `~/.arm3_token` (0600, outside this public repository).
- `bash network/wake.sh --check` passes with no side effects.

## Measured costs (not estimated — see METHODS.md)
| stage | at 12,368 genes / n=74 |
|:--|:--|
| Reference TOM + blockwiseModules | ~786 s (~13 min) |
| Test TOM inside modulePreservation | ~786 s (recomputed every call — no API to inject) |
| modulePreservation, per permutation | ~7.8 s |
| Full draw (ref TOM + test TOM + 50 perms) | ~1,962 s (~33 min) |
| 20 draws (per row, lith only) | ~10.9 h |

K=500 permutations ≈ 65 min for bp_nolith. Peak RSS ~5–6 GB at full gene set.
**Never run two experiments concurrently.**

## Known scientific constraints (decided before any run — see METHODS.md)
- Scale-free topology **fails in this data**: R² never reaches 0.80 at any power in
  1:20, in any group, and is negative at low powers. `powerEstimate` returns 1
  (control, bp_nolith) or 2 (bp_lith) — the function failing, not a threshold.
  **Confirmed by base-001 scan: max signed-R² = 0.955 at power 2 (positive slope,
  not scale-free). Use power 12 (WGCNA FAQ signed n>40 table). See METHODS.md.**
- Draw count revised to 20 (from spec 100) for all `draws=100` rows. Per-draw cost
  ~1,962 s; 100 draws ≈ 54.5 h (10× oversubscribed). Decision in METHODS.md.
- Power 12 confirmed: earlier draft said 14 (wrong row of FAQ table; n=74 → n>40 row).
  METHODS.md second audit corrected this before any network was built.

## Machinery
- `network/run_next.sh` — driver. `--peek` inspects the next row without claiming.
  Exit 0 ran, 2 needs a script, 3 queue empty. No status is absorbing except BLOCKED.
- `network/wake.sh` — session launcher. 4 h cap, token validation, auth-failure alert,
  zero-throughput watchdog. `--check` for side-effect-free testing.
- `network/RESUME_PROMPT.md` — standing instruction.
- Scheduler: launchd `com.aaylab.networkarm`, 18000 s (5 h), RunAtLoad true.

## Done (experiments with final results)
*(none yet)*

## Running
| id | started | stage | notes |
|:---|:--------|:------|:------|
| base-001 | 2026-10-09T15:16:23 | reference network TOM (building) | SFT scan done; ref network ~50% |

## Scripts written, pending in queue
| id | family | script | notes |
|:---|:-------|:-------|:------|
| base-002-ceil | baseline | ✓ | split-half ceiling; needs base-001 first to calibrate |
| inp-021 | input | ✓ | technical residualised (age+sex+rin+plate+seqpc1-3) |
| inp-022 | input | ✓ | **KEY QUESTION** composition residualised (ILR b1-b4) |
| inp-023 | input | ✓ | **KEY QUESTION** technical+composition residualised |
| pow-009 | power | ✓ | scale-free R² curves all 3 groups (fast, ~15 min) |
| gene-027 | geneset | ✓ | circularity check: network minus 547 LM22 genes |

## Next unwritten scripts (by queue order)
| id | family | spec (truncated) |
|:---|:-------|:-----------------|
| perm-031 | perm | group-label permutation null K=500 |
| trait-033 | modtrait | eigengene ~ lithium, dx, age, sex, RIN, plate |
| trait-034 | modtrait | eigengene ~ 5 lineage fractions |
| trait-035 | modtrait | eigengene ~ ILR balances b1-b4 |
| pres-030 | preserve | density + connectivity components separately |
| samp-043 | samplestruct | eigengene-space separation by group |
| ann-041 | annot | module overlap with lithium + bipolar DEG lists |
| ann-042 | annot | module overlap with LM22 markers |
| ref-002 | reference | ref=bp_nolith (alternative reference) |
| ref-003 | reference | ref=bp_lith (alternative reference) |
| ref-004 | reference | ref=consensus(all3) |
| corr-005 | corr | pearson instead of bicor |
| corr-006 | corr | spearman instead of bicor |
| net-007 | nettype | unsigned network |
| net-008 | nettype | signed_hybrid network |
| pow-010..012 | power | fixed power 6, 8, 12 sensitivity |
| mod-013..016 | moddetect | deepSplit / minModuleSize sensitivity |

## Blocked
See `network/BLOCKED.md`. Nothing outstanding that prevents the run.

## Scope
Write only inside `network/`. Stage with `git add network/`, never `-A`.
**This repository is PUBLIC** — never commit a credential or private path.
`build_index.py --query` is allowed (read-only); a bare run is not.
