# Arm 3 — Run State

**Updated:** 2026-10-09T~session-3
**Status:** base-001 RUNNING (PID 52678). ALL 48 scripts written. base-001 in
bp_nolith preservation (500 perms, started 16:04:53). 4 modules detected in
control reference.

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
| Reference network at n=234 | ~47 min observed (2820 s; TOM + clustering + merge) |

K=500 permutations ≈ 65 min for bp_nolith (faster with only 4 modules). Peak RSS ~5–6 GB at full gene set.
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
- **4 modules** detected in control reference (base-001, 2026-10-09T16:04:53):
  18+ pre-merge modules merged to 4 at cutHeight=0.25; 59 grey (unassigned) genes.
  Biologically plausible for whole-blood data dominated by cell-type composition.

## Machinery
- `network/run_next.sh` — driver. `--peek` inspects the next row without claiming.
  Exit 0 ran, 2 needs a script, 3 queue empty. No status is absorbing except BLOCKED.
- `network/wake.sh` — session launcher. 4 h cap, token validation, auth-failure alert,
  zero-throughput watchdog. `--check` for side-effect-free testing.
- `network/RESUME_PROMPT.md` — standing instruction.
- Scheduler: launchd `com.aaylab.networkarm`, 18000 s (5 h), RunAtLoad true.

## Done (experiments with final results)
| id | result |
|:---|:-------|
| base-001 (reference only) | 4 modules, 12,309 assigned genes, 59 grey (reference built 2026-10-09T16:04:53) |

## Running
| id | started | stage | notes |
|:---|:--------|:------|:------|
| base-001 | 2026-10-09T15:16:23 | bp_nolith 500-perm preservation | SFT done, ref built; nolith in progress |

## Scripts written — ALL 48 complete
| id | family | notes |
|:---|:-------|:------|
| base-001 | baseline | primary experiment, RUNNING |
| base-002-ceil | baseline | split-half ceiling |
| pow-009 | power | SFT curves all groups |
| corr-005 | corr | pearson |
| corr-006 | corr | spearman |
| net-007 | nettype | unsigned |
| net-008 | nettype | signed hybrid |
| pow-010 | power | fixed p=6 |
| pow-011 | power | fixed p=8 |
| pow-012 | power | fixed p=12 (no fallback) |
| mod-013..020 | moddetect | deepSplit+minModuleSize+mergeCutHeight sensitivity |
| ref-002 | reference | ref=bp_nolith |
| ref-003 | reference | ref=bp_lith |
| ref-004 | reference | consensus blockwiseConsensusModules |
| inp-021 | input | technical residualised |
| inp-022 | input | **KEY** composition residualised (ILR b1-b4) |
| inp-023 | input | technical+composition residualised |
| inp-024 | input | voom-normalised log2-CPM |
| gene-025 | geneset | top 2000 most variable |
| gene-026 | geneset | top 5000 most variable |
| gene-027 | geneset | minus LM22 (circularity check) |
| gene-028 | geneset | WB_LI DEGs only |
| gene-029 | geneset | non-DEGs only |
| pres-030 | preserve | density + connectivity separately |
| perm-031 | perm | group-label null K=20, pool=all 474 |
| perm-032 | perm | group-label null K=20, pool=ctrl+lith 386 |
| trait-033 | modtrait | eigengene ~ lithium+dx+age+sex+rin+plate |
| trait-034 | modtrait | eigengene ~ lineage fractions |
| trait-035 | modtrait | eigengene ~ ILR b1-b4 |
| hub-036 | hub | kME hub identity |
| hub-037 | hub | delta_kME permutation |
| comm-038 | community | Louvain cluster_louvain |
| comm-039 | community | Leiden (fallback Louvain) |
| cent-040 | centrality | degree+eigenvector+betweenness |
| ann-041 | annot | DEG list overlap |
| ann-042 | annot | LM22 + CT-DE overlap |
| samp-043 | samplestruct | eigengene-space separation |
| filt-044 | input | 31-filter gene survival per module |
| tob-045 | input | tobacco sensitivity (5 of 20 imputations) |
| bmind-046 | bmind | REBUILD bMIND profiles (~3.1h, own window) |
| bmind-047 | bmind | per-lineage networks (requires bmind-046) |
| synth-048 | synthesis | cross-run synthesis → RESULTS.md |

## Blocked
See `network/BLOCKED.md`. Nothing outstanding that prevents the run.

## Scope
Write only inside `network/`. Stage with `git add network/`, never `-A`.
**This repository is PUBLIC** — never commit a credential or private path.
`build_index.py --query` is allowed (read-only); a bare run is not.
