# Arm 3 — Live Result Ranking

> **GENERATED FILE — do not hand-edit.** Rebuild with `Rscript network/build_index.R`.
> Edit the per-finding files in `network/findings/` instead. On a merge conflict here,
> regenerate rather than merge. See `CONTRIBUTING.md`.

**Findings:** 6   **Updated:** 2026-10-14

Two rankings. **Table A (Evidence, 0–25)** governs the write-up; it rewards evidential
strength so a well-powered null can score top. **Table B (Interest, 0–10)** records how
striking a result is. **Δ = (Interest/10 − Evidence/25) × 10**: positive means the story
is ahead of the proof.

## Table A — Evidence

| # | Evidence | Interest | Δ | finding | components | source |
|--:|--:|--:|--:|:--|:--|:--|
| 1 | **22/25** | 7/10 | -1.8 | [The merge threshold is decisive: 33 modules vs 3](findings/2026-10-09-merge-threshold-decisive.md) | E=5 R=5 B=3 D=5 C=4 | network/runs/mod-000-resolution |
| 2 | **21/25** | 9/10 | +0.6 | [A lithium-associated co-expression dimension exists — and it is 73% cell composition](findings/2026-10-14-lithium-dimension-is-composition.md) |  | network/analyses/002-effective-dimensionality |
| 3 | **20/25** | 9/10 | +1.0 | [No module reliably beats a random gene set (degenerate decomposition)](findings/2026-10-09-gold-baseline-no-module-beats-random.md) | E=5 R=5 B=3 D=5 C=2 | network/runs/base-001 |
| 4 | **18/25** | 6/10 | -1.2 | [No module eigengene associates with lithium — and the modules are not independent](findings/2026-10-14-module-trait-no-lithium-association.md) |  | network/analyses/001-module-trait-ds4 |
| 5 | **15/25** | 7/10 | +1.0 | [Scale-free topology fails at every power and both sample sizes](findings/2026-10-09-scale-free-fails.md) | E=4 R=2 B=3 D=4 C=2 | network/runs/base-001 (sft curve) + Phase 0 probes |
| 6 | **6/25** | 1/10 | -1.4 | [Cost scaling, measured](findings/2026-10-09-cost-scaling.md) | E=4 R=0 B=1 D=1 C=0 | Phase 0 probes + base-001 draw timings |

## Table B — Interest

| # | Interest | Evidence | Δ | finding |
|--:|--:|--:|--:|:--|
| 1 | **9/10** | 21/25 | +0.6 | [A lithium-associated co-expression dimension exists — and it is 73% cell composition](findings/2026-10-14-lithium-dimension-is-composition.md) |
| 2 | **9/10** | 20/25 | +1.0 | [No module reliably beats a random gene set (degenerate decomposition)](findings/2026-10-09-gold-baseline-no-module-beats-random.md) |
| 3 | **7/10** | 22/25 | -1.8 | [The merge threshold is decisive: 33 modules vs 3](findings/2026-10-09-merge-threshold-decisive.md) |
| 4 | **7/10** | 15/25 | +1.0 | [Scale-free topology fails at every power and both sample sizes](findings/2026-10-09-scale-free-fails.md) |
| 5 | **6/10** | 18/25 | -1.2 | [No module eigengene associates with lithium — and the modules are not independent](findings/2026-10-14-module-trait-no-lithium-association.md) |
| 6 | **1/10** | 6/25 | -1.4 | [Cost scaling, measured](findings/2026-10-09-cost-scaling.md) |
