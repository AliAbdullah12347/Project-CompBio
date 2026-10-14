# 001 — Module-trait relationships at the usable decomposition

**What was asked.** Does any module eigengene associate with lithium, and does any such
association survive adjustment for cell composition? This is the arm's question posed at
the module level rather than through preservation statistics.

**What was done.** Module assignments were taken from `mod-000-resolution`'s usable cell
(`deepSplit = 4, mergeCutHeight = 0`, 33 modules, largest 11.68% of genes) — *not* the
degenerate 3–4 module structure `base-001` used. Eigengenes were computed on the 226 BP1
samples only, since the lithium contrast is BP1-only (`IMAC_README` §5; `assessment group`
is constant there and was excluded as inestimable). Three models per module:
marginal, composition-adjusted (ILR b1–b4), and fully adjusted (+ age, sex, RIN, seqPC1–3).
BH within each trait; `de_v2/scripts/mtc.R` reused per §6.8.

**What came out.**

| model | modules with BH < 0.05 for lithium |
|:--|--:|
| marginal | **0 / 33** |
| composition-adjusted | **0 / 33** |
| fully adjusted | **0 / 33** |

Smallest raw p was 0.0055 (`darkturquoise`, BH 0.138). Eight modules reach raw p < 0.05,
twelve reach p < 0.10, and 25 of 33 betas are negative — a weak, broadly consistent
downward shift that does not survive correction.

**Composition explains about half of it.** Median attenuation of the lithium coefficient
after adjusting for ILR b1–b4 is **47.9%** (range 44–61% across the top modules). So the
already-sub-threshold signal is roughly halved by cell composition — directionally exactly
what this project exists to test, but on an effect too weak to call.

**Two things worth more than the headline.**

1. **The 33 eigengenes correlate at median |r| = 0.822** (mean inter-module r = 0.769).
   These are not 33 independent tests. The same density that made `mergeCutHeight = 0.25`
   cascade is still present at the eigengene level in the *usable* decomposition. Any
   multiple-testing correction that assumes 33 independent modules is mis-specified here,
   in the conservative direction.
2. **The modules are not cell-type proxies.** Only **1 of 33** correlates above |r| = 0.5
   with any lineage fraction (`darkturquoise`, r = 0.566). The intuition that whole-blood
   WGCNA modules are essentially lineage signals does not hold at this decomposition —
   which is also why composition adjustment attenuates but does not abolish the effect.

**What would change the answer.** A larger BP1 sample; a decomposition with genuinely
independent eigengenes (unclear that one exists for this data); or the permutation-based
preservation analysis at these parameters (`base-003-remerge`), which asks a different
question and may be more sensitive than eigengene correlation.

**Not a new finding:** lithium is confounded with age (on-lithium patients are 5.3 years
younger, p = 7.2e-4). `de_analysis/METHODS.md` §6 already documents this — *"Age is the
one real imbalance and is adjusted for"* — and reports it changes the DE answer little
(1,500 unadjusted vs 1,426 adjusted). The fully-adjusted model here agrees: still 0/33.

**Files.** `config.json`, `results/module_trait_marginal.csv` (33 modules x 16 traits),
`results/lithium_composition_adjusted.csv` (three models per module),
`results/module_eigengenes_bp1.rds`, `results/summary.rds`.
