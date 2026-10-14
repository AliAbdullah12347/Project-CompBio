# No module eigengene associates with lithium — and the modules are not independent

**Date:** 2026-10-14
**Analysis:** `network/analyses/001-module-trait-ds4`
**Evidence:** 18/25 · **Interest:** 6/10 · **Δ:** −1.2

## Claim
At the usable decomposition (33 modules, `deepSplit=4, mergeCutHeight=0`), **no module
eigengene associates with lithium status in BP1 (n=226) after BH correction** — 0 of 33 in
the marginal, composition-adjusted and fully-adjusted models alike. Smallest raw p =
0.0055 (BH 0.138).

Cell composition (ILR b1–b4) attenuates the lithium coefficient by a **median 47.9%**,
so roughly half of the sub-threshold signal is compositional.

## Supporting
- 25 of 33 betas negative, 8 raw p < 0.05, 12 raw p < 0.10 — a weak consistent direction
  that does not survive correction.
- Only **1 of 33** modules correlates above |r| = 0.5 with any lineage fraction, so these
  modules are not simply cell-type proxies.
- **Median inter-eigengene |r| = 0.822.** The 33 modules are far from independent; BH over
  33 is mis-specified in the conservative direction.

## Scoring
E=4 (full model panel, three specifications, validated MTC library, deterministic — but
one decomposition and one sample). R=5 (this is the arm's question, asked directly).
B=2 (one decomposition, one parameter set, not yet replicated at other merge thresholds).
D=4 (a clear, correctable null with a quantified composition contribution).
C=3 (config, CSVs, README; no permutation null).

Interest 6: a null is the expected outcome for an exploratory arm, and the eigengene
inter-correlation is the more surprising part. Δ −1.2 — better supported than it is
striking.

## Caveats
Eigengene correlation is a weaker instrument than preservation for detecting *structural*
differences; `base-003-remerge` asks a different question at the same parameters and may
be more sensitive. The age confound is real but previously documented and adjusted for
(`de_analysis/METHODS.md` §6).
