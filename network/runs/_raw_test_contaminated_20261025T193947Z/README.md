# Archived 20261025T193947Z — raw test matrices scored against offset-removed references

32 experiment scripts built their TEST expression with `t(d$logtpm[, sids])`,
bypassing `expr_for()` and therefore the per-sample offset removal, while the
REFERENCE came from an offset-removed cache. 36 such uses were fixed.

- **base-002-ceil** — worst case: its reference was built from raw expression too
  (base-002-ceil.R:160), yet cached under a key ending `_offsetrm`. The cache file
  here is that mislabelled object. This is why its Zsummary tracked module size at
  r = +0.911 (the dense raw-graph signature) while bp_nolith showed -0.036.
- **base-001 lith draws** — reference and bp_nolith were clean (both via expr_for);
  only the 20 bp_lith test matrices were raw. nolith_preservation.rds was kept.
- **base-003-remerge** — pure derivative of base-001's mixed results.

Not valid as results. `network/lint_expr.sh` now fails any experiment that
reaches for d$logtpm as expression input.
