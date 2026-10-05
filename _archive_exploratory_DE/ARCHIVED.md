# ARCHIVED — exploratory / preliminary work only

**Status: archived 30 September 2026. Not part of the confirmatory analysis.**

Everything in this folder is **exploratory**. It was produced by varying
thresholds, methods and models deliberately, and several analyses were chosen
*after* seeing earlier results. No p-value here is corrected for the number of
analyses run.

## What it may be used for

- **Preliminary results.** Reporting a finding as preliminary/exploratory is
  legitimate and is what this work is.
- Motivating the pre-specified analysis in `../de_analysis/`.
- Method development — the DE engine, the spike-in calibration design and the
  ILR machinery were all built here.

## What it must NOT be used for

- Confirmatory claims.
- Quoting a p-value or DEG count as though it were a planned test.
- Any hypothesis in `../de_analysis/config.yaml`. Those were written against a
  frozen protocol and committed before results existed.

The confirmatory analysis in `../de_analysis/` does not read this folder, does
not import from it, and does not cite its numbers.

See `FINDINGS.md` for what was found and `README.md` for how it was built.
