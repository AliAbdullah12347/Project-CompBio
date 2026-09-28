#!/usr/bin/env python3
"""Impute the 30 missing tobacco values by stratified Beta-Bernoulli draw.

Method
------
For each of m = 20 completed datasets, and within each diagnosis stratum:

    p*  ~  Beta(1 + smokers, 1 + non-smokers)      posterior under a Jeffreys-ish
                                                   uniform prior on the rate
    t_i ~  Bernoulli(p*)                           for each subject missing tobacco

Drawing p* afresh for every dataset is what makes this *proper* multiple
imputation. Reusing a single point estimate of the rate would produce m datasets
that differ only by Bernoulli noise, understating the variance of the rate
itself and giving standard errors that are too small.

Why stratify on diagnosis and nothing else
------------------------------------------
Diagnosis is the only covariate in this dataset that predicts smoking. Measured
by 10x5-fold cross-validation on the 444 subjects with observed status:
diagnosis alone gives AUC 0.573, age and sex together give 0.491 -- worse than
chance -- and all nine demographic covariates together reach only 0.593. So a
richer imputation model would be a model to defend with nothing to show for it.

The prior also handles the one real numerical problem. Bipolar II has 0 smokers
in 13 observed subjects, and one of the 30 subjects needing imputation is
bipolar II. Beta(1, 14) gives that subject a rate near 6.7% -- small, non-zero,
and bounded. A logistic regression with diagnosis as a predictor instead sends
that coefficient past 17 and pins the subject at probability zero.

No expression data is read. There is consequently no circularity to defend:
imputed tobacco is not a function of the matrix that is the outcome in every
analysis where tobacco appears as a covariate.

Assumption
----------
MAR within diagnosis stratum. MCAR is already falsified -- all 30 missing
subjects are outside bipolar I. MAR itself is untestable, which is why
`--delta` runs the sensitivity analysis and why the assumption-free bounds are
reported next to every pooled estimate.

    python data/impute_tobacco.py
    python data/impute_tobacco.py --delta 0.5
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
COHORT = HERE / "cohort_474"
OUT = HERE / "imputed_tobacco"

# Move to config.yaml when the project scaffold is rebuilt; single source of truth.
SEED = 481
M = 20            # completed datasets. FMI for control statistics is ~29/234
                  # = 12%, and the usual rule sets m ~ 100 x FMI, so 20 is ample.


# ----------------------------------------------------------------- Student t
# Implemented here rather than pulled from scipy, which this project does not
# install. Verified against published values in _selftest() below.
def _betacf(a: float, b: float, x: float) -> float:
    TINY, EPS, MAXIT = 1e-300, 3e-16, 300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < TINY:
        d = TINY
    d = 1.0 / d
    h = d
    for i in range(1, MAXIT + 1):
        i2 = 2 * i
        aa = i * (b - i) * x / ((qam + i2) * (a + i2))
        d = 1.0 + aa * d
        c = 1.0 + aa / c
        if abs(d) < TINY:
            d = TINY
        if abs(c) < TINY:
            c = TINY
        d = 1.0 / d
        h *= d * c
        aa = -(a + i) * (qab + i) * x / ((a + i2) * (qap + i2))
        d = 1.0 + aa * d
        c = 1.0 + aa / c
        if abs(d) < TINY:
            d = TINY
        if abs(c) < TINY:
            c = TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < EPS:
            break
    return h


def _betainc(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lb = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(lb + a * math.log(x) + b * math.log1p(-x)) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lb + b * math.log1p(-x) + a * math.log(x)) * _betacf(b, a, 1 - x) / b


def t_ppf(q: float, df: float) -> float:
    """Quantile of Student's t by bisection on its CDF."""
    def cdf(t: float) -> float:
        p = 0.5 * _betainc(df / 2.0, 0.5, df / (df + t * t))
        return 1.0 - p if t > 0 else p
    lo, hi = -400.0, 400.0
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if cdf(mid) < q:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def _selftest() -> None:
    for df, want in ((10, 2.228), (30, 2.042), (1e7, 1.960)):
        got = t_ppf(0.975, df)
        assert abs(got - want) < 5e-4, f"t_ppf(.975,{df}) = {got}, expected {want}"


# ------------------------------------------------------------------ Rubin
def rubin(estimates: np.ndarray, variances: np.ndarray,
          df_complete: float) -> dict[str, float]:
    """Pool m estimates by Rubin's rules, with Barnard-Rubin degrees of freedom."""
    m = len(estimates)
    q_bar = float(estimates.mean())
    w_bar = float(variances.mean())                       # within-imputation
    b = float(estimates.var(ddof=1))                      # between-imputation
    t = w_bar + (1.0 + 1.0 / m) * b                       # total
    se = math.sqrt(t)

    # lambda: proportion of total variance attributable to the missing data.
    lam = ((1.0 + 1.0 / m) * b) / t if t > 0 else 0.0
    lam = min(max(lam, 1e-12), 1 - 1e-12)

    df_old = (m - 1) / lam ** 2
    df_obs = ((df_complete + 1.0) / (df_complete + 3.0)) * df_complete * (1.0 - lam)
    df = 1.0 / (1.0 / df_old + 1.0 / df_obs)

    crit = t_ppf(0.975, df)
    return {
        "estimate": q_bar, "se": se, "within_var": w_bar, "between_var": b,
        "lambda": lam, "fmi": (lam + 2.0 / (df + 3.0)) / (1.0 + lam), "df": df,
        "ci_low": q_bar - crit * se, "ci_high": q_bar + crit * se,
    }


# ------------------------------------------------------------------ imputation
def impute(meta: pd.DataFrame, rng: np.random.Generator, m: int = M,
           delta: float = 0.0) -> pd.DataFrame:
    """Return an n x m frame of completed tobacco columns.

    `delta` shifts the drawn log-odds for imputed subjects only, as an MNAR
    sensitivity analysis. delta = 0 is the MAR primary.
    """
    tob = meta["tobacco"].to_numpy(dtype=float)
    dx = meta["diagnosis"].to_numpy()
    missing = np.isnan(tob)

    draws = np.tile(tob, (m, 1))
    rates: list[dict] = []
    for stratum in pd.unique(dx):
        in_s = dx == stratum
        obs = in_s & ~missing
        need = in_s & missing
        smokers = float(tob[obs].sum())
        non = float(obs.sum() - smokers)
        p_star = rng.beta(1.0 + smokers, 1.0 + non, size=m)
        if delta:
            odds = p_star / (1.0 - p_star)
            p_star = (odds * math.exp(delta)) / (1.0 + odds * math.exp(delta))
        rates.append({"stratum": stratum, "observed_smokers": int(smokers),
                      "observed_non_smokers": int(non), "to_impute": int(need.sum()),
                      "observed_rate": smokers / max(obs.sum(), 1),
                      "mean_p_star": float(p_star.mean())})
        if need.any():
            draws[:, need] = (rng.random((m, int(need.sum())))
                              < p_star[:, None]).astype(float)

    out = pd.DataFrame(draws.T, index=meta["title"],
                       columns=[f"imp_{k+1:02d}" for k in range(m)])
    out.attrs["strata"] = pd.DataFrame(rates)
    return out


# ------------------------------------------------------------------ estimands
def control_rate(tob: np.ndarray, dx: np.ndarray) -> tuple[float, float]:
    """Control smoking prevalence, with its complete-data sampling variance."""
    y = tob[dx == "Control"]
    p, n = float(y.mean()), len(y)
    return p, p * (1.0 - p) / n


def bp1_vs_control_logor(tob: np.ndarray, dx: np.ndarray) -> tuple[float, float]:
    """log odds ratio for smoking, bipolar I against control (Haldane corrected)."""
    a = float(tob[dx == "BP1"].sum())
    b = float((dx == "BP1").sum() - a)
    c = float(tob[dx == "Control"].sum())
    d = float((dx == "Control").sum() - c)
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return math.log((a * d) / (b * c)), 1 / a + 1 / b + 1 / c + 1 / d


def bounds(meta: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """Assumption-free range over every possible assignment of the missing values."""
    tob = meta["tobacco"].to_numpy(dtype=float)
    dx = meta["diagnosis"].to_numpy()
    lo, hi = tob.copy(), tob.copy()
    lo[np.isnan(lo)] = 0.0
    hi[np.isnan(hi)] = 1.0
    return {
        "control_rate": (control_rate(lo, dx)[0], control_rate(hi, dx)[0]),
        # Filling every unknown control with a smoker raises the control rate,
        # which LOWERS the BP1-vs-control odds ratio; hence the swap.
        "bp1_vs_control_logor": (bp1_vs_control_logor(hi, dx)[0],
                                 bp1_vs_control_logor(lo, dx)[0]),
    }


ESTIMANDS = {"control_rate": control_rate,
             "bp1_vs_control_logor": bp1_vs_control_logor}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--delta", type=float, default=0.0,
                    help="MNAR log-odds shift applied to imputed subjects only")
    ap.add_argument("-m", type=int, default=M, help="number of completed datasets")
    args = ap.parse_args()

    _selftest()
    OUT.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(COHORT / "metadata_474.csv", dtype={"title": str})
    meta = meta.rename(columns={"bipolar disorder diagnosis": "diagnosis",
                                "tobacco use": "tobacco"})
    meta["tobacco"] = pd.to_numeric(meta["tobacco"], errors="coerce")
    n_missing = int(meta["tobacco"].isna().sum())
    print(f"cohort {len(meta)} samples, {n_missing} missing tobacco, "
          f"m = {args.m}, seed = {SEED}"
          + (f", MNAR delta = {args.delta:+.2f}" if args.delta else ""))

    rng = np.random.default_rng(SEED)
    imp = impute(meta, rng, m=args.m, delta=args.delta)

    print("\nstrata (Beta posterior per draw):")
    print(imp.attrs["strata"].to_string(index=False,
          float_format=lambda v: f"{v:.4f}"))

    dx = meta["diagnosis"].to_numpy()
    bnd = bounds(meta)
    df_complete = float(len(meta) - 2)

    print(f"\npooled by Rubin's rules over {args.m} datasets:")
    rows = []
    for name, fn in ESTIMANDS.items():
        qs, us = [], []
        for col in imp.columns:
            q, u = fn(imp[col].to_numpy(), dx)
            qs.append(q)
            us.append(u)
        r = rubin(np.array(qs), np.array(us), df_complete)
        r["estimand"] = name
        r["bound_low"], r["bound_high"] = bnd[name]
        rows.append(r)
        print(f"\n  {name}")
        print(f"    estimate      {r['estimate']:+.4f}   se {r['se']:.4f}")
        print(f"    95% CI        [{r['ci_low']:+.4f}, {r['ci_high']:+.4f}]   "
              f"df {r['df']:.1f}")
        print(f"    FMI           {r['fmi']:.3f}   "
              f"(within {r['within_var']:.2e}, between {r['between_var']:.2e})")
        print(f"    bound         [{r['bound_low']:+.4f}, {r['bound_high']:+.4f}]"
              "   assumption-free")

    # A sensitivity run must never clobber the primary. Suffix its outputs.
    tag = "" if args.delta == 0.0 else f"_delta{args.delta:+.3f}"
    imp.to_csv(OUT / f"tobacco_imputations{tag}.csv")
    pooled = pd.DataFrame(rows)[["estimand", "estimate", "se", "ci_low", "ci_high",
                                 "df", "fmi", "within_var", "between_var",
                                 "bound_low", "bound_high"]]
    pooled.to_csv(OUT / f"pooled_estimates{tag}.csv", index=False)

    imputed_only = imp.loc[meta.loc[meta["tobacco"].isna(), "title"]]
    per_subject = pd.DataFrame({
        "diagnosis": meta.set_index("title").loc[imputed_only.index, "diagnosis"],
        "n_smoker_draws": imputed_only.sum(axis=1).astype(int),
        "posterior_mean": imputed_only.mean(axis=1),
    })
    per_subject.to_csv(OUT / f"imputed_subjects{tag}.csv")

    print(f"\nthe {n_missing} imputed subjects, share of draws assigned smoker:")
    print(per_subject.groupby("diagnosis")["posterior_mean"]
          .agg(["count", "mean", "min", "max"]).round(3).to_string())
    print(f"\nwrote {OUT.name}/: tobacco_imputations.csv, pooled_estimates.csv, "
          "imputed_subjects.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
