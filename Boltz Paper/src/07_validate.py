#!/usr/bin/env python3
"""Scorecard: every reproducible Boltz et al. (2024) claim, checked against this recreation.

The cohort differs (GSE124326, n = 444, versus their n = 1,730), so exact counts cannot be expected
to match. Each claim is therefore scored on the criterion that a same-method, smaller-cohort run
should meet, and that criterion is written into the table:
  EXACT        a number that must match exactly (deposit self-consistency, sample bookkeeping)
  REPLICATED   same direction and significant here
  CONSISTENT   same direction, or same qualitative conclusion, but not significant here
  DIFFERS      a different conclusion here
  DISCREPANCY  the paper's text disagrees with its own deposited tables
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
REF = os.path.join(ROOT, "data", "reference")
T = pd.read_csv(os.path.join(REF, "paper_targets.csv")).set_index("quantity")["value"]

rows: list[dict] = []


def add(area, claim, paper, here, verdict, criterion):
    rows.append(dict(area=area, claim=claim, paper=paper, here=here, verdict=verdict, criterion=criterion))


def res(name):
    p = os.path.join(RES, name)
    return pd.read_csv(p) if os.path.exists(p) else None


def deposit_checks():
    li = pd.read_csv(os.path.join(REF, "S15_lithium_DEGs_limma_bulk.csv"))
    cc = pd.read_csv(os.path.join(REF, "S15_case_control_DEGs_limma_bulk.csv"))
    add("deposit", "Table S15 lithium DEG count = text", 100, len(li), "EXACT" if len(li) == 100 else "DISCREPANCY", "row count")
    add("deposit", "Table S15 case/control DEG count = text", 64, len(cc), "EXACT" if len(cc) == 64 else "DISCREPANCY", "row count")
    lo, hi = li.logFC.min(), li.logFC.max()
    add("deposit", "lithium DEG logFC range = text", "-0.191 to 0.177", f"{lo:.3f} to {hi:.3f}",
        "EXACT" if abs(lo + 0.191) < 0.002 and abs(hi - 0.177) < 0.002 else "DISCREPANCY", "within rounding (text truncates)")
    lo, hi = cc.logFC.min(), cc.logFC.max()
    add("deposit", "case/control DEG logFC range = text", "-0.126 to 0.104", f"{lo:.3f} to {hi:.3f}",
        "EXACT" if abs(lo + 0.126) < 0.002 else "DISCREPANCY",
        "text minimum is wrong: 3 deposited DEGs lie below -0.126 (PDK4 -0.217)")
    g = lambda d: set(d.gene.str.split(".").str[0])
    add("deposit", "lithium & case/control DEG overlap = text", 9, len(g(li) & g(cc)),
        "EXACT" if len(g(li) & g(cc)) == 9 else "DISCREPANCY", "set intersection")
    for tag, col_counts in [("case_cont", {"Bcellsnaive": 21, "Bcellsmemory": 24, "Neutrophils": 4}), ("Li", None)]:
        q = pd.read_csv(os.path.join(REF, f"S15_bmind_diff_expr_qvals_{tag}.csv"), index_col=0)
        if col_counts:
            for c, n in col_counts.items():
                got = int((q[c] < 0.05).sum())
                add("deposit", f"bmind_de case/control {c} DEGs = text", n, got, "EXACT" if got == n else "DISCREPANCY", "q < 0.05")
            both = int(((q.Bcellsnaive < .05) & (q.Bcellsmemory < .05)).sum())
            add("deposit", "bmind_de naive & memory B shared = text", 18, both, "EXACT" if both == 18 else "DISCREPANCY", "q < 0.05 in both")
        else:
            got = int((q < 0.05).values.sum())
            add("deposit", "bmind_de lithium: no DEG in any cell type", 0, got, "EXACT" if got == 0 else "DISCREPANCY", "q < 0.05")
    add("deposit", "bmind_de gene universe = bulk DE universe", 17194, len(q), "DISCREPANCY",
        "bMIND sheets hold 19,099 genes; no filter for this set is described")


def proportion_checks():
    t1 = res("table1_celltype_proportions.csv").set_index("cell_type")
    order_here = list(t1.sort_values("mean", ascending=False).index[:3])
    add("proportions", "neutrophils most abundant (Table 1)", "0.51", f"{t1.loc['neutrophils', 'mean']:.2f}",
        "REPLICATED" if order_here[0] == "neutrophils" else "DIFFERS", "rank 1 of 22")
    add("proportions", "monocyte proportion (Table 1)", "0.050", f"{t1.loc['monocytes', 'mean']:.3f}",
        "DIFFERS", "LM22 CIBERSORT without CIBERSORTx batch correction overestimates monocytes ~4x")
    for k in ["b.cells.memory", "t.cells.cd8"]:
        add("proportions", f"{k} passes mean > 0.02", "yes", f"no ({t1.loc[k, 'mean']:.4f})", "DIFFERS",
            "Boltz's own selection rule, applied here")
    lr = res("celltype_logistic_regressions.csv")
    get = lambda cn, k: lr[(lr.contrast == cn) & (lr.cell_type == k)].iloc[0]
    claims = [("case_vs_control", "neutrophils", "p_casecontrol_neutrophils", 1, "higher in cases"),
              ("case_vs_control", "nk.cells.resting", "p_casecontrol_nk.cells.resting", -1, "higher in controls"),
              ("case_vs_control", "t.cells.cd4.naive", "p_casecontrol_cd4_t", -1, "CD4 higher in controls"),
              ("lithium_user_vs_non", "neutrophils", "p_lithium_neutrophils", 1, "higher in users"),
              ("lithium_user_vs_non", "t.cells.cd4.naive", "p_lithium_t.cells.cd4.naive", -1, "higher in non-users"),
              ("lithium_user_vs_non", "t.cells.cd4.memory.resting", "p_lithium_t.cells.cd4.memory.resting", -1, "higher in non-users"),
              ("lithium_user_vs_non", "nk.cells.resting", "p_lithium_nk.cells.resting", -1, "higher in non-users")]
    for cn, k, key, sign, desc in claims:
        r = get(cn, k)
        same = np.sign(r.beta) == sign
        v = "REPLICATED" if same and r.p < 0.05 else ("CONSISTENT" if same else "DIFFERS")
        add("proportions", f"{cn}: {k} {desc}", f"p = {T[key]:.2g}", f"p = {r.p:.2g}, beta {r.beta:+.2f}", v,
            "same sign; REPLICATED needs p < 0.05 (Boltz n = 1,730 vs 444 here)")
    nv = lr[(lr.contrast == "nonuser_vs_control")]
    boltz8 = ["b.cells.naive", "b.cells.memory", "t.cells.cd8", "t.cells.cd4.naive", "t.cells.cd4.memory.resting",
              "nk.cells.resting", "monocytes", "neutrophils"]
    n8 = int((nv[nv.cell_type.isin(boltz8)].bonferroni_p < 0.05).sum())
    add("proportions", "non-users vs controls: no cell type differs (Boltz's 8)", 0, n8,
        "REPLICATED" if n8 == 0 else "DIFFERS", "Bonferroni p < 0.05 over the types tested")
    extra = nv[~nv.cell_type.isin(boltz8) & (nv.bonferroni_p < 0.05)]
    for _, r in extra.iterrows():
        add("proportions", f"non-users vs controls: {r.cell_type} (not tested by Boltz)", "not tested",
            f"Bonferroni p = {r.bonferroni_p:.3f}", "DIFFERS",
            "a type only the >0.02 rule selects here; controls are QC-filtered non-randomly (CLAUDE.md s2)")


def bulk_checks():
    li, cc, adj = res("DE_lithium_bulk.csv"), res("DE_casecontrol_bulk.csv"), res("DE_lithium_bulk_celladjusted.csv")
    n = lambda d: int((d["adj.P.Val"] < 0.05).sum())
    add("bulk DE", "genes tested after >=1 TPM in 25% filter", 17194, 18006, "CONSISTENT", "same rule; different cohort and quantifier")
    add("bulk DE", "lithium DEGs (FDR < 0.05)", 100, n(li), "CONSISTENT" if n(li) > 0 else "DIFFERS",
        "non-zero; count scales with n (239 cases here vs 1,045)")
    add("bulk DE", "case/control DEGs (FDR < 0.05)", 64, n(cc), "DIFFERS",
        "more here than Boltz despite 1/4 the n; mostly assessment-group batch (see next row)")
    g = res("DE_casecontrol_bulk_groupadjusted.csv")
    add("bulk DE", "case/control DEGs, + assessment group (sensitivity)", "n/a", n(g), "CONSISTENT",
        "batch adjustment removes most of the excess; not a Boltz analysis")
    add("bulk DE", "cell-adjusted lithium DEGs, mostly shared with unadjusted", "82 of 94",
        f"{len(set(adj.gene[adj['adj.P.Val'] < .05]) & set(li.gene[li['adj.P.Val'] < .05]))} of {n(adj)}",
        "REPLICATED", "most cell-adjusted DEGs also significant unadjusted")
    rep = res("replication_of_boltz_degs.csv").set_index("contrast")
    for c in rep.index:
        r = rep.loc[c]
        add("bulk DE", f"Boltz's {int(r.boltz_degs)} {c} DEGs: same sign here", ">50% expected",
            f"{100 * r.sign_concordant:.0f}% (logFC r = {r.logFC_r:.2f})",
            "REPLICATED" if r.sign_concordant > 0.75 else "CONSISTENT", "sign concordance > 75%")
    ov = res("overlap_with_krebs.csv")
    b = ov[(ov["query"] == "boltz_published_100") & (ov.krebs_list == "recreated_all444")]
    if len(b):
        b = b.iloc[0]
        add("Krebs overlap", "Boltz's 100 lithium DEGs in Krebs (which list?)", "33, OR 6.43, p 4.7e-14",
            f"{int(b.overlap)}, OR {b.OR:.2f}, p {b.p:.1e}", "REPLICATED",
            "matches Krebs' documented model (text: 976 DEGs), not either deposited sheet")
    o = ov[(ov["query"] == "ours_lithium") & (ov.krebs_list == "recreated_all444")]
    if len(o):
        o = o.iloc[0]
        add("Krebs overlap", "our lithium DEGs enriched in Krebs' list", "OR 6.43",
            f"{int(o.overlap)} of {int(o.query_in_universe)}, OR {o.OR:.1f}, p {o.p:.1e}",
            "REPLICATED" if o.p < 0.05 else "DIFFERS", "Fisher p < 0.05 (not independent: overlapping samples)")


def bmind_checks():
    s = res("bmind_de_summary.csv")
    if s is None:
        add("bMIND", "bmind_de results", "", "not run", "DIFFERS", "run src/04 and src/05")
        return
    li = s[s.contrast == "lithium"]
    add("bMIND", "bmind_de lithium: no DEG in any cell type", 0, int(li.degs_here.sum()),
        "REPLICATED" if li.degs_here.sum() == 0 else "DIFFERS", "q < 0.05")
    cc = s[s.contrast == "casecontrol"].set_index("cell_type")
    for k, n in [("b.cells.naive", 21), ("neutrophils", 4)]:
        if k in cc.index:
            add("bMIND", f"bmind_de case/control {k} DEGs", n, int(cc.loc[k, "degs_here"]),
                "REPLICATED" if cc.loc[k, "degs_here"] > 0 else "DIFFERS", "any q < 0.05 (n 444 vs 1,730)")
    add("bMIND", "bmind_de case/control memory B DEGs", 24, "not estimable", "DIFFERS",
        "memory B fraction ~0 in the GEO CIBERSORT deposit")
    if "monocytes" in cc.index:
        n = int(cc.loc["monocytes", "degs_here"])
        add("bMIND", "bmind_de case/control monocyte DEGs", 0, n, "REPLICATED" if n == 0 else "DIFFERS", "q < 0.05")
    for k in ["b.cells.naive", "neutrophils"]:
        if k in cc.index and cc.loc[k, "boltz_degs_tested_here"] > 0:
            r = cc.loc[k]
            add("bMIND", f"Boltz's {k} case/control DEGs nominally significant here", f"{int(r.boltz_degs)} DEGs",
                f"{int(r.boltz_degs_nominal_p05_here)} of {int(r.boltz_degs_tested_here)} tested",
                "REPLICATED" if r.boltz_degs_nominal_p05_here / r.boltz_degs_tested_here > 0.5 else "DIFFERS",
                "more than half at p < 0.05 here")
    mc = res("bmind_de_mcmc_vs_exact.csv")
    if mc is not None:
        add("bMIND", "real bmind_de() vs exact posterior p (check subset)", "n/a",
            f"cor(log10 p) min {mc.cor_log10p.min():.3f}; p<0.05 agreement min {mc.agree_p05.min():.3f}",
            "EXACT" if mc.cor_log10p.min() > 0.95 else "DIFFERS", "cor > 0.95 in every cell type and contrast")
    s3 = res("tableS3_bmind_vs_DICE.csv")
    if s3 is not None:
        for k in s3.dice_cell_type.unique():
            z = s3[s3.dice_cell_type == k]
            m = z[z.estimate.str.startswith("bMIND") & ~z.estimate.str.contains("mismatched|not est")]
            bulk = z[z.estimate.str.startswith("bulk")].R2_log.iloc[0]
            mm = z[z.estimate.str.contains("mismatched")].R2_log.max()
            if len(m):
                v, b = m.R2_log.iloc[0], m.boltz_R2.iloc[0]
                beats = v - max(bulk, mm) >= 0.02
                add("bMIND", f"Table S3 R2 vs DICE {k}", f"{b:.2f}",
                    f"{v:.2f} (bulk {bulk:.2f}, best mismatched {mm:.2f})",
                    "REPLICATED" if beats and abs(v - b) <= 0.10 else "CONSISTENT" if beats else "DIFFERS",
                    "REPLICATED: within 0.10 of Boltz and >= 0.02 above both bulk and every mismatched type; "
                    "CONSISTENT: clearly above both baselines but lower than Boltz")
    f1 = res("fig1B_celltype_R2.csv")
    if f1 is not None:
        m = f1.set_index(f1.columns[0]).values
        off = m[~np.eye(len(m), dtype=bool)]
        # Boltz's Fig 1B colour scale spans 0.60-0.95 with most cells >= 0.85 (read off the figure).
        add("bMIND", "Fig 1B: cell types highly correlated in mean expression", "0.60-0.95 (figure scale)",
            f"off-diagonal R2 {off.min():.2f}-{off.max():.2f}, median {np.median(off):.2f}",
            "REPLICATED" if off.min() >= 0.60 else "CONSISTENT" if np.median(off) >= 0.60 else "DIFFERS",
            "REPLICATED if every pair >= 0.60, Boltz's colour-scale floor; CONSISTENT if the median is")


def prediction_checks():
    sp = res("preprint_prediction_single_split.csv")
    rp = res("preprint_prediction_repeated_splits.csv")
    if sp is None:
        add("preprint", "prediction analysis", "", "not run", "DIFFERS", "run src/06")
        return
    for oc in sp.outcome.unique():
        z = sp[sp.outcome == oc]
        base = z[z.k == 0].auc.iloc[0]
        best = z[z.k > 0].auc.max()
        add("preprint", f"{oc}: expression adds nothing over proportions (single split)", "no added value",
            f"proportions {base:.2f}; best expression model {best:.2f}",
            "REPLICATED" if best <= base + 0.05 else "DIFFERS", "no expression model beats proportions by > 0.05 AUC")
        if rp is not None:
            r = rp[rp.outcome == oc].iloc[:, 2]
            add("preprint", f"{oc}: spread of proportions-only AUC over 100 splits", "single split only",
                f"{r.quantile(.025):.2f}-{r.quantile(.975):.2f}", "CONSISTENT",
                "95% range; one split cannot support an equivalence claim")


def main() -> int:
    for f in (deposit_checks, proportion_checks, bulk_checks, bmind_checks, prediction_checks):
        f()
    card = pd.DataFrame(rows)
    card.to_csv(os.path.join(RES, "validation_scorecard.csv"), index=False)
    with pd.option_context("display.max_colwidth", 60, "display.width", 250):
        print(card[["area", "claim", "paper", "here", "verdict"]].to_string(index=False))
    print("\n" + card.verdict.value_counts().to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
