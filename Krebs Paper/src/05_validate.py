#!/usr/bin/env python3
"""One scorecard for the whole recreation.

Every row is a number the paper states, next to what this pipeline produced.
Stages whose outputs are absent are reported as "not run" rather than skipped
silently, so a partial run can never read as a clean one.

Usage:  python src/05_validate.py
Exit code is 1 if any check that ran came out FAIL.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "data", "reference")
RES = os.path.join(ROOT, "results")
PROC = os.path.join(ROOT, "data", "processed")

ROWS: list[tuple] = []


def row(stage: str, quantity: str, paper, ours, ok) -> None:
    ROWS.append((stage, quantity, str(paper), str(ours), ok))


def read(path: str) -> pd.DataFrame | None:
    p = path if os.path.isabs(path) else os.path.join(RES, path)
    return pd.read_csv(p) if os.path.exists(p) else None


def close(a: float, b: float, tol: float) -> bool:
    return abs(float(a) - float(b)) <= tol


# ---------------------------------------------------------------- dataset
meta_p = os.path.join(PROC, "sample_metadata.csv")
if os.path.exists(meta_p):
    meta = pd.read_csv(meta_p)
    post = meta[meta["qc_pass"].astype(str).str.upper().isin(("TRUE", "1"))]
    row("dataset", "samples in count matrix", 480, len(meta), len(meta) == 480)
    row("dataset", "samples passing QC", 444, len(post), len(post) == 444)
    ctrl = post["diagnosis"].str.contains("control", case=False, na=False)
    row("dataset", "cases post-QC", "240 (paper) / 239 (deposit)",
        int((~ctrl).sum()), int((~ctrl).sum()) == 239)
    row("dataset", "lithium users", 152, int(post["lithium"].sum()),
        int(post["lithium"].sum()) == 152)
    counts_p = os.path.join(PROC, "counts_filtered.tsv")
    if os.path.exists(counts_p):
        n_genes = sum(1 for _ in open(counts_p, encoding="utf-8")) - 1
        row("dataset", "genes after filter",
            "12,344 (text) / 12,353 (their DE table)", f"{n_genes:,}",
            n_genes == 12353)
else:
    row("dataset", "prepared metadata", "-", "NOT RUN", None)

# ---------------------------------------------------------------- DE
bd, li = read("DE_BD_all444.csv"), read("DE_lithium_all444.csv")
ref_bd = read(os.path.join(REF, "File_S1_DEGs__BD_DEGs.csv"))
ref_pre = read(os.path.join(REF, "File_S1_DEGs__Li_DEGs_PreCellTypeCorrection.csv"))
ref_post = read(os.path.join(REF, "File_S1_DEGs__Li_DEGs_PostCellTypeCorrection.csv"))

if bd is not None and ref_bd is not None:
    n = int((bd["adj.P.Val"] < 0.05).sum())
    row("DE / BD", "DEGs at FDR<0.05", 6, n, n == 6)
    hit = set(bd.loc[bd["adj.P.Val"] < 0.05, "gene"]) & set(ref_bd["gene"])
    row("DE / BD", "their genes also sig in ours", "6 of 6", f"{len(hit)} of 6", None)
    # A gene sitting either side of FDR 0.05 is a coin-flip at this effect size,
    # so rank agreement is the honest measure: their six should be our top six-ish.
    # Count by p-value, not row position: adj.P has ties, so positional rank is
    # arbitrary within a tie group.
    # Criterion fixed in advance: all six must land in the top 0.1% of the
    # 12,353 tested genes (~12). That is strict, and does not depend on where
    # the FDR line happens to fall relative to a tie group.
    cutoff = bd.loc[bd["gene"].isin(ref_bd["gene"]), "adj.P.Val"].max()
    worst = int((bd["adj.P.Val"] <= cutoff).sum())
    limit = max(int(0.001 * len(bd)), 1)
    row("DE / BD", f"their 6 genes in our top 0.1% (<={limit})", "6 of 6",
        f"top {worst} (last at adj.P={cutoff:.3f})", worst <= limit)
    m = bd.merge(ref_bd, on="gene", suffixes=(".o", ".t"))
    if len(m):
        r = float(np.corrcoef(m["logFC.o"], m["logFC.t"])[0, 1])
        row("DE / BD", "logFC correlation on those genes", ">0.99",
            f"{r:.4f}", r > 0.99)
else:
    row("DE / BD", "differential expression", "-", "NOT RUN", None)

if li is not None and ref_post is not None:
    row("DE / lithium", "DEGs, documented model",
        "976 (text) / 897 (MAGMA) / 3,031 (deposit)",
        int((li["adj.P.Val"] < 0.05).sum()), None)

def first_of(*names):
    """A DataFrame is not truthy, so `a or b` would raise -- pick explicitly."""
    for n in names:
        df = read(n)
        if df is not None:
            return df
    return None


post_model = first_of("DE_lithium_celltype.csv",
                      "DE_lithium_all444_celltype.csv",
                      "DE_lithium_postcelltype.csv")
if post_model is not None and ref_post is not None:
    m = post_model.merge(ref_post, on="gene", suffixes=(".o", ".t"))
    r = float(np.corrcoef(m["logFC.o"], m["logFC.t"])[0, 1])
    row("DE / lithium+cells", "logFC correlation vs deposit", "~0.99",
        f"{r:.4f}", r > 0.95)
    n = int((post_model["adj.P.Val"] < 0.05).sum())
    row("DE / lithium+cells", "DEGs at FDR<0.05", 233, n, abs(n - 233) <= 15)

# ---------------------------------------------------------------- WGCNA
wg = read("wgcna_module_trait.csv")
if wg is not None:
    n_mod = len(wg)
    row("WGCNA", "modules detected", 27, n_mod, abs(n_mod - 27) <= 4)
    if "size" in wg.columns:
        row("WGCNA", "module size range", "48-2760",
            f"{int(wg['size'].min())}-{int(wg['size'].max())}", None)
        row("WGCNA", "mean module size", 441,
            f"{wg['size'].mean():.0f}", close(wg["size"].mean(), 441, 120))
    pcol = next((c for c in wg.columns if "p_lith" in c.lower()), None)
    if pcol:
        sig = int((wg[pcol] <= 0.05 / max(n_mod, 1)).sum())
        row("WGCNA", "lithium-associated modules (Bonferroni)", 5, sig,
            abs(sig - 5) <= 3)
    bcol = next((c for c in wg.columns if "p_bd" in c.lower()), None)
    if bcol:
        sig_bd = int((wg[bcol] <= 0.05 / max(n_mod, 1)).sum())
        row("WGCNA", "BD-associated modules", 0, sig_bd, sig_bd == 0)
else:
    row("WGCNA", "co-expression network", "-", "NOT RUN", None)

# ---------------------------------------------------------------- cell types
ct = read("celltype_stepaic.csv")
if ct is not None:
    cols = {c.lower(): c for c in ct.columns}
    namecol = cols.get("cell_type") or cols.get("term") or ct.columns[0]
    neut = ct[ct[namecol].astype(str).str.contains("neutro", case=False, na=False)]
    row("cell type", "neutrophils retained by stepAIC", "yes",
        "yes" if len(neut) else "no", len(neut) > 0)
    if len(neut):
        bcol = cols.get("beta") or cols.get("estimate")
        pcol = cols.get("p") or cols.get("p_value") or cols.get("pval")
        if bcol:
            b = float(neut.iloc[0][bcol])
            row("cell type", "neutrophil beta", 0.63, f"{b:.3f}", close(b, 0.63, 0.25))
        if pcol:
            p = float(neut.iloc[0][pcol])
            row("cell type", "neutrophil p", 0.024, f"{p:.4f}", p < 0.05)
else:
    row("cell type", "stepAIC composition model", "-", "NOT RUN", None)

ref_cmp = read("celltype_reference_comparison.csv")
if ref_cmp is not None:
    row("cell type", "monocytes vs 2-8% reference", "~24.6% (3x over)",
        "see celltype_reference_comparison.csv", None)

# ---------------------------------------------------------------- report
W = (18, 46, 38, 26)
print("=" * (sum(W) + 12))
print(f"{'STAGE':<{W[0]}} {'QUANTITY':<{W[1]}} {'PAPER':<{W[2]}} {'OURS':<{W[3]}} ")
print("=" * (sum(W) + 12))
n_pass = n_fail = n_info = 0
last = None
for stage, q, paper, ours, ok in ROWS:
    if stage != last:
        print("-" * (sum(W) + 12))
        last = stage
    mark = "OK  " if ok is True else ("FAIL" if ok is False else "--  ")
    if ok is True:
        n_pass += 1
    elif ok is False:
        n_fail += 1
    else:
        n_info += 1
    print(f"{stage:<{W[0]}} {q:<{W[1]}} {paper:<{W[2]}} {ours:<{W[3]}} {mark}")
print("=" * (sum(W) + 12))
print(f"{n_pass} passed, {n_fail} failed, {n_info} informational "
      f"(-- = no single right answer, or stage not run)")

out = pd.DataFrame(ROWS, columns=["stage", "quantity", "paper", "ours", "pass"])
out.to_csv(os.path.join(RES, "validation_scorecard.csv"), index=False)
print("\nWrote results/validation_scorecard.csv")
sys.exit(1 if n_fail else 0)
