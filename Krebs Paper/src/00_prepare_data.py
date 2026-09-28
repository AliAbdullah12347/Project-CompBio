#!/usr/bin/env python3
"""Parse GSE124326 into analysis-ready tables for the R pipeline.

Kept in Python because the series matrix is a ragged key:value format that is
fiddly to parse safely in R; everything statistical downstream is R, matching
the paper. Outputs are written to data/processed/ as plain CSV/TSV so the R
scripts have no parsing logic of their own.

Krebs et al. 2020, Psychological Medicine 50(15), 2575-2586.
"""
from __future__ import annotations

import gzip
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
OUT = os.path.join(ROOT, "data", "processed")
os.makedirs(OUT, exist_ok=True)

SERIES = os.path.join(RAW, "GSE124326_series_matrix.txt.gz")
COUNTS = os.path.join(RAW, "GSE124326_count_matrix.txt.gz")

# The 22 LM22 cell types as deposited, in CIBERSORT's canonical order.
LM22 = [
    "b.cells.naive", "b.cells.memory", "plasma.cells", "t.cells.cd8",
    "t.cells.cd4.naive", "t.cells.cd4.memory.resting",
    "t.cells.cd4.memory.activated", "t.cells.follicular.helper",
    "t.cells.regulatory", "t.cells.gamma.delta", "nk.cells.resting",
    "nk.cells.activated", "monocytes", "macrophages.m0", "macrophages.m1",
    "macrophages.m2", "dendritic.cells.resting", "dendritic.cells.activated",
    "mast.cells.resting", "mast.cells.activated", "eosinophils", "neutrophils",
]


def parse_series(path: str) -> pd.DataFrame:
    """Read the series matrix into one row per sample.

    Characteristics lines are ragged: a given field does not occupy the same
    line index for every sample (Sex is present for 478 of 480). So every
    'field: value' cell is assigned by its own key, never by row position.
    """
    simple, characteristics = {}, []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key, vals = parts[0], [p.strip('"') for p in parts[1:]]
            if key == "!Sample_characteristics_ch1":
                characteristics.append(vals)
            elif key in ("!Sample_title", "!Sample_geo_accession"):
                simple[key] = vals

    n = len(simple["!Sample_title"])
    meta = pd.DataFrame({
        "title": simple["!Sample_title"],
        "gsm": simple["!Sample_geo_accession"],
    })

    fields: dict[str, list] = {}
    for row in characteristics:
        for i, cell in enumerate(row):
            if ":" not in cell:
                continue
            field, value = cell.split(":", 1)
            fields.setdefault(field.strip(), [None] * n)[i] = value.strip()

    for field, values in fields.items():
        meta[field] = values
    return meta


def tidy(meta: pd.DataFrame) -> pd.DataFrame:
    """Rename to analysis names and coerce types."""
    df = pd.DataFrame({
        "title": meta["title"],
        "gsm": meta["gsm"],
        "diagnosis": meta["bipolar disorder diagnosis"],
        "age": pd.to_numeric(meta["age"], errors="coerce"),
        "sex": meta["Sex"],
        "lithium": pd.to_numeric(
            meta["lithium use (non-user=0, user = 1)"], errors="coerce"),
        "tobacco": meta["tobacco use"],
        "group": meta["assessment group"],
        "rin": pd.to_numeric(meta["rin"], errors="coerce"),
        "plate": meta["sequencing plate"],
        "seqpc1": pd.to_numeric(meta["sequencing metric pc1"], errors="coerce"),
        "seqpc2": pd.to_numeric(meta["sequencing metric pc2"], errors="coerce"),
        "seqpc3": pd.to_numeric(meta["sequencing metric pc3"], errors="coerce"),
        "qc_pass": meta["included in final analysis"].str.strip().str.upper()
                   == "TRUE",
    })
    for ct in LM22:
        df[ct] = pd.to_numeric(meta[ct], errors="coerce")
    return df


def load_counts(path: str) -> pd.DataFrame:
    """Read the count matrix; columns are sample titles with a .counts suffix."""
    counts = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    counts.columns = [c[:-7] if c.endswith(".counts") else c
                      for c in counts.columns]
    return counts


def check(label: str, got, want, results: list) -> None:
    ok = got == want
    results.append(ok)
    flag = "OK  " if ok else "FAIL"
    detail = "" if ok else f"   <-- expected {want}"
    print(f"  [{flag}] {label}: {got}{detail}")


def main() -> int:
    print("Parsing series matrix ...")
    meta = tidy(parse_series(SERIES))
    print("Reading count matrix (16 MB gz) ...")
    counts = load_counts(COUNTS)

    results: list[bool] = []
    print("\n--- Dataset facts ---")
    check("count matrix genes", counts.shape[0], 57820, results)
    check("count matrix samples", counts.shape[1], 480, results)
    check("metadata samples", len(meta), 480, results)

    # Ensembl IDs carry version suffixes; they are unique once stripped.
    stripped = counts.index.str.replace(r"\.\d+$", "", regex=True)
    check("duplicate IDs after stripping version", int(stripped.duplicated().sum()), 0, results)

    # ENSGR rows are pseudoautosomal-region Y duplicates of ENSG00000... genes.
    ensgr = counts.index.str.startswith("ENSGR")
    check("ENSGR rows", int(ensgr.sum()), 47, results)
    check("ENSGR total counts", int(counts.loc[ensgr].values.sum()), 1518, results)

    # Join on the title key, never on column order.
    check("titles matching metadata", int(counts.columns.isin(meta["title"]).sum()), 480, results)
    check("column order == metadata order",
          list(counts.columns) == list(meta["title"]), True, results)
    check("titles starting with a digit",
          int(meta["title"].str[0].str.isdigit().sum()), 480, results)

    counts = counts.loc[~ensgr]
    check("genes after dropping ENSGR", counts.shape[0], 57773, results)

    print("\n--- Sample QC ---")
    check("QC-passing samples", int(meta["qc_pass"].sum()), 444, results)
    post = meta[meta["qc_pass"]]
    dx = post["diagnosis"].value_counts()
    pre = meta["diagnosis"].value_counts()
    print(f"    pre-QC diagnosis: {dict(pre)}")
    print(f"    post-QC diagnosis: {dict(dx)}")

    bp1 = [k for k in dx.index if "1" in k or "I" == k.strip()[-1:]]
    print(f"    (BP1 label resolved as: {bp1})")

    print("\n--- Krebs Table 1 (paper claims 240 cases / 204 controls) ---")
    # The deposit disagrees with the publication by exactly one subject, in both
    # the case/control split and the non-lithium case count. Reported, not
    # "fixed": the deposit is what we can actually analyse.
    ctrl = post["diagnosis"].str.contains("control", case=False, na=False)
    n_ctrl, n_case = int(ctrl.sum()), int((~ctrl).sum())
    n_li = int(post.loc[~ctrl, "lithium"].sum())
    n_noli = int((post.loc[~ctrl, "lithium"] == 0).sum())
    for label, got, paper in [
        ("controls post-QC", n_ctrl, 204),
        ("cases post-QC", n_case, 240),
        ("non-lithium cases", n_noli, 88),
    ]:
        note = "" if got == paper else f"   [deposit {got} vs paper {paper}: known 1-subject discrepancy]"
        print(f"  [note] {label}: {got}{note}")
    check("lithium users (cases)", n_li, 152, results)
    check("lithium users among controls", int(post.loc[ctrl, "lithium"].sum()), 0, results)

    print("\n--- Covariate completeness post-QC ---")
    cov = ["age", "sex", "tobacco", "group", "rin", "plate",
           "seqpc1", "seqpc2", "seqpc3"]
    miss = {c: int(post[c].isna().sum()) for c in cov}
    check("covariates with missing values post-QC",
          sum(1 for v in miss.values() if v), 0, results)
    print(f"    sex: {dict(post['sex'].value_counts())}")
    print(f"    tobacco: {dict(post['tobacco'].value_counts())}")
    print(f"    group: {dict(post['group'].value_counts())}")

    print("\n--- Gene filter (>10 counts in >=90% of 444) ---")
    post_counts = counts[post["title"].tolist()]
    keep = (post_counts > 10).sum(axis=1) >= 0.90 * post_counts.shape[1]
    n_keep = int(keep.sum())
    # The paper's text says 12,344, but the authors' own deposited DE table has
    # 12,353 rows -- which is exactly what this filter produces. 12,344 is the
    # WGCNA gene count; 9 genes are dropped between the DE and network stages.
    check("genes passing filter (matches their DE table)", n_keep, 12353, results)
    print("    note: paper's text says 12,344 -- that is the WGCNA count, not the DE count")

    print("\n--- Deposited CIBERSORT fractions ---")
    frac = meta[LM22]
    # Deposited to ~15 significant figures but rounded, so the row sums miss 1 by
    # up to ~6e-5. That is deposit rounding, not a renormalisation error; a 1e-6
    # tolerance would wrongly flag 106 of 480 samples.
    dev = float(np.abs(frac.sum(axis=1) - 1).max())
    check("fractions sum to 1 (max abs deviation < 1e-4)", bool(dev < 1e-4), True, results)
    print(f"    max |rowsum - 1| = {dev:.2e}")
    zeros = [c for c in LM22 if (frac[c] == 0).all()]
    check("structural zero cell types", len(zeros), 4, results)
    print(f"    structural zeros: {zeros}")

    print("\n--- Writing processed tables ---")
    meta.to_csv(os.path.join(OUT, "sample_metadata.csv"), index=False)
    post_counts.loc[keep].to_csv(
        os.path.join(OUT, "counts_filtered.tsv"), sep="\t")
    counts[post["title"].tolist()].to_csv(
        os.path.join(OUT, "counts_qc_allgenes.tsv.gz"), sep="\t",
        compression="gzip")
    for name, shape in [
        ("sample_metadata.csv", meta.shape),
        ("counts_filtered.tsv", (n_keep, post_counts.shape[1])),
        ("counts_qc_allgenes.tsv.gz", (counts.shape[0], post_counts.shape[1])),
    ]:
        print(f"    {name}: {shape[0]:,} x {shape[1]:,}")

    failed = results.count(False)
    print(f"\n{len(results) - failed}/{len(results)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
