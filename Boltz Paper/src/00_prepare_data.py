#!/usr/bin/env python3
"""Parse GSE124326 into the inputs Boltz et al. (2024) worked from: TPM, cell fractions, covariates.

Boltz et al. quantified 1,730 low-coverage samples with kallisto (GENCODE v33lift37) and worked in
TPM throughout. Their data (dbGaP phs002856) is controlled-access, so this recreation runs their
methods on the public GSE124326 cohort (Krebs et al. 2020, same Utrecht cohort, same lab). GEO
provides HTSeq gene counts, not transcript TPM, so TPM is rebuilt here from counts and the GENCODE
v19 union-exon length of each gene -- the standard gene-level approximation. It is an
approximation of kallisto TPM (no effective-length or multi-mapping correction), recorded as a
deviation in REPRODUCIBILITY.md.

Outputs (data/processed/):
  sample_metadata.csv   444 QC-passing samples, tidy covariates + 22 LM22 fractions + group flags
  gene_annotation.csv   57,773 genes: id, symbol, biotype, union-exon length
  tpm_filtered.tsv.gz   genes passing Boltz's filter (>= 1 TPM in >= 25% of samples) x 444
"""
from __future__ import annotations

import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "raw")
EXT = os.path.join(ROOT, "data", "external")
OUT = os.path.join(ROOT, "data", "processed")
os.makedirs(OUT, exist_ok=True)

SERIES = os.path.join(RAW, "GSE124326_series_matrix.txt.gz")
COUNTS = os.path.join(RAW, "GSE124326_count_matrix.txt.gz")
GTF = os.path.join(EXT, "gencode_v19", "gencode.v19.annotation.gtf.gz")

LM22 = [
    "b.cells.naive", "b.cells.memory", "plasma.cells", "t.cells.cd8",
    "t.cells.cd4.naive", "t.cells.cd4.memory.resting",
    "t.cells.cd4.memory.activated", "t.cells.follicular.helper",
    "t.cells.regulatory", "t.cells.gamma.delta", "nk.cells.resting",
    "nk.cells.activated", "monocytes", "macrophages.m0", "macrophages.m1",
    "macrophages.m2", "dendritic.cells.resting", "dendritic.cells.activated",
    "mast.cells.resting", "mast.cells.activated", "eosinophils", "neutrophils",
]

# Boltz: ">= 1 TPM in at least 436 individuals (about 25% of the total 1,730)". 436/1730 = 25.2%;
# the rule is carried over as a proportion, since the absolute 436 is specific to their cohort.
TPM_MIN = 1.0
TPM_FRAC = 0.25


def parse_series(path: str) -> pd.DataFrame:
    """One row per sample. Characteristics lines are ragged, so each 'field: value' is assigned by
    its own key, never by line position (Sex sits on different lines for different samples)."""
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
    meta = pd.DataFrame({"title": simple["!Sample_title"], "gsm": simple["!Sample_geo_accession"]})
    fields: dict[str, list] = {}
    for row in characteristics:
        for i, cell in enumerate(row):
            if ":" in cell:
                field, value = cell.split(":", 1)
                fields.setdefault(field.strip(), [None] * n)[i] = value.strip()
    for field, values in fields.items():
        meta[field] = values
    return meta


def tidy(meta: pd.DataFrame) -> pd.DataFrame:
    df = pd.DataFrame({
        "title": meta["title"],
        "gsm": meta["gsm"],
        "diagnosis": meta["bipolar disorder diagnosis"],
        "age": pd.to_numeric(meta["age"], errors="coerce"),
        "sex": meta["Sex"],
        "lithium": pd.to_numeric(meta["lithium use (non-user=0, user = 1)"], errors="coerce"),
        "rin": pd.to_numeric(meta["rin"], errors="coerce"),
        # Not a Boltz covariate; kept for a sensitivity check, because every case is in group A
        # while controls span A and B, so group is a collection batch confounded with diagnosis.
        "group": meta["assessment group"],
        "qc_pass": meta["included in final analysis"].str.strip().str.upper() == "TRUE",
    })
    for ct in LM22:
        df[ct] = pd.to_numeric(meta[ct], errors="coerce")
    return df


def gene_lengths(path: str) -> pd.DataFrame:
    """Union-exon length per gene: overlapping exons of all transcripts merged, then summed. This is
    the length HTSeq union-mode counts are distributed over, so it is the right TPM denominator."""
    exons: dict[str, list[tuple[int, int]]] = {}
    info: dict[str, tuple[str, str]] = {}
    gid_re = re.compile(r'gene_id "([^"]+)"')
    name_re = re.compile(r'gene_name "([^"]+)"')
    type_re = re.compile(r'gene_type "([^"]+)"')
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] == "gene":
                gid = gid_re.search(f[8]).group(1)
                info[gid] = (name_re.search(f[8]).group(1), type_re.search(f[8]).group(1))
            elif f[2] == "exon":
                gid = gid_re.search(f[8]).group(1)
                exons.setdefault(gid, []).append((int(f[3]), int(f[4])))
    rows = []
    for gid, iv in exons.items():
        iv.sort()
        total, (s0, e0) = 0, iv[0]
        for s, e in iv[1:]:
            if s <= e0 + 1:
                e0 = max(e0, e)
            else:
                total += e0 - s0 + 1
                s0, e0 = s, e
        total += e0 - s0 + 1
        rows.append((gid, info[gid][0], info[gid][1], total))
    return pd.DataFrame(rows, columns=["gene_id", "symbol", "biotype", "length"]).set_index("gene_id")


def check(label, got, want, results):
    ok = got == want
    results.append(ok)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {label}: {got}" + ("" if ok else f"   <-- expected {want}"))


def main() -> int:
    results: list[bool] = []
    print("Parsing series matrix ...")
    meta = tidy(parse_series(SERIES))
    print("Reading count matrix ...")
    counts = pd.read_csv(COUNTS, sep="\t", index_col=0, compression="gzip")
    counts.columns = [c[:-7] if c.endswith(".counts") else c for c in counts.columns]

    print("\n--- Dataset facts (CLAUDE.md section 2) ---")
    check("count matrix shape", counts.shape, (57820, 480), results)
    check("titles matching metadata", int(counts.columns.isin(meta["title"]).sum()), 480, results)
    ensgr = counts.index.str.startswith("ENSGR")
    check("ENSGR (PAR-Y duplicate) rows", int(ensgr.sum()), 47, results)
    counts = counts.loc[~ensgr]
    check("genes after dropping ENSGR", counts.shape[0], 57773, results)
    check("QC-passing samples", int(meta["qc_pass"].sum()), 444, results)

    post = meta[meta["qc_pass"]].reset_index(drop=True)
    # Boltz's cases were "all individuals diagnosed with BP or SCZ"; GSE124326 has no SCZ, so cases
    # are BP-I plus BP-II. Their lithium contrast used "only diagnosed individuals".
    post["case"] = (~post["diagnosis"].str.contains("control", case=False)).astype(int)
    print(f"    post-QC diagnosis: {post['diagnosis'].value_counts().to_dict()}")
    check("cases (BP-I + BP-II) post-QC", int(post["case"].sum()), 239, results)
    check("lithium users among cases", int(post.loc[post.case == 1, "lithium"].sum()), 152, results)
    check("lithium users among controls", int(post.loc[post.case == 0, "lithium"].sum()), 0, results)
    check("missing age/sex/RIN post-QC", int(post[["age", "sex", "rin"]].isna().sum().sum()), 0, results)

    print("\n--- Gene lengths from GENCODE v19 ---")
    ann = gene_lengths(GTF)
    check("GTF genes", len(ann), 57820, results)
    ann = ann.loc[counts.index]  # raises if any count-matrix gene is missing from the GTF
    check("count-matrix genes with a length", int(ann["length"].notna().sum()), 57773, results)

    print("\n--- TPM ---")
    # Per-sample operation, so computing it once on all genes before filtering is leak-free and
    # matches kallisto, whose TPM is normalised over the whole transcriptome.
    x = counts[post["title"]].to_numpy(dtype=float)
    rate = x / (ann["length"].to_numpy()[:, None] / 1e3)
    tpm = rate / rate.sum(axis=0, keepdims=True) * 1e6
    tpm = pd.DataFrame(tpm, index=counts.index, columns=post["title"])
    check("TPM columns sum to 1e6", bool(np.allclose(tpm.sum(axis=0), 1e6)), True, results)
    print(f"    max TPM = {tpm.values.max():,.0f}  (Boltz: 'largest expression value > 50', so log2 before bMIND)")

    need = int(np.ceil(TPM_FRAC * tpm.shape[1]))
    keep = (tpm >= TPM_MIN).sum(axis=1) >= need
    print(f"    filter: >= {TPM_MIN} TPM in >= {need} of {tpm.shape[1]} samples -> {int(keep.sum()):,} genes "
          f"(Boltz, n = 1,730: 17,194)")

    print("\n--- Deposited CIBERSORT fractions ---")
    dev = float(np.abs(post[LM22].sum(axis=1) - 1).max())
    check("fractions sum to 1 (tolerance 1e-4; deposit is rounded)", dev < 1e-4, True, results)

    print("\n--- Writing ---")
    post.to_csv(os.path.join(OUT, "sample_metadata.csv"), index=False)
    ann.to_csv(os.path.join(OUT, "gene_annotation.csv"))
    tpm.loc[keep].to_csv(os.path.join(OUT, "tpm_filtered.tsv.gz"), sep="\t", float_format="%.6g",
                         compression="gzip")
    print(f"    sample_metadata.csv  {post.shape[0]} x {post.shape[1]}")
    print(f"    gene_annotation.csv  {ann.shape[0]:,} genes")
    print(f"    tpm_filtered.tsv.gz  {int(keep.sum()):,} x {tpm.shape[1]}")

    failed = results.count(False)
    print(f"\n{len(results) - failed}/{len(results)} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
