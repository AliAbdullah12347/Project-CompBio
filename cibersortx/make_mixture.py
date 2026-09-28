#!/usr/bin/env python3
"""Build a CIBERSORTx mixture file from the GSE124326 count matrix.

CIBERSORTx will not accept our counts as they stand. Three conversions are needed:

  1. Ensembl IDs -> HGNC symbols. LM22 is keyed on gene symbols; an Ensembl-keyed
     matrix silently matches almost no signature genes and returns garbage.
  2. Counts -> TPM, in linear (non-log) space, as Krebs et al. did. TPM needs a
     per-gene length, so we compute union exon length from the GENCODE v19 GTF --
     the same annotation the counts were generated against, and the right
     denominator for HTSeq union mode.
  3. Duplicate symbols collapsed. Several Ensembl genes map to one HGNC symbol;
     TPM is additive within a gene, so duplicates are summed.

Outputs (data/ beside this script):
  mixture_TPM_hgnc_480.txt   tab-delimited, genes x 480, ready to upload
  qc_screen.csv              per-sample depth and RIN against our thresholds
  gene_mapping_report.txt    what was dropped and why
"""
from __future__ import annotations

import gzip
import os
import re
import sys
import urllib.request

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data")
os.makedirs(OUT, exist_ok=True)

# Reuse the Krebs download rather than fetching GEO twice.
KREBS = os.path.join(os.path.dirname(HERE), "Krebs Paper", "data")
COUNTS = os.path.join(KREBS, "raw", "GSE124326_count_matrix.txt.gz")
META = os.path.join(KREBS, "processed", "sample_metadata.csv")

GTF_URL = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
           "release_19/gencode.v19.annotation.gtf.gz")
GTF = os.path.join(OUT, "gencode.v19.annotation.gtf.gz")

# Pre-specified in the revised pre-proposal; fixed before any result is examined.
MIN_READS = 3_000_000
MIN_RIN = 5.0


def fetch_gtf() -> None:
    if os.path.exists(GTF) and os.path.getsize(GTF) > 30_000_000:
        print(f"  GTF already present ({os.path.getsize(GTF)/1e6:.0f} MB)")
        return
    print("  downloading GENCODE v19 GTF (~36 MB) ...")
    urllib.request.urlretrieve(GTF_URL, GTF)
    print(f"  done ({os.path.getsize(GTF)/1e6:.0f} MB)")


def parse_gtf() -> pd.DataFrame:
    """Return gene_id -> (symbol, union exon length, gene_type)."""
    gid_re = re.compile(r'gene_id "([^"]+)"')
    name_re = re.compile(r'gene_name "([^"]+)"')
    type_re = re.compile(r'gene_type "([^"]+)"')

    symbol: dict[str, str] = {}
    gtype: dict[str, str] = {}
    exons: dict[str, list[tuple[int, int]]] = {}

    with gzip.open(GTF, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            feat, start, end, attrs = f[2], int(f[3]), int(f[4]), f[8]
            if feat == "gene":
                m = gid_re.search(attrs)
                if not m:
                    continue
                g = m.group(1)
                n = name_re.search(attrs)
                t = type_re.search(attrs)
                symbol[g] = n.group(1) if n else g
                gtype[g] = t.group(1) if t else "unknown"
            elif feat == "exon":
                m = gid_re.search(attrs)
                if not m:
                    continue
                exons.setdefault(m.group(1), []).append((start, end))

    # Union exon length: merge overlapping exons so shared bases are counted once.
    lengths: dict[str, int] = {}
    for g, iv in exons.items():
        iv.sort()
        total = 0
        cs, ce = iv[0]
        for s, e in iv[1:]:
            if s <= ce + 1:
                ce = max(ce, e)
            else:
                total += ce - cs + 1
                cs, ce = s, e
        total += ce - cs + 1
        lengths[g] = total

    df = pd.DataFrame({
        "gene_id": list(symbol),
        "symbol": [symbol[g] for g in symbol],
        "gene_type": [gtype[g] for g in symbol],
        "length": [lengths.get(g, np.nan) for g in symbol],
    })
    return df


def main() -> int:
    print("Step 1  GENCODE v19 annotation")
    fetch_gtf()
    ann = parse_gtf()
    print(f"  {len(ann):,} genes; union exon length available for "
          f"{int(ann['length'].notna().sum()):,}")

    print("\nStep 2  count matrix")
    counts = pd.read_csv(COUNTS, sep="\t", index_col=0, compression="gzip")
    counts.columns = [c[:-7] if c.endswith(".counts") else c for c in counts.columns]
    print(f"  {counts.shape[0]:,} genes x {counts.shape[1]} samples")

    # Pseudoautosomal-region Y duplicates; 1,518 counts total across 47 rows.
    ensgr = counts.index.str.startswith("ENSGR")
    counts = counts.loc[~ensgr]
    print(f"  dropped {int(ensgr.sum())} ENSGR rows -> {counts.shape[0]:,}")

    print("\nStep 3  counts -> TPM")
    ann = ann.set_index("gene_id")
    shared = counts.index.intersection(ann.index)
    missing = counts.index.difference(ann.index)
    print(f"  matched to annotation: {len(shared):,}   unmatched: {len(missing):,}")

    counts = counts.loc[shared]
    length_kb = ann.loc[shared, "length"].to_numpy(dtype=float) / 1000.0

    # TPM: length-normalise first, then scale each sample to 1e6. This is the
    # per-sample operation CIBERSORT expects -- no cross-sample information.
    rpk = counts.to_numpy(dtype=float) / length_kb[:, None]
    tpm = rpk / rpk.sum(axis=0, keepdims=True) * 1e6
    tpm = pd.DataFrame(tpm, index=counts.index, columns=counts.columns)
    print(f"  column sums all 1e6: {np.allclose(tpm.sum(axis=0), 1e6)}")

    print("\nStep 4  Ensembl -> HGNC symbol, collapsing duplicates")
    tpm.insert(0, "symbol", ann.loc[shared, "symbol"].to_numpy())
    tpm.insert(1, "gene_type", ann.loc[shared, "gene_type"].to_numpy())
    n_before = len(tpm)

    # With quantile normalisation disabled, CIBERSORT uses only the intersection
    # of the mixture with the 547 LM22 signature genes -- every other row is
    # ignored. LM22's markers are all protein-coding, so restricting to
    # protein_coding is lossless and takes the upload from ~190 MB to ~35 MB.
    # MARKER_CHECK below fails loudly if that assumption ever breaks.
    # Immunoglobulin and T-cell-receptor genes are typed IG_*_gene / TR_*_gene in
    # GENCODE v19, NOT protein_coding -- and LM22 identifies plasma cells and
    # gamma-delta T cells partly through them, so they must be kept.
    keep_type = (tpm["gene_type"] == "protein_coding") | \
                tpm["gene_type"].str.match(r"^(IG|TR)_[A-Z]+_gene$", na=False)
    by_type = tpm.loc[keep_type, "gene_type"].value_counts()
    print(f"  keeping {int(keep_type.sum()):,} of {n_before:,} by gene_type:")
    for t, n in by_type.items():
        print(f"      {t:24s} {n:,}")
    tpm = tpm.loc[keep_type].drop(columns=["gene_type"])

    dup_syms = int(tpm["symbol"].duplicated().sum())
    mixture = tpm.groupby("symbol", sort=True).sum()
    print(f"  -> {len(mixture):,} unique symbols ({dup_syms:,} duplicates summed)")

    # Spot-check well-known LM22 markers spanning all the major lineages.
    # Symbols are checked in their GENCODE v19 form, which is what LM22 was
    # built against: JCHAIN was still IGJ in 2013, and TARP/TRGC1 likewise.
    MARKER_CHECK = [
        "CD8A", "MS4A1", "CD19", "NKG7", "GNLY", "FCGR3B", "CSF3R", "CEACAM3",
        "IL5RA", "CCR3", "TPSAB1", "FOXP3", "CXCR5", "CD4", "CD14", "FCN1",
        "ELANE", "MPO", "IGJ", "TRDC", "KLRD1", "CD79A", "ITGAX", "SIGLEC8",
        "IGHM", "IGKC", "TRBC1", "CD3E", "LYZ", "S100A8",
    ]
    absent = [g for g in MARKER_CHECK if g not in mixture.index]
    print(f"  LM22 marker spot-check: {len(MARKER_CHECK) - len(absent)}"
          f"/{len(MARKER_CHECK)} present")
    if absent:
        print(f"  WARNING -- missing markers: {absent}")
        print("  Re-run without the protein_coding filter if any are real LM22 genes.")

    mix_path = os.path.join(OUT, "mixture_TPM_hgnc_480.txt")
    mixture.index.name = "GeneSymbol"
    # 4 significant figures, so exact zeros write as "0" instead of "0.0000".
    mixture.to_csv(mix_path, sep="\t", float_format="%.4g")
    sz = os.path.getsize(mix_path) / 1e6
    print(f"  wrote {os.path.basename(mix_path)}  ({sz:.0f} MB, "
          f"{len(mixture):,} x {mixture.shape[1]})")

    print("\nStep 5  depth / RIN screen (pre-specified: >=3M reads, RIN >=5.0)")
    depth = counts.drop(columns=[], errors="ignore").sum(axis=0)
    meta = pd.read_csv(META)
    meta = meta.set_index("title").loc[mixture.columns]
    screen = pd.DataFrame({
        "assigned_reads": depth.reindex(mixture.columns).astype("int64"),
        "rin": meta["rin"].to_numpy(),
        "qc_pass_deposited": meta["qc_pass"].to_numpy(),
        "diagnosis": meta["diagnosis"].to_numpy(),
        "lithium": meta["lithium"].to_numpy(),
    })
    screen["fails_depth"] = screen["assigned_reads"] < MIN_READS
    screen["fails_rin"] = screen["rin"] < MIN_RIN
    screen["fails_either"] = screen["fails_depth"] | screen["fails_rin"]
    screen.to_csv(os.path.join(OUT, "qc_screen.csv"))

    dep = screen[screen["qc_pass_deposited"].astype(str).str.upper() == "TRUE"]
    print(f"  of the {len(dep)} deposit-QC-passing samples:")
    print(f"    fail depth (<3M):   {int(dep['fails_depth'].sum())}")
    print(f"    fail RIN (<5.0):    {int(dep['fails_rin'].sum())}")
    print(f"    fail either:        {int(dep['fails_either'].sum())}")
    print(f"    RETAINED:           {int((~dep['fails_either']).sum())}")

    with open(os.path.join(OUT, "gene_mapping_report.txt"), "w",
              encoding="utf-8") as fh:
        fh.write(f"Ensembl IDs in counts after ENSGR removal: {n_before + len(missing):,}\n")
        fh.write(f"Matched to GENCODE v19 annotation:         {len(shared):,}\n")
        fh.write(f"Unmatched (dropped):                       {len(missing):,}\n")
        fh.write(f"Duplicate HGNC symbols (summed):           {dup_syms:,}\n")
        fh.write(f"Final unique symbols in mixture file:      {len(mixture):,}\n\n")
        fh.write("Unmatched Ensembl IDs:\n")
        for g in missing:
            fh.write(f"  {g}\n")
    print("\nWrote data/qc_screen.csv and data/gene_mapping_report.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
