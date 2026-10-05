#!/usr/bin/env python3
"""Build a full-transcriptome TPM mixture for the 474-sample cohort.

Unlike the filtered mixture that `build_cohort.py` derives, this keeps every
gene in the annotation rather than only protein-coding plus IG/TR segments.

Why you might want it. CIBERSORT's support vector regression uses only the
intersection with LM22's 547 signature genes, so for a plain run the extra rows
change nothing. But CIBERSORTx's **B-mode batch correction operates on the whole
mixture matrix** against the signature's source profiles, so a larger gene set
gives that correction more to work with. Krebs et al. also fed CIBERSORT their
unfiltered TPM matrix, so this is the closer reproduction of their run.

Why you might not. The file is several times larger and may exceed the web
portal's upload limit, in which case use the filtered mixture or the Docker
version of CIBERSORTx.

TPM needs a per-gene length. We compute the union of each gene's exons from the
GENCODE v19 GTF -- the same annotation the counts were generated against, and
the correct denominator for HTSeq union mode, which assigns a read to a gene
only if it overlaps that gene's exons and no other's.

    python data/make_mixture_full.py
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
COHORT = HERE / "cohort_474"
GTF = HERE.parent / "cibersortx" / "data" / "gencode.v19.annotation.gtf.gz"

GID = re.compile(r'gene_id "([^"]+)"')
GNAME = re.compile(r'gene_name "([^"]+)"')


def union_exon_lengths() -> pd.DataFrame:
    """gene_id -> (symbol, union exon length in bases).

    Exons of different transcripts overlap, so their lengths cannot simply be
    summed; the intervals are merged first. Summing them instead inflates the
    length of multi-transcript genes and deflates their TPM.
    """
    if not GTF.exists():
        raise SystemExit(f"GENCODE v19 GTF not found at {GTF}")
    print(f"  parsing {GTF.name} ...", flush=True)
    exons: dict[str, list[tuple[int, int]]] = {}
    names: dict[str, str] = {}
    with gzip.open(GTF, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[2] != "exon":
                continue
            attrs = f[8]
            gid = GID.search(attrs).group(1).split(".")[0]
            exons.setdefault(gid, []).append((int(f[3]), int(f[4])))
            if gid not in names:
                m = GNAME.search(attrs)
                names[gid] = m.group(1) if m else gid

    out = {}
    for gid, spans in exons.items():
        spans.sort()
        total, cur_s, cur_e = 0, spans[0][0], spans[0][1]
        for s, e in spans[1:]:
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                total += cur_e - cur_s + 1
                cur_s, cur_e = s, e
        out[gid] = total + cur_e - cur_s + 1
    print(f"  {len(out):,} genes with exon models")
    return pd.DataFrame({"symbol": pd.Series(names), "length": pd.Series(out)})


def main() -> int:
    print("Reading the 474-sample count matrix ...", flush=True)
    counts = pd.read_csv(COHORT / "counts_474.tsv.gz", sep="\t", index_col=0)
    print(f"  {counts.shape[0]:,} genes x {counts.shape[1]} samples")

    # Pseudoautosomal-region Y duplicates double-count 1,518 reads across 47
    # rows. They must go before TPM, or they steal a share of the 1e6.
    counts = counts.loc[~counts.index.str.startswith("ENSGR")]
    counts.index = counts.index.str.replace(r"\.\d+$", "", regex=True)
    print(f"  {counts.shape[0]:,} after dropping ENSGR rows")

    ann = union_exon_lengths()
    shared = counts.index.intersection(ann.index)
    missing = len(counts) - len(shared)
    if missing:
        print(f"  {missing} genes absent from the GTF, dropped")
    counts = counts.loc[shared]
    ann = ann.loc[shared]

    print("Converting to TPM ...", flush=True)
    rate = counts.to_numpy(dtype=np.float64) / (ann["length"].to_numpy()[:, None] / 1000.0)
    tpm = pd.DataFrame(rate / rate.sum(axis=0) * 1e6,
                       index=ann["symbol"].to_numpy(), columns=counts.columns)

    dups = int(tpm.index.duplicated().sum())
    # Several Ensembl genes share one HGNC symbol. TPM is additive within a
    # gene, so the duplicates are summed rather than averaged or dropped.
    tpm = tpm.groupby(level=0, sort=True).sum()
    print(f"  {len(tpm):,} unique symbols ({dups:,} duplicate rows summed)")

    col = tpm.sum(axis=0)
    print(f"  column sums: {col.min():,.1f} to {col.max():,.1f}  (target 1,000,000)")

    tpm.index.name = "GeneSymbol"
    path = COHORT / f"mixture_TPM_hgnc_full_{tpm.shape[1]}.txt"
    tpm.to_csv(path, sep="\t", float_format="%.4g")
    size = path.stat().st_size / 1e6
    print(f"\nWrote {path.name}  ({size:.0f} MB, {tpm.shape[0]:,} x {tpm.shape[1]})")
    if size > 100:
        print("  NOTE: this may exceed the CIBERSORTx web upload limit. If the")
        print("  portal rejects it, use mixture_TPM_hgnc_474.txt or the Docker")
        print("  version of CIBERSORTx.")

    zero = int((tpm.sum(axis=1) == 0).sum())
    print(f"  genes zero in all {tpm.shape[1]} samples: {zero:,} (kept, as Krebs did)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
