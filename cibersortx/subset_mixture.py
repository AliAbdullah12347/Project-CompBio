#!/usr/bin/env python3
"""Write a CIBERSORTx mixture file for a chosen subset of samples.

Why this exists: whether you can drop samples *after* running CIBERSORTx depends
entirely on one setting.

  batch correction OFF  ->  each sample is deconvolved independently against the
                            fixed signature matrix, so dropping samples later
                            leaves every remaining sample's fractions untouched.
  batch correction ON   ->  the correction is estimated across the whole mixture
                            matrix. Drop a sample and the correction shifts, so
                            every other sample's fractions change. You must
                            re-run on the exact set you intend to analyse.

The same applies to quantile normalisation, which is cross-sample -- we disable it.

So: for any batch-corrected run, fix the sample set first and generate its own
mixture file here.

Usage:
  python subset_mixture.py analysis      # 444 deposit-QC, minus our depth/RIN screen
  python subset_mixture.py qc444         # the deposited QC flag only
  python subset_mixture.py bp1           # bipolar I only (the lithium contrast group)
"""
from __future__ import annotations

import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
MIX = os.path.join(DATA, "mixture_TPM_hgnc_480.txt")
SCREEN = os.path.join(DATA, "qc_screen.csv")


def main() -> int:
    which = sys.argv[1] if len(sys.argv) > 1 else "analysis"
    screen = pd.read_csv(SCREEN, index_col=0)
    qc = screen["qc_pass_deposited"].astype(str).str.upper() == "TRUE"

    sets = {
        "qc444":    (qc, "deposited QC flag only"),
        "analysis": (qc & ~screen["fails_either"],
                     "deposited QC flag AND our depth/RIN screen"),
        "bp1":      (qc & ~screen["fails_either"] & (screen["diagnosis"] == "BP1"),
                     "bipolar I only, after both screens"),
    }
    if which not in sets:
        print(f"unknown set '{which}'; choose from {sorted(sets)}", file=sys.stderr)
        return 2

    mask, desc = sets[which]
    keep = [s for s in screen.index[mask]]
    print(f"Sample set '{which}': {len(keep)} samples ({desc})")

    print("Reading the 480-sample mixture ...")
    mix = pd.read_csv(MIX, sep="\t", index_col=0)
    missing = [s for s in keep if s not in mix.columns]
    if missing:
        print(f"  ERROR: {len(missing)} samples not in mixture file", file=sys.stderr)
        return 1

    out = mix[keep]
    path = os.path.join(DATA, f"mixture_TPM_hgnc_{which}_{len(keep)}.txt")
    out.index.name = "GeneSymbol"
    out.to_csv(path, sep="\t", float_format="%.4g")
    print(f"Wrote {os.path.basename(path)}  "
          f"({os.path.getsize(path)/1e6:.0f} MB, {out.shape[0]:,} x {out.shape[1]})")

    # Composition of the retained set, so the run is documented alongside its file.
    sub = screen.loc[keep]
    print("\nComposition:")
    for dx, n in sub["diagnosis"].value_counts().items():
        print(f"  {dx:10s} {n}")
    li = sub[sub["diagnosis"] != "Control"]["lithium"]
    print(f"  lithium users among cases: {int(li.sum())} / {len(li)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
