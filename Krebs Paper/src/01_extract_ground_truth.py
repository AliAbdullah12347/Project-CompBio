#!/usr/bin/env python3
"""Extract the paper's published result tables into CSV for validation.

These are the authors' own outputs, so they are the reference the recreation is
checked against. xlsx is parsed straight from its OOXML rather than via a
spreadsheet library, to avoid adding a dependency for two files.

Supplementary File S1 (DEG lists)   -> sup003.xlsx
Supplementary File S2 (WGCNA modules) -> sup001.xlsx
"""
from __future__ import annotations

import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

import pandas as pd

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUP = os.path.join(ROOT, "data", "supplementary")
OUT = os.path.join(ROOT, "data", "reference")
os.makedirs(OUT, exist_ok=True)


def _shared_strings(z: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    out = []
    for si in root.findall(NS + "si"):
        # A string can be split across several <t> runs; concatenate them.
        out.append("".join(t.text or "" for t in si.iter(NS + "t")))
    return out


def _col_index(ref: str) -> int:
    """'BC12' -> zero-based column index."""
    letters = re.match(r"([A-Z]+)", ref).group(1)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_sheets(path: str) -> dict[str, pd.DataFrame]:
    z = zipfile.ZipFile(path)
    strings = _shared_strings(z)

    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels}

    sheets = {}
    for sh in wb.iter(NS + "sheet"):
        name = sh.get("name")
        tgt = target[sh.get(REL + "id")].lstrip("/")
        if not tgt.startswith("xl/"):
            tgt = "xl/" + tgt

        root = ET.fromstring(z.read(tgt))
        rows = []
        for row in root.iter(NS + "row"):
            cells: dict[int, object] = {}
            for c in row.findall(NS + "c"):
                ref, ctype = c.get("r"), c.get("t")
                v = c.find(NS + "v")
                if ctype == "s":                       # shared string
                    val = strings[int(v.text)] if v is not None else None
                elif ctype == "inlineStr":
                    is_ = c.find(NS + "is")
                    val = "".join(t.text or "" for t in is_.iter(NS + "t")) if is_ is not None else None
                elif v is None:
                    val = None
                else:
                    try:
                        val = float(v.text)
                    except (TypeError, ValueError):
                        val = v.text
                cells[_col_index(ref)] = val
            if cells:
                width = max(cells) + 1
                rows.append([cells.get(i) for i in range(width)])

        if not rows:
            sheets[name] = pd.DataFrame()
            continue
        width = max(len(r) for r in rows)
        rows = [r + [None] * (width - len(r)) for r in rows]
        header = [str(h) if h is not None else f"col{i}"
                  for i, h in enumerate(rows[0])]
        sheets[name] = pd.DataFrame(rows[1:], columns=header)
    return sheets


def main() -> int:
    files = {
        "S0033291719002745sup003.xlsx": "File_S1_DEGs",
        "S0033291719002745sup001.xlsx": "File_S2_WGCNA_modules",
    }
    summary = []
    for fname, label in files.items():
        path = os.path.join(SUP, fname)
        if not os.path.exists(path):
            print(f"MISSING: {path}", file=sys.stderr)
            return 1
        print(f"\n=== {label}  ({fname}) ===")
        for sheet, df in read_sheets(path).items():
            safe = re.sub(r"[^A-Za-z0-9_]+", "_", sheet).strip("_")
            dest = os.path.join(OUT, f"{label}__{safe}.csv")
            df.to_csv(dest, index=False)
            print(f"  {sheet:34s} {df.shape[0]:6,} x {df.shape[1]:2d}  -> {os.path.basename(dest)}")
            print(f"      columns: {list(df.columns)[:8]}")
            summary.append((label, sheet, df.shape))
    print(f"\n{len(summary)} sheets written to data/reference/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
