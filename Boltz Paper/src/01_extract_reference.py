#!/usr/bin/env python3
"""Extract Boltz et al.'s published results into CSV: the validation target for everything downstream.

Two sources:
  1. Supplementary Tables S1-S3 and S13-S15 (xlsx), parsed straight from OOXML so no spreadsheet
     library is needed. S4-S12 (FUSION/TWAS, 75 MB) are kept as downloaded but not extracted: TWAS
     needs genotypes, so nothing here can be compared against them.
  2. Numbers stated only in the main text (Table 1, proportion p-values, DE counts), hand-entered
     from the PDF with page references into paper_targets.csv.
"""
from __future__ import annotations

import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SUP = os.path.join(ROOT, "data", "supplementary")
OUT = os.path.join(ROOT, "data", "reference")
os.makedirs(OUT, exist_ok=True)

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

FILES = {
    "Table_S1_DeconeQTL_eGenes (mmc2).xlsx": "S1_deconeqtl",
    "Table_S2_replication_bulk_vs_celltype (mmc3).xlsx": "S2_replication",
    "Table_S3_scRNAseq_expression_correlation (mmc4).xlsx": "S3_scrna_correlation",
    "Table_S13_lithium_interaction_eQTL (mmc14).xlsx": "S13_li_interaction",
    "Table_S14_Li-eGene_summary (mmc15).xlsx": "S14_li_egene_summary",
    "Table_S15_limma_DE_results (mmc16).xlsx": "S15",
}

# Main-text numbers. (quantity, value, where in the paper). Proportion p-values are from glm
# logistic regressions, one cell type at a time, adjusted for age, sex, RNA concentration and RIN.
TARGETS = [
    ("n_total", 1730, "p324 cohort"), ("n_bp", 1045, "p324"), ("n_scz", 84, "p324"),
    ("n_control", 601, "p324"), ("n_li_user", 709, "p330"), ("n_li_nonuser", 336, "p330"),
    ("prop_mean_b.cells.naive", 0.025, "Table 1"), ("prop_sd_b.cells.naive", 0.020, "Table 1"),
    ("prop_mean_b.cells.memory", 0.020, "Table 1"), ("prop_sd_b.cells.memory", 0.014, "Table 1"),
    ("prop_mean_t.cells.cd8", 0.025, "Table 1"), ("prop_sd_t.cells.cd8", 0.025, "Table 1"),
    ("prop_mean_t.cells.cd4.naive", 0.15, "Table 1"), ("prop_sd_t.cells.cd4.naive", 0.042, "Table 1"),
    ("prop_mean_t.cells.cd4.memory.resting", 0.066, "Table 1"),
    ("prop_sd_t.cells.cd4.memory.resting", 0.034, "Table 1"),
    ("prop_mean_nk.cells.resting", 0.066, "Table 1"), ("prop_sd_nk.cells.resting", 0.029, "Table 1"),
    ("prop_mean_monocytes", 0.050, "Table 1"), ("prop_sd_monocytes", 0.039, "Table 1"),
    ("prop_mean_neutrophils", 0.51, "Table 1"), ("prop_sd_neutrophils", 0.094, "Table 1"),
    # "CD4 T cells" in the case/control sentence is not further specified in the paper.
    ("p_casecontrol_cd4_t", 1.8e-7, "p330, higher in controls"),
    ("p_casecontrol_nk.cells.resting", 1.2e-7, "p330, higher in controls"),
    ("p_casecontrol_neutrophils", 2.3e-8, "p330, higher in cases"),
    ("p_lithium_t.cells.cd4.naive", 8e-4, "p330, higher in non-users"),
    ("p_lithium_t.cells.cd4.memory.resting", 4e-4, "p330, higher in non-users"),
    ("p_lithium_nk.cells.resting", 3e-4, "p330, higher in non-users"),
    ("p_lithium_neutrophils", 1.5e-9, "p330, higher in users"),
    ("n_sig_nonuser_vs_control", 0, "p330, no cell type significant"),
    ("n_genes_tested", 17194, "p326/p332, >=1 TPM in >=436 samples"),
    ("n_deg_lithium", 100, "p332, FDR<0.05"), ("li_logfc_min", -0.191, "p332"),
    ("li_logfc_max", 0.177, "p332"),
    ("n_deg_lithium_in_krebs", 33, "p333"), ("fisher_or_krebs", 6.43, "p333"),
    ("fisher_p_krebs", 4.74e-14, "p333"),
    ("n_deg_casecontrol", 64, "p332"), ("cc_logfc_min", -0.126, "p332"), ("cc_logfc_max", 0.104, "p332"),
    ("n_deg_overlap_cc_li", 9, "p332"),
    ("n_deg_lithium_celladj", 94, "p333"), ("n_deg_lithium_celladj_shared", 82, "p333"),
    ("bmind_cc_neutrophils", 4, "p333"), ("bmind_cc_b.cells.memory", 24, "p333"),
    ("bmind_cc_b.cells.naive", 21, "p333"), ("bmind_cc_bmem_bnaive_shared", 18, "p333"),
    ("bmind_li_any", 0, "p333, no DEG in any cell type"),
    ("r2_bmind_schmiedel_monocytes", 0.82, "Table S3"),
    ("r2_bmind_schmiedel_t.cells.cd4.naive", 0.84, "Table S3"),
    ("r2_bmind_schmiedel_b.cells.naive", 0.61, "Table S3"),
]


def _shared_strings(z: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    return ["".join(t.text or "" for t in si.iter(NS + "t")) for si in root.findall(NS + "si")]


def _col_index(ref: str) -> int:
    n = 0
    for ch in re.match(r"([A-Z]+)", ref).group(1):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_sheets(path: str) -> dict[str, list[list]]:
    """Every sheet as a list of rows (raw cells). Headers are handled per table, because several
    Boltz sheets put a title sentence in row 1 and the real header in row 2 or 3."""
    z = zipfile.ZipFile(path)
    strings = _shared_strings(z)
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels}
    sheets = {}
    for sh in wb.iter(NS + "sheet"):
        tgt = target[sh.get(REL + "id")].lstrip("/")
        tgt = tgt if tgt.startswith("xl/") else "xl/" + tgt
        rows = []
        for row in ET.fromstring(z.read(tgt)).iter(NS + "row"):
            cells = {}
            for c in row.findall(NS + "c"):
                v, t = c.find(NS + "v"), c.get("t")
                if t == "s":
                    val = strings[int(v.text)] if v is not None else None
                elif t == "inlineStr":
                    val = "".join(x.text or "" for x in c.iter(NS + "t"))
                elif v is None:
                    val = None
                else:
                    try:
                        val = float(v.text)
                    except ValueError:
                        val = v.text
                cells[_col_index(c.get("r"))] = val
            if cells:
                rows.append([cells.get(i) for i in range(max(cells) + 1)])
        width = max((len(r) for r in rows), default=0)
        sheets[sh.get("name")] = [r + [None] * (width - len(r)) for r in rows]
    return sheets


def to_frame(rows: list[list]) -> pd.DataFrame:
    """Use the first row whose first cell names a column ('gene', 'SNP', 'Cell Type', ...) as the
    header. Title and note rows above it are dropped; they are captions, not data."""
    for i, r in enumerate(rows):
        filled = [x for x in r if x not in (None, "")]
        if len(filled) >= 2 and all(isinstance(x, str) for x in filled) and len(str(r[0] or "")) < 40:
            header = [str(h) if h not in (None, "") else f"col{j}" for j, h in enumerate(r)]
            df = pd.DataFrame(rows[i + 1:], columns=header)
            return df.dropna(axis=1, how="all").dropna(axis=0, how="all")
    return pd.DataFrame(rows)


def main() -> int:
    for fname, label in FILES.items():
        path = os.path.join(SUP, fname)
        if not os.path.exists(path):
            print(f"MISSING {path} -- run download_data.sh", file=sys.stderr)
            return 1
        for sheet, rows in read_sheets(path).items():
            df = to_frame(rows)
            # S15 is the only workbook with several sheets; name each by its sheet.
            name = label if label != "S15" else "S15_" + re.sub(r"[^A-Za-z0-9]+", "_", sheet).strip("_")
            df.to_csv(os.path.join(OUT, f"{name}.csv"), index=False)
            print(f"  {name:42s} {df.shape[0]:6,} x {df.shape[1]:2d}   {list(df.columns)[:5]}")

    pd.DataFrame(TARGETS, columns=["quantity", "value", "source"]).to_csv(
        os.path.join(OUT, "paper_targets.csv"), index=False)
    print(f"  {'paper_targets':42s} {len(TARGETS):6,} main-text numbers")
    return 0


if __name__ == "__main__":
    sys.exit(main())
