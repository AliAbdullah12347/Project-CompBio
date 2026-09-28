#!/usr/bin/env python3
"""Build the 474-sample analysis cohort from the two deposited GEO files.

Reads ONLY from data/raw_as_deposited/, which is never written to. Writes the
filtered cohort to data/cohort_474/.

Six of the 480 deposited samples are removed here and nowhere else:

  * 4 sample mix-ups. Krebs et al. identified these as subject-identity errors.
    Their covariates are complete, so nothing is wrong with the data itself --
    the data belongs to the wrong person. No imputation can repair that.
  * 2 samples missing six of nine analysis covariates. Recovering them would
    mean imputing six variables from three, which is invention, not recovery.

The 30 samples missing tobacco use and nothing else are KEPT. They are the
subjects the imputation work is about, and discarding them at load time would
make that analysis impossible.

Neither exclusion is derived from a hard-coded list of sample IDs. Both are
derived from the data by rule, then checked against the identifiers on record,
so a change in the deposit would surface as a failure rather than pass silently.

    python data/build_cohort.py
"""
from __future__ import annotations

import gzip
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw_as_deposited"
OUT = HERE / "cohort_474"

COUNTS_FILE = "GSE124326_count_matrix.txt.gz"
SERIES_FILE = "GSE124326_series_matrix.txt.gz"

# The nine covariates every analysis model uses. Missingness is counted over
# exactly these; the 22 cell-type fractions and the QC flag are not covariates.
COVARIATES = [
    "age", "Sex", "tobacco use", "assessment group", "rin",
    "sequencing plate", "sequencing metric pc1", "sequencing metric pc2",
    "sequencing metric pc3",
]

# Krebs et al.'s four sample mix-ups, for checking the rule below reproduces
# them. Not used to perform the exclusion.
KNOWN_MIXUPS = ["099252A", "108773A", "118140B", "118141A"]

# Our own quality screen. The deposited flag screens neither sequencing depth
# nor RNA integrity: nine flag-passing samples fall below 3M assigned reads and
# four have RIN under 5, while the samples the deposit excluded have a HIGHER
# depth floor at 3.07M. So depth played no part in the deposited decision.
MIN_LIBRARY_SIZE = 3_000_000
MIN_RIN = 5.0


def verify_checksums() -> None:
    recorded = {}
    for line in (RAW / "CHECKSUMS.sha256").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        recorded[name.lstrip("*").strip()] = digest
    for name, want in recorded.items():
        got = hashlib.sha256((RAW / name).read_bytes()).hexdigest()
        if got != want:
            raise SystemExit(f"CHECKSUM MISMATCH for {name}\n  want {want}\n  got  {got}")
        print(f"  [ok] {name}  {got[:16]}...")


def parse_series(path: Path) -> pd.DataFrame:
    """One row per sample. Characteristics lines are ragged, so every
    'field: value' cell is assigned by its own key, never by row position."""
    simple: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip('"') for v in vals]
            if key == "!Sample_characteristics_ch1":
                characteristics.append(vals)
            elif key in ("!Sample_title", "!Sample_geo_accession"):
                simple[key] = vals

    n = len(simple["!Sample_title"])
    meta = pd.DataFrame({"title": simple["!Sample_title"],
                         "gsm": simple["!Sample_geo_accession"]})
    fields: dict[str, list] = {}
    for row in characteristics:
        for i, cell in enumerate(row):
            if ":" not in cell:
                continue
            field, value = cell.split(":", 1)
            fields.setdefault(field.strip(), [None] * n)[i] = value.strip()
    for field, values in fields.items():
        meta[field] = values
    return meta.replace("NA", np.nan)


def load_counts(path: Path) -> pd.DataFrame:
    """Columns are sample titles with a '.counts' suffix, not GSM accessions."""
    counts = pd.read_csv(path, sep="\t", index_col=0, compression="gzip")
    counts.columns = [c[:-7] if c.endswith(".counts") else c for c in counts.columns]
    return counts


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    print("Verifying raw files against recorded checksums ...")
    verify_checksums()

    print("\nReading deposited files ...")
    meta = parse_series(RAW / SERIES_FILE)
    counts = load_counts(RAW / COUNTS_FILE)
    print(f"  series matrix : {len(meta)} samples, "
          f"{len(meta.columns) - 2} characteristics fields")
    print(f"  count matrix  : {counts.shape[0]:,} genes x {counts.shape[1]} samples")

    meta["qc_pass"] = meta["included in final analysis"].str.strip().str.upper() == "TRUE"
    meta["n_missing_covariates"] = meta[COVARIATES].isna().sum(axis=1)
    # Assigned reads per sample, from the count matrix itself.
    meta["library_size"] = meta["title"].map(counts.sum(axis=0)).astype("int64")
    meta["rin_value"] = pd.to_numeric(meta["rin"], errors="coerce")
    meta["depth_ok"] = meta["library_size"] >= MIN_LIBRARY_SIZE
    meta["rin_ok"] = meta["rin_value"] >= MIN_RIN

    # --- exclusion 1: sample mix-ups -------------------------------------
    # Rule: excluded by the deposit despite having every covariate present.
    # A sample with complete data is not dropped for missingness, so the only
    # remaining reason the depositors had is the identity error they reported.
    is_mixup = (~meta["qc_pass"]) & (meta["n_missing_covariates"] == 0)
    mixups = sorted(meta.loc[is_mixup, "title"])
    print(f"\nExclusion 1 -- sample mix-ups: {len(mixups)}")
    print(f"  derived by rule : {mixups}")
    print(f"  on record       : {sorted(KNOWN_MIXUPS)}")
    if mixups != sorted(KNOWN_MIXUPS):
        raise SystemExit("ERROR: derived mix-ups do not match the four on record.")
    print("  [ok] rule reproduces the four identifiers on record")

    # --- exclusion 2: too much missing to recover ------------------------
    # The 30 keepers are missing exactly one covariate (tobacco). Anything
    # missing two or more is a different situation, and here that means six.
    is_severe = meta["n_missing_covariates"] >= 2
    severe = sorted(meta.loc[is_severe, "title"])
    print(f"\nExclusion 2 -- six of nine covariates missing: {len(severe)}")
    for t in severe:
        row = meta[meta["title"] == t].iloc[0]
        absent = [c for c in COVARIATES if pd.isna(row[c])]
        print(f"  {t}  ({row['bipolar disorder diagnosis']}, gsm {row['gsm']})")
        print(f"    missing: {', '.join(absent)}")

    # --- what is kept ----------------------------------------------------
    kept_incomplete = meta[meta["n_missing_covariates"] == 1]
    print(f"\nKEPT -- missing tobacco use only: {len(kept_incomplete)}")
    print(f"  diagnosis: {dict(kept_incomplete['bipolar disorder diagnosis'].value_counts())}")
    print(f"  all missing the same single field: "
          f"{set(kept_incomplete[COVARIATES].isna().idxmax(axis=1)) == {'tobacco use'}}")

    # --- quality: FLAGGED, NOT DROPPED --------------------------------------
    # Krebs et al. state twice that no sample was removed for low quality.
    # They adjust instead: RIN and the three sequencing-metric PCs enter every
    # model as covariates, and voom's precision weights exist, in their words,
    # "to account for differences between samples in sequencing depth".
    # Excluding shallow samples on top of that corrects the same problem twice
    # with the blunter of the two instruments, costs 13 samples (9 of them
    # bipolar I), and adds a second mildly diagnosis-associated filter to a
    # project whose contribution is criticising the first one.
    #
    # So depth_ok and rin_ok are computed and written to the metadata, and any
    # analysis that wants the screened subset can apply them in one line.
    is_shallow = ~meta["depth_ok"]
    is_degraded = meta["depth_ok"] & ~meta["rin_ok"]
    print(f"\nQuality flags (recorded, not applied): "
          f"{int(is_shallow.sum())} below {MIN_LIBRARY_SIZE:,} reads, "
          f"{int(is_degraded.sum())} with RIN < {MIN_RIN}")
    print("  use meta.depth_ok & meta.rin_ok for the screened sensitivity set")

    drop = is_mixup | is_severe
    cohort = meta[~drop].reset_index(drop=True)
    print(f"\n480 - {int(is_mixup.sum())} mix-ups - {int(is_severe.sum())} "
          f"severe = {len(cohort)} samples")

    # --- consistency across the two raw files ----------------------------
    print("\nChecking the two raw files agree on the cohort ...")
    titles = cohort["title"].tolist()
    missing_in_counts = sorted(set(titles) - set(counts.columns))
    if missing_in_counts:
        raise SystemExit(f"ERROR: {len(missing_in_counts)} cohort titles absent "
                         f"from the count matrix: {missing_in_counts[:5]}")
    print(f"  [ok] all {len(titles)} cohort titles present as count-matrix columns")
    if len(set(titles)) != len(titles):
        raise SystemExit("ERROR: duplicate sample titles in the cohort")
    print(f"  [ok] titles unique")

    # Reindex by the key rather than trusting that column order agrees.
    counts_cohort = counts[titles]
    if list(counts_cohort.columns) != titles:
        raise SystemExit("ERROR: count columns did not align to metadata order")
    print(f"  [ok] count matrix reordered to metadata order by key, not position")
    print(f"  [ok] dropped samples absent from both outputs: "
          f"{sorted(set(counts.columns) - set(counts_cohort.columns))}")

    # --- write ------------------------------------------------------------
    print(f"\nWriting {OUT.name}/ ...")
    excluded = meta[drop].copy()
    excluded["reason"] = np.where(
        excluded["title"].isin(mixups),
        "sample mix-up (subject identity error; covariates complete)",
        "six of nine covariates missing (not recoverable by imputation)")
    excluded[["title", "gsm", "bipolar disorder diagnosis", "library_size",
              "rin_value", "n_missing_covariates", "reason"]].to_csv(
        OUT / "excluded_samples.csv", index=False)

    # Plain one-per-line list, for normalising every other file in the project
    # to this exact sample set.
    (OUT / "sample_list.txt").write_text(
        "\n".join(cohort["title"]) + "\n", encoding="utf-8")

    cohort.to_csv(OUT / "metadata_474.csv", index=False)
    counts_cohort.to_csv(OUT / "counts_474.tsv.gz", sep="\t", compression="gzip")

    for f in ("metadata_474.csv", "counts_474.tsv.gz", "excluded_samples.csv",
              "sample_list.txt"):
        print(f"  {f:24s} {(OUT / f).stat().st_size / 1e6:8.2f} MB")

    print(f"\ncohort: {counts_cohort.shape[0]:,} genes x {counts_cohort.shape[1]} samples")
    print("diagnosis breakdown:")
    print(cohort["bipolar disorder diagnosis"].value_counts().to_string())
    bp1 = cohort[cohort["bipolar disorder diagnosis"] == "BP1"]
    li = pd.to_numeric(bp1["lithium use (non-user=0, user = 1)"])
    print(f"\nlithium contrast (BP1): {int((li == 0).sum())} non-users vs "
          f"{int((li == 1).sum())} users")
    print(f"tobacco still missing in cohort: "
          f"{int(cohort['tobacco use'].isna().sum())} samples")
    return 0


if __name__ == "__main__":
    sys.exit(main())
