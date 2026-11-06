#!/bin/sh
# One worker of the bipolar-status contest.
#
# No PVC: the inputs are 15 MB, so each pod pulls its own copy from S3 into a
# pod-local emptyDir. Scratch volumes here are reaped with their claim -- one
# already vanished mid-run -- and at this size the volume buys nothing.
#
# Each worker takes the worklist rows where (row number mod NW) equals its own
# index, so the 375 items divide evenly with no coordination. Results go
# straight back to S3 after each item, so an interrupted worker loses at most
# the item it was on.
set -eu

NW="${NW:-6}"
IDX="${VC_TASK_INDEX:-${HOSTNAME##*-}}"
# dash has no 10# base prefix; strip leading zeros portably.
IDX=$(echo "$IDX" | sed 's/^0*//'); [ -n "$IDX" ] || IDX=0

W=/work
mkdir -p "$W/prediction/cache" "$W/prediction/splits" "$W/prediction/results/bpd" "$W/scripts"
export R_LIBS_USER="${R_LIBS_USER:-/usr/local/lib/R/site-library}"

echo "worker $IDX of $NW on $(hostname), $(nproc) cores"

echo "=== fetching inputs from S3 ==="
IN=prediction/bpd/inputs
python3 /cfg/s3.py get "$IN/bpd_inputs.rds"              "$W/prediction/cache/bpd_inputs.rds"
python3 /cfg/s3.py get "$IN/splits_bpd_groupA.csv"       "$W/prediction/splits/splits_bpd_groupA.csv"
python3 /cfg/s3.py get "$IN/splits_bpd_full.csv"         "$W/prediction/splits/splits_bpd_full.csv"
python3 /cfg/s3.py get "$IN/splits_bpd_groupA_plate.csv" "$W/prediction/splits/splits_bpd_groupA_plate.csv"
python3 /cfg/s3.py get "$IN/worklist_bpd.tsv"            "$W/worklist.tsv"
python3 /cfg/s3.py get "$IN/11_run_bpd_split.R"          "$W/scripts/11_run_bpd_split.R"

echo "=== R packages ==="
Rscript -e '
  need <- c("edgeR","glmnet","randomForest")
  miss <- need[!sapply(need, requireNamespace, quietly=TRUE)]
  if (length(miss)) {
    # Posit binaries for this distro: minutes faster than compiling from source,
    # and six pods each compiling edgeR is minutes nobody needs to spend.
    options(repos = c(CRAN = "https://p3m.dev/cran/__linux__/jammy/latest"))
    if ("edgeR" %in% miss) {
      if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
      BiocManager::install("edgeR", ask=FALSE, update=FALSE)
    }
    cran <- setdiff(miss, "edgeR")
    if (length(cran)) install.packages(cran)
  }
  for (p in need) cat(sprintf("  %-14s %s\n", p, requireNamespace(p, quietly=TRUE)))
' 2>&1 | tail -6

echo "=== items already in S3 (skipped) ==="
python3 /cfg/s3.py ls prediction/bpd/results/ > "$W/done.txt" 2>/dev/null || : > "$W/done.txt"
echo "  $(wc -l < "$W/done.txt") already present"

echo "=== work items ==="
N=0
DONE=0
# Tab-separated, and IFS is a literal tab: the models and features fields
# contain commas, and splitting on comma tore "enet,rf" into two fields -- the
# whole first cluster run produced zero-row files that way.
TAB=$(printf '	')
# Skip the header, then take every NWth row.
tail -n +2 "$W/worklist.tsv" | while IFS="$TAB" read -r COHORT PERM SPLIT MODELS FEATURES; do
  N=$((N + 1))
  [ $(( (N - 1) % NW )) -eq "$IDX" ] || continue
  KEY="prediction/bpd/results/${COHORT}_perm$(printf '%02d' "$PERM")_split$(printf '%03d' "$SPLIT").csv"
  OUT="$W/prediction/results/bpd/$(basename "$KEY")"
  if grep -qx "$KEY" "$W/done.txt" 2>/dev/null; then
    echo "--- item $N: $KEY already in S3, skipping ---"; continue
  fi
  echo "--- item $N: $COHORT perm=$PERM split=$SPLIT [$MODELS / $FEATURES] ---"
  if Rscript "$W/scripts/11_run_bpd_split.R" \
       --split "$SPLIT" --cohort "$COHORT" --perm "$PERM" \
       --models "$MODELS" --features "$FEATURES" \
       --impl "$W" --outdir "$W/prediction/results/bpd" 2>&1 | grep -vE "^Loading|^Warning|built under"; then
    if [ -s "$OUT" ] && [ "$(wc -l < "$OUT")" -gt 1 ]; then
      python3 /cfg/s3.py put "$OUT" "$KEY" >/dev/null && echo "    -> $KEY"
      DONE=$((DONE + 1))
    else
      echo "    !! no usable output ($(wc -l < "$OUT" 2>/dev/null || echo 0) lines); not uploading"
    fi
  else
    echo "    !! item failed, continuing"
  fi
done

echo "worker $IDX finished"
