#!/bin/sh
# Run a cohort's 75 splits locally, N at a time. Used when the cluster quota is
# held by another job; the cluster worker does the same work from S3.
COHORT="${1:-groupA}"; PAR="${2:-5}"
N=0
for S in $(seq 1 75); do
  OUT="results/bpd/${COHORT}_perm00_split$(printf '%03d' $S).csv"
  [ -s "$OUT" ] && continue
  Rscript scripts/11_run_bpd_split.R --split "$S" --cohort "$COHORT" --impl ".." \
    --outdir results/bpd > "results/bpd/.log_${COHORT}_$S" 2>&1 &
  N=$((N+1))
  [ $((N % PAR)) -eq 0 ] && wait
done
wait
echo "$COHORT done: $(ls results/bpd/${COHORT}_perm00_split*.csv 2>/dev/null | wc -l) of 75"
