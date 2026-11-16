#!/bin/sh
# One worker of the cell-type expression stage (proposal inputs 4 and 5).
#
# Each worker takes the worklist rows where (row number mod NW) equals its
# index. One row is one split: bMIND on the training fold, then separately on
# the test fold, cached to S3 as a single RDS the prediction runner reads.
#
# This stage is slow -- bMIND is MCMC, and the DE arm's single run over 474
# samples took 187 minutes -- so an item that is already in the bucket is
# skipped and a re-run resumes rather than restarting.
set -eu

NW="${NW:-12}"
IDX="${VC_TASK_INDEX:-${HOSTNAME##*-}}"
IDX=$(echo "$IDX" | sed 's/^0*//'); [ -n "$IDX" ] || IDX=0

W=/work
mkdir -p "$W/prediction/cache/ct" "$W/prediction/results/ct" "$W/prediction/splits" "$W/scripts"
# Which permutations to score for inputs 4 and 5. 0 is the real run;
# 1..50 are the null and are applied to the first 3 splits only.
PERMS="${PERMS:-$(seq 0 50)}"
export R_LIBS_USER="${R_LIBS_USER:-/usr/local/lib/R/site-library}"
export BMIND_CORES="${BMIND_CORES:-2}"

echo "worker $IDX of $NW on $(hostname), $(nproc) cores, BMIND_CORES=$BMIND_CORES"

echo "=== fetching inputs ==="
IN=prediction/bpd/inputs
python3 /cfg/s3.py get "$IN/lithium_inputs.rds"     "$W/prediction/cache/lithium_inputs.rds"
python3 /cfg/s3.py get "$IN/bpd_inputs.rds"         "$W/prediction/cache/bpd_inputs.rds"
python3 /cfg/s3.py get "$IN/splits_bpd_lithium.csv" "$W/prediction/splits/splits_bpd_lithium.csv"
python3 /cfg/s3.py get "$IN/splits_bpd_groupA.csv"  "$W/prediction/splits/splits_bpd_groupA.csv"
python3 /cfg/s3.py get "$IN/${WORKLIST:-worklist_bmind.tsv}" "$W/worklist.tsv"
python3 /cfg/s3.py get "$IN/19_bmind_split.R"       "$W/scripts/19_bmind_split.R"
python3 /cfg/s3.py get "$IN/21_ct_predict.R"        "$W/scripts/21_ct_predict.R"

echo "=== R packages ==="
# MIND is NOT on CRAN -- it is randel/MIND on GitHub, and the first smoke run
# failed with "there is no package called MIND" after BiocManager looked for it
# in the wrong place. Its dependencies ARE on CRAN and Bioconductor, so those
# are installed first from Posit binaries (fast) and only the package itself
# comes from source.
# Output goes to a file and the exit status is checked directly. Piping
# Rscript into `tail` makes the pipeline report TAIL's status, so the first
# version of this guard never fired and the worker walked into 75 identical
# failures.
#
# errexit is suspended around it for the second half of that bug: with
# `set -e` in force, a non-zero Rscript kills the script on the spot, so
# PKGOK=$? never runs and the diagnostics below never print. The symptom is a
# worker whose log simply stops after "=== R packages ===".
set +e
Rscript -e '
  options(repos = c(CRAN = "https://p3m.dev/cran/__linux__/jammy/latest"),
          Ncpus = max(1L, as.integer(Sys.getenv("BMIND_CORES", "2"))))
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  # dependencies=TRUE on the GitHub install is what the first attempt got
  # wrong: MIND imports Biobase and edgeR from Bioconductor, and with
  # dependencies=FALSE its compile fails with only "non-zero exit status".
  bioc <- c("edgeR","Biobase")
  bmiss <- bioc[!sapply(bioc, requireNamespace, quietly=TRUE)]
  if (length(bmiss)) BiocManager::install(bmiss, ask=FALSE, update=FALSE)
  # BisqueRNA'\''s own chain is listed explicitly because it is installed from a
  # tarball below, and install.packages(repos=NULL) resolves NO dependencies --
  # it fails one at a time, a cluster round-trip per missing package.
  # glmnet and randomForest are for 21_ct_predict.R, which runs in this same
  # worker right after each deconvolution. Leaving them out meant every
  # profile computed correctly and every prediction died on "there is no
  # package called glmnet" -- a 15-hour job that would have produced no
  # results at all.
  deps <- c("nnls","doParallel","foreach","MCMCglmm","Matrix","matrixcalc","remotes",
            "limSolve","quadprog","lpSolve","MASS","glmnet","randomForest")
  miss <- deps[!sapply(deps, requireNamespace, quietly=TRUE)]
  if (length(miss)) install.packages(miss)
  # BisqueRNA has been ARCHIVED from CRAN, so a plain install.packages finds
  # nothing and MIND then fails to compile with only "dependency BisqueRNA is
  # not available". Installed from the CRAN archive at a pinned version, which
  # is also better for reproducibility than whatever HEAD happens to be.
  if (!requireNamespace("BisqueRNA", quietly=TRUE))
    install.packages(
      "https://cran.r-project.org/src/contrib/Archive/BisqueRNA/BisqueRNA_1.0.5.tar.gz",
      repos = NULL, type = "source")
  for (p in c(deps, bioc))
    if (!requireNamespace(p, quietly=TRUE)) cat("  MISSING DEP:", p, "\n")
  if (!requireNamespace("MIND", quietly=TRUE))
    remotes::install_github("randel/MIND", upgrade="never", dependencies=TRUE)
  ok <- requireNamespace("MIND", quietly=TRUE)
  cat(sprintf("  MIND %s version %s | edgeR %s\n", ok,
              if (ok) as.character(packageVersion("MIND")) else "-",
              requireNamespace("edgeR", quietly=TRUE)))
  pred <- c("glmnet","randomForest")
  pmiss <- pred[!sapply(pred, requireNamespace, quietly=TRUE)]
  if (length(pmiss)) cat("  MISSING PREDICTION DEP:", paste(pmiss, collapse=", "), "
")
  if (!ok || length(pmiss)) quit(status = 1)
' > "$W/pkg.log" 2>&1
PKGOK=$?
set -e
tail -5 "$W/pkg.log"
if [ "$PKGOK" -ne 0 ]; then
  echo "!! MIND unavailable, aborting worker"
  # A pod and its logs are deleted a minute after the job ends, so the log has
  # to leave the pod to be readable at all. The previous failure was
  # undiagnosable for exactly this reason.
  python3 /cfg/s3.py put "$W/pkg.log" "prediction/bpd/logs/pkg_${IDX}.log" >/dev/null 2>&1 \
    && echo "   install log -> prediction/bpd/logs/pkg_${IDX}.log"
  grep -iE "error|cannot|unable|non-zero|MISSING DEP|there is no package" "$W/pkg.log" | tail -20
  exit 1
fi

echo "=== already done ==="
python3 /cfg/s3.py ls prediction/bpd/ct/        >  "$W/done.txt" 2>/dev/null || : > "$W/done.txt"
python3 /cfg/s3.py ls prediction/bpd/ctresults/ >> "$W/done.txt" 2>/dev/null || :
echo "  $(wc -l < "$W/done.txt") present"

N=0
TAB=$(printf '	')
tail -n +2 "$W/worklist.tsv" | while IFS="$TAB" read -r COHORT SPLIT NGENES; do
  N=$((N + 1))
  [ $(( (N - 1) % NW )) -eq "$IDX" ] || continue
  # The gene count is in the key. Without it, a 400-gene smoke run and a
  # 4,000-gene production run collide: the production worker finds the smoke
  # file, skips the item, and one split in 75 silently carries a different
  # gene universe from the other 74. That happened on the first dispatch.
  KEY="prediction/bpd/ct/ct_${COHORT}_g${NGENES}_split$(printf '%03d' "$SPLIT").rds"
  OUT="$W/prediction/cache/ct/$(basename "$KEY")"
  if grep -qx "$KEY" "$W/done.txt" 2>/dev/null; then
    echo "--- item $N: $KEY present, skipping ---"; continue
  fi
  echo "--- item $N: $COHORT split=$SPLIT ngenes=$NGENES ---"
  if Rscript "$W/scripts/19_bmind_split.R" \
       --split "$SPLIT" --cohort "$COHORT" --ngenes "$NGENES" \
       --impl "$W" --outdir "$W/prediction/cache/ct" 2>&1 | grep -vE "^Loading|^Warning|built under"; then
    # Size guard: an RDS under 10 KB cannot hold a samples x (genes x types)
    # matrix, so it means the run produced nothing usable.
    if [ -s "$OUT" ] && [ "$(wc -c < "$OUT")" -gt 10240 ]; then
      python3 /cfg/s3.py put "$OUT" "$KEY" >/dev/null && echo "    -> $KEY"

      # Score inputs 4 and 5 HERE, in the same worker, while the profile is
      # still on local disk. Deconvolution is the expensive step and the
      # prediction on top of it is seconds, so splitting them into two jobs
      # would mean a second job that has to wait for this one -- and therefore
      # something outside the cluster to start it. This way the arm finishes
      # on its own.
      for P in $PERMS; do
        # Permutations only need the prediction step: shuffling labels does
        # not change a deconvolution. They run on the first 3 splits only,
        # matching the 50 x 3 null used for inputs 1-3.
        if [ "$P" -gt 0 ] && [ "$SPLIT" -gt 3 ]; then continue; fi
        RKEY="prediction/bpd/ctresults/ct_${COHORT}_g${NGENES}_perm$(printf '%02d' "$P")_split$(printf '%03d' "$SPLIT").csv"
        ROUT="$W/prediction/results/ct/$(basename "$RKEY")"
        if grep -qx "$RKEY" "$W/done.txt" 2>/dev/null; then continue; fi
        if Rscript "$W/scripts/21_ct_predict.R" \
             --split "$SPLIT" --cohort "$COHORT" --ngenes "$NGENES" --perm "$P" \
             --models enet,rf --impl "$W" --ctdir "$W/prediction/cache/ct" \
             --outdir "$W/prediction/results/ct" 2>&1 | grep -vE "^Loading|^Warning|built under"; then
          if [ -s "$ROUT" ] && [ "$(wc -l < "$ROUT")" -gt 1 ]; then
            python3 /cfg/s3.py put "$ROUT" "$RKEY" >/dev/null && echo "    -> $RKEY"
          else
            echo "    !! perm $P produced no usable rows"
          fi
        else
          echo "    !! perm $P failed, continuing"
        fi
        rm -f "$ROUT"
      done

      rm -f "$OUT"            # 75 of these would fill the 20Gi node disk cap
    else
      echo "    !! output too small ($(wc -c < "$OUT" 2>/dev/null || echo 0) bytes); not uploading"
    fi
  else
    echo "    !! item failed, continuing"
  fi
done

echo "worker $IDX finished"
