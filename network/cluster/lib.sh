# Shared helpers for the Arm 3 cluster jobs. Sourced by every job's command.
#
# S3 is the only durable store: a pod's files and logs are deleted a minute
# after the job ends. curl --aws-sigv4 is used rather than the AWS CLI so no
# pip install is needed; Content-Type is set on every PUT because Ceph reads a
# body without one as a browser form and refuses it.
set -euo pipefail

: "${BUCKET:?BUCKET not set}" "${PREFIX:?PREFIX not set}"
S3="${AWS_ENDPOINT_URL}/${BUCKET}"
CA=/etc/s3-ca/ca-bundle.pem

s3_put () {   # s3_put <localfile> <key-under-PREFIX>
  curl --fail --silent --show-error --retry 3 --retry-delay 5 --cacert "$CA" \
    --aws-sigv4 "aws:amz:${AWS_REGION}:s3" \
    --user "${AWS_ACCESS_KEY_ID}:${AWS_SECRET_ACCESS_KEY}" \
    -H "Content-Type: application/octet-stream" -T "$1" "${S3}/${PREFIX}/$2"
}
s3_put_abs () {   # s3_put_abs <localfile> <absolute key>
  curl --fail --silent --show-error --retry 3 --retry-delay 5 --cacert "$CA" \
    --aws-sigv4 "aws:amz:${AWS_REGION}:s3" \
    --user "${AWS_ACCESS_KEY_ID}:${AWS_SECRET_ACCESS_KEY}" \
    -H "Content-Type: application/octet-stream" -T "$1" "${S3}/$2"
}
s3_get () {   # s3_get <absolute key> <localfile>   (non-fatal: returns 1 if absent)
  curl --fail --silent --show-error --retry 2 --cacert "$CA" \
    --aws-sigv4 "aws:amz:${AWS_REGION}:s3" \
    --user "${AWS_ACCESS_KEY_ID}:${AWS_SECRET_ACCESS_KEY}" \
    -o "$2" "${S3}/$1"
}

# Sparse checkout: the whole repo is ~420 MB and almost all of it is results
# this job does not read. These paths are what network/00_load.R and the
# experiment scripts actually open.
clone_repo () {
  git clone --filter=blob:none --no-checkout --depth 1 \
    https://github.com/AliAbdullah12347/Project-CompBio.git /work/repo
  cd /work/repo
  git sparse-checkout init --no-cone
  git sparse-checkout set \
    "/data/analysis_matrices/*" "/data/LM22.txt" "/de_analysis/data/*" "/de_analysis/scripts/*" \
    "/network/*" "/de_v2/scripts/*" "/de_final/results/*"
  git checkout
  mkdir -p network/cache network/runs
}

# One pod builds the R library and uploads it; the rest untar it in ~1 min
# instead of compiling WGCNA from source for ten.
RLIB=/work/rlibs
setup_r () {
  mkdir -p "$RLIB"
  if s3_get "${PREFIX}/env/rlibs.tar.gz" /work/rlibs.tar.gz 2>/dev/null; then
    echo "-- reusing prebuilt R library from S3"
    tar xzf /work/rlibs.tar.gz -C "$RLIB"
  else
    echo "-- building R library (first pod only)"
    Rscript -e 'options(Ncpus=16); BiocManager::install(c("WGCNA","impute","preprocessCore","GO.db","AnnotationDbi","org.Hs.eg.db","igraph","fastcluster"), ask=FALSE, update=FALSE)'
    tar czf /work/rlibs.tar.gz -C "$RLIB" .
    s3_put /work/rlibs.tar.gz env/rlibs.tar.gz || true
  fi
  # The cached tarball predates inp-024's need for limma. Top up what is
  # missing rather than invalidating the whole tarball for one package.
  Rscript -e 'need <- setdiff(c("WGCNA","GO.db","org.Hs.eg.db","igraph","fastcluster",
                                "limma","edgeR","statmod"),
                              rownames(installed.packages()))
              if (length(need)) { cat("-- installing missing:", need, fill=TRUE)
                BiocManager::install(need, ask=FALSE, update=FALSE) }
              # bMIND. Not on CRAN or Bioconductor - de_analysis/scripts/06_bmind_profiles.R
              # does library(MIND), so the project already depends on it; only the
              # cluster image lacked it. Non-fatal: every row except bmind-046/047
              # runs fine without it.
              if (!requireNamespace("MIND", quietly=TRUE)) {
                if (!requireNamespace("remotes", quietly=TRUE))
                  install.packages("remotes", repos="https://cloud.r-project.org")
                # MIND is pure R but its dependencies are not all pulled in by
                # install_github; install them first and leave the log visible,
                # because quiet=TRUE hid the real compile error on the first try.
                BiocManager::install(setdiff(c("MCMCpack","nnls","matrixStats","Matrix",
                                               "doParallel","foreach","gtools","coda"),
                                             rownames(installed.packages())),
                                     ask=FALSE, update=FALSE)
                # BisqueRNA is archived on CRAN, so MIND cannot resolve it from the
                # normal repos - that is the real reason the build failed, not a
                # compile error. Take it from the archive first.
                # limSolve is what BisqueRNA needs and nothing pulls in, so the
                # chain is MIND -> BisqueRNA -> limSolve and each layer only
                # names the next once the one below it is present.
                if (!requireNamespace("limSolve", quietly=TRUE))
                  try(install.packages("limSolve", repos="https://cloud.r-project.org"))
                if (!requireNamespace("BisqueRNA", quietly=TRUE))
                  try(install.packages(
                    "https://cran.r-project.org/src/contrib/Archive/BisqueRNA/BisqueRNA_1.0.5.tar.gz",
                    repos=NULL, type="source"))
                try(remotes::install_github("randel/MIND", upgrade="never", quiet=FALSE))
              }'
  Rscript -e 'for (p in c("WGCNA","GO.db","org.Hs.eg.db","igraph","limma","edgeR","MIND")) cat(sprintf("%-14s %s\n", p, requireNamespace(p, quietly=TRUE)))'
}

# The base-001 reference modules. Seven experiment scripts quit(status=1) when
# this file is absent, and it is gitignored, so it must come from S3.
fetch_cache () {
  if s3_get "${PREFIX}/cache/ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds" \
            network/cache/ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds; then
    echo "-- reference cache fetched"
  else
    echo "!! reference cache missing in S3 — run arm3-prep first" >&2; return 1
  fi
  s3_get "${PREFIX}/cache/project.rds" network/cache/project.rds || true
}

# Push a finished row's outputs. Checkpoints (.rds/.RData) go up too so a
# re-run resumes instead of recomputing.
# A sliced row runs on several pods at once, each holding a different part of
# the work, and every one of them would PUT to runs/<id>.tar.gz -- last writer
# wins, and nine slices of ten vanish with nothing in any log to say so. Give
# each slice its own key; collect.py unpacks them all into one directory,
# which is what the row's own aggregation step expects to find.
run_key () {   # run_key <run_id>
  if [ "${SLICE_ID:-0}" -gt 0 ] 2>/dev/null; then
    printf 'runs/%s.slice%02d.tar.gz' "$1" "${SLICE_ID}"
  else
    printf 'runs/%s.tar.gz' "$1"
  fi
}

upload_run () {   # upload_run <run_id>
  local id="$1" dir="network/runs/$1" key
  [ -d "$dir" ] || { echo "-- $id: no output dir"; return 0; }
  key=$(run_key "$id")
  tar czf "/work/${id}.tar.gz" -C network/runs "$id"
  s3_put "/work/${id}.tar.gz" "$key"
  echo "-- $id: uploaded $(du -h "/work/${id}.tar.gz" | cut -f1) to ${key}"
}

# ---------------------------------------------------------------------------
# Correct detectCores() for the cgroup.
#
# parallel::detectCores() reports the HOST's cores - 32 on these nodes - not
# the 4 this pod may use. WGCNA's enableWGCNAThreads() then takes nCores-1, so
# every row of the 8 Nov run logged "up to 31 working processes" while holding
# a 4-CPU quota and spent its time context-switching: 17-25 minutes per
# modulePreservation call, against the ~80 s per 50 permutations recorded in
# gene-027's header. Reporting the cgroup's allocation is what detectCores()
# should do; R has no way to see it, so we supply it.
#
# R reads .Rprofile from the working directory at startup, before the script
# body, so this is in place before library(WGCNA) runs.
# ---------------------------------------------------------------------------
write_rprofile () {   # write_rprofile <cpu-limit>
  cat > /work/repo/.Rprofile <<RPROF
local({
  n  <- ${1}L
  ok <- try({
    ns <- asNamespace("parallel")
    unlockBinding("detectCores", ns)
    assign("detectCores", function(all.tests = FALSE, logical = TRUE) n + 1L, envir = ns)
    lockBinding("detectCores", ns)
    TRUE
  }, silent = TRUE)
  cat(sprintf("[.Rprofile] detectCores -> %d so enableWGCNAThreads picks %d: %s\n",
              n + 1L, n, if (isTRUE(ok)) "ok" else "FAILED, threads left at host count"))
})
RPROF
  echo "-- .Rprofile written, threads pinned to ${1}"
}

# ---------------------------------------------------------------------------
# Restore a row's earlier partial output before re-running it.
#
# Every long row checkpoints per draw (lith_draw_NNN.rds) and skips draws whose
# checkpoint already exists. The 8 Nov run timed out inside the draw loop but
# uploaded what it had, so restoring first turns a 20-draw re-run into a 2- to
# 16-draw one. Without this a resubmission silently repeats hours of finished
# work and times out in the same place.
# ---------------------------------------------------------------------------
restore_run () {   # restore_run <run_id>
  local id="$1" n d
  if s3_get "${PREFIX}/$(run_key "$id")" "/work/${id}.prev.tar.gz" 2>/dev/null; then
    tar xzf "/work/${id}.prev.tar.gz" -C network/runs || true
    # find, not ls+glob: with set -e, `d=$(ls nomatch*)` is a FAILED assignment
    # and kills the whole job. Rows that had draw checkpoints survived it and
    # rows that did not were killed the instant they were restored, which is
    # exactly the 11 rows that kept failing on 9 Nov. find exits 0 on no match.
    n=$(find "network/runs/${id}" -maxdepth 1 -type f 2>/dev/null | wc -l)
    d=$(find "network/runs/${id}" -maxdepth 1 -name '*draw_*.rds' 2>/dev/null | wc -l)
    echo "-- ${id}: restored ${n} files, ${d} draw checkpoints will be skipped"
  else
    echo "-- ${id}: nothing in S3, starting clean"
  fi
}

# filt-044 reads the 474-sample raw count matrix: 16 MB, untracked, so it
# travels through S3 the same way the reference cache does.
fetch_counts () {
  mkdir -p data/cohort_474
  if s3_get "${PREFIX}/data/counts_474.tsv.gz" data/cohort_474/counts_474.tsv.gz; then
    echo "-- raw counts fetched ($(du -h data/cohort_474/counts_474.tsv.gz | cut -f1))"
  else
    echo "!! counts_474.tsv.gz missing in S3" >&2; return 1
  fi
}

# bmind-047 reads per-lineage deconvolved expression. The source object is
# de_analysis/results/bmind_profiles.rds, 441 MB and untracked; 047 uses only
# its posterior mean, so network/cluster/prep_bmind.R reshapes that into a
# 209 MB per-lineage file and this ships it. bmind-046 also looks for it, and
# adopts it instead of spending 3.1 h redrawing a Gibbs sampler whose output
# the DE arm has already published results from.
fetch_bmind () {
  mkdir -p de_analysis/results
  if s3_get "${PREFIX}/data/bmind_lineage_expr.rds" \
            de_analysis/results/bmind_lineage_expr.rds; then
    echo "-- bMIND lineage profiles fetched ($(du -h de_analysis/results/bmind_lineage_expr.rds | cut -f1))"
  else
    echo "!! bmind_lineage_expr.rds missing in S3 — run network/cluster/prep_bmind.R and upload" >&2
    return 1
  fi
}

# perm-054 refits the de_final whole-blood model, which reads prep.rds: the
# filtered count matrix, covariates and lineage fractions, 74 MB and
# untracked. Only ens2sym.tsv of de_analysis/data is in git, so the sparse
# checkout gives the directory without the file and readRDS fails two seconds
# in. Same route as the reference cache and the bMIND profiles.
fetch_prep () {
  mkdir -p de_analysis/data
  if s3_get "${PREFIX}/data/prep.rds" de_analysis/data/prep.rds; then
    echo "-- prep.rds fetched ($(du -h de_analysis/data/prep.rds | cut -f1))"
  else
    echo "!! prep.rds missing in S3 - collect.py --put de_analysis/data/prep.rds data/prep.rds" >&2
    return 1
  fi
}
