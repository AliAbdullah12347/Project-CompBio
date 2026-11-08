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
    "/data/analysis_matrices/*" "/data/LM22.txt" "/de_analysis/data/*" \
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
  Rscript -e 'for (p in c("WGCNA","GO.db","org.Hs.eg.db","igraph")) cat(sprintf("%-14s %s\n", p, requireNamespace(p, quietly=TRUE)))'
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
upload_run () {   # upload_run <run_id>
  local id="$1" dir="network/runs/$1"
  [ -d "$dir" ] || { echo "-- $id: no output dir"; return 0; }
  tar czf "/work/${id}.tar.gz" -C network/runs "$id"
  s3_put "/work/${id}.tar.gz" "runs/${id}.tar.gz"
  echo "-- $id: uploaded $(du -h "/work/${id}.tar.gz" | cut -f1)"
}
