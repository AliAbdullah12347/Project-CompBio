#!/usr/bin/env Rscript
# ===========================================================================
# trait-033  Module eigengene ~ clinical and batch traits
#
# Spec: eigengene ~ lithium, dx, age, sex, RIN, plate; BH across modules
#
# For each base-001 module: fit lm(eigengene ~ lithium + dx + age + sex + rin + plate)
# on all 474 samples (with group membership as covariates). Report t/F statistics
# and BH-adjusted p-values across modules for each predictor.
#
# Dependency: requires base-001 reference cache to exist.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "trait-033"
stopifnot(id == "trait-033")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

REF_CACHE_KEY <- "control_all_bicor_signed_p12"
REF_CACHE     <- file.path("network/cache", paste0("ref_", REF_CACHE_KEY, ".rds"))
SEED          <- 20261009L
FORMULA       <- ~ lithium + dx + age + sex + rin + plate

# ---- dependency check -------------------------------------------------------
if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY MISSING: ", REF_CACHE, " — base-001 must complete first")
  msg("Exiting with status=1 so driver requeues.")
  quit(save = "no", status = 1L)
}

# ---- config ----------------------------------------------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  scl <- function(v) {
    if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
    if (is.integer(v) || is.numeric(v)) return(as.character(v[1]))
    paste0('"', gsub('"', '\\"', as.character(v[1]), fixed = TRUE), '"')
  }
  jsn <- function(nm, v) paste0('  "', nm, '": ', scl(v))
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD 2>/dev/null", intern=TRUE)[1]),
                  error = function(e) "unknown")
  cfg <- list(id=id, family="modtrait",
    spec="eigengene ~ lithium + dx + age + sex + rin + plate; BH across modules",
    ref_cache=REF_CACHE_KEY,
    formula=deparse(FORMULA),
    samples="all 474",
    correction="BH across modules, per predictor",
    seed=SEED,
    git_sha=if(length(sha)&&!is.na(sha))sha else "unknown",
    r_version=paste(R.version$major,R.version$minor,sep="."),
    wgcna_version=as.character(packageVersion("WGCNA")))
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY=TRUE)
  writeLines(c("{", paste(lines, collapse=",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

source("network/00_load.R")
d <- load_project()

# ---- load reference modules ------------------------------------------------
ref_net      <- readRDS(REF_CACHE)
moduleColors <- ref_net$moduleColors
modules      <- sort(unique(moduleColors))
modules      <- modules[modules != "grey"]
msg(length(modules), " modules (excl grey)")

# ---- compute module eigengenes for ALL 474 samples -------------------------
E <- t(d$logtpm)   # 474 x 12368
MEs <- moduleEigengenes(E, colors = moduleColors, excludeGrey = TRUE)$eigengenes
rownames(MEs) <- colnames(d$logtpm)   # sample ids
msg("Module eigengenes computed: ", nrow(MEs), " samples x ", ncol(MEs), " modules")

# ---- merge with metadata ---------------------------------------------------
meta <- d$meta   # 474 rows, columns include lithium, dx, age, sex, rin, plate
stopifnot(all(meta$title == colnames(d$logtpm)))   # alignment check
df <- cbind(MEs, meta[, c("lithium","dx","age","sex","rin","plate"), drop=FALSE])

# ---- fit lm per module per predictor ---------------------------------------
predictors <- c("lithium","dx","age","sex","rin","plate")
results <- list()
for (mod in modules) {
  me_col <- paste0("ME", mod)
  if (!me_col %in% colnames(df)) me_col <- paste0("me", mod)
  if (!me_col %in% colnames(df)) {
    msg("WARNING: no eigengene column for module ", mod); next
  }
  fit <- lm(df[[me_col]] ~ lithium + dx + age + sex + rin + plate, data=df)
  cs  <- summary(fit)$coefficients
  for (pred in predictors) {
    row_nm <- pred
    if (!row_nm %in% rownames(cs)) {
      # try with 'lithium' → 'lithiumyes' etc. (factor encoding)
      matches <- rownames(cs)[grepl(paste0("^", pred), rownames(cs), ignore.case=TRUE)]
      if (!length(matches)) next
      row_nm <- matches[1]
    }
    results[[length(results)+1]] <- data.frame(
      module=mod, predictor=pred, coef=cs[row_nm,"Estimate"],
      tstat=cs[row_nm,"t value"], pval=cs[row_nm,"Pr(>|t|)"],
      stringsAsFactors=FALSE)
  }
}
res <- do.call(rbind, results)
for (pred in predictors) {
  idx <- res$predictor == pred
  if (sum(idx) > 0) res$padj[idx] <- p.adjust(res$pval[idx], method="BH")
}
res$padj[is.na(res$padj)] <- NA_real_

write.csv(res, file.path(run_dir,"results_modtrait.csv"), row.names=FALSE)
write.csv(MEs, file.path(run_dir,"module_eigengenes_474.csv"))
n_sig <- sum(res$padj < 0.05, na.rm=TRUE)
msg(n_sig, " module-predictor pairs significant at BH 5%")
msg("=== ", id, " COMPLETE ===")
