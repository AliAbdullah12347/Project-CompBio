#!/usr/bin/env Rscript
# ===========================================================================
# ann-044  Module annotation on the CURRENT reference decomposition
#
# Spec: re-run the module annotation on the current network. analyses/012-full-
#       module-annotation ran GO over-representation on the RETRACTED 33-module
#       decomposition (pre-offset-fix); every annotation in it, the
#       darkturquoise=ribosomal finding included, describes modules that no
#       longer exist. WGCNA names colours by size rank, so the 71-module
#       reference re-uses the same names for different gene sets (old
#       darkturquoise 140 genes, new 153). Repeat 012 exactly against
#       runs/base-001 moduleColors: Fisher exact one-sided, FDR<0.05 within
#       module+namespace, universe = all 12368 network genes, min 10 annotated
#       genes, top 3 terms per module per namespace, plus top 10 kME hub genes
#       per module on the 234 controls. ALSO REQUIRED: write
#       reference_gene_modules.csv (gene,module) for base-001. Then report the
#       gene overlap between each old 33-module annotation and its same-named
#       new module.
#
# Two things differ from 012, both deliberate:
#
#   1. 012 built its expression matrix as t(d$logtpm[, ctrl_idx]) - the raw
#      matrix, before the per-sample offset is removed. The current reference
#      is built on expr_for(d, "control"), which removes it. Annotating the
#      current modules with the old matrix would reintroduce exactly the
#      difference this row exists to correct, so kME is computed on expr_for.
#
#   2. The spec asks for gene-level overlap against the old 33 modules. That is
#      not computable: the old labels file was overwritten when the
#      decomposition was rebuilt, and no run saved per-gene assignments - which
#      is the reason reference_gene_modules.csv is required here. What survives
#      is the annotation itself, so the comparison is made on GO terms and
#      module size, and the output says so rather than implying gene overlap
#      was measured. Writing reference_gene_modules.csv means the next
#      decomposition CAN be compared gene by gene.
#
# Dependency: base-001 reference cache.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "ann-044"
stopifnot(id == "ann-044")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

suppressPackageStartupMessages({
  library(WGCNA)
  library(AnnotationDbi)
  library(GO.db)
  library(org.Hs.eg.db)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

REF_CACHE_KEY <- "control_all_bicor_signed_p12_ds4_mch00_offsetrm"
REF_CACHE     <- file.path("network/cache", paste0("ref_", REF_CACHE_KEY, ".rds"))
OLD_GO        <- "network/analyses/012-full-module-annotation/results/top3_go_per_module.csv"
MIN_TERM      <- 10L     # annotated genes in the universe
FDR_CUT       <- 0.05
TOP_PER_NS    <- 3L
TOP_HUBS      <- 10L
SEED          <- 20261009L
set.seed(SEED)

if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY MISSING: ", REF_CACHE, " - base-001 must complete first")
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
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE)[1]),
                  error = function(e) "unknown")
  cfg <- list(id = id, family = "annot",
    spec = "re-run module annotation (GO over-representation + kME hubs) on the current base-001 decomposition; write per-gene module assignments",
    ref_cache = REF_CACHE_KEY, expr_mode = "offsetrm", kme_group = "control",
    test = "fisher exact one-sided", fdr_within = "module x namespace",
    fdr_cut = FDR_CUT, min_genes_per_term = MIN_TERM,
    top_terms_per_namespace = TOP_PER_NS, top_hubs_per_module = TOP_HUBS,
    seed = SEED,
    git_sha = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version = as.character(packageVersion("WGCNA")),
    orgdb_version = as.character(packageVersion("org.Hs.eg.db")))
  writeLines(c("{", paste(mapply(jsn, names(cfg), cfg, SIMPLIFY = TRUE),
                          collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

# ---- data ------------------------------------------------------------------
source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)
all_genes <- rownames(d$logtpm)

ref_net <- readRDS(REF_CACHE)
colors  <- ref_net$moduleColors
# The cache stores colours positionally: base-001 built them with
# labels2colors(cutreeDynamic(...)), which returns an unnamed vector, and the
# preservation rows consume them that way. expr_for() only ever subsets
# columns, so rownames(d$logtpm) is exactly the order they are in.
if (is.null(names(colors))) names(colors) <- all_genes
stopifnot(length(colors) == length(all_genes))

modules <- sort(unique(colors[colors != "grey"]))
msg(length(all_genes), " genes, ", length(modules), " non-grey modules, ",
    sum(colors == "grey"), " grey")

# ---- REQUIRED OUTPUT: per-gene module assignments --------------------------
# No run has ever saved these, which is why the old annotation cannot be
# checked gene by gene. Written first so it exists even if GO fails.
sym <- suppressMessages(AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = all_genes, column = "SYMBOL",
  keytype = "ENSEMBL", multiVals = "first"))
gene_modules <- data.frame(gene = all_genes, symbol = unname(sym),
                           module = unname(colors[all_genes]))
write.csv(gene_modules, file.path(run_dir, "reference_gene_modules.csv"),
          row.names = FALSE)
msg("wrote reference_gene_modules.csv (", nrow(gene_modules), " genes, ",
    sum(!is.na(gene_modules$symbol)), " with a symbol)")

# ---- GO universe -----------------------------------------------------------
eg2go <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = all_genes, columns = c("GOALL", "ONTOLOGYALL"),
  keytype = "ENSEMBL"))
eg2go <- eg2go[!is.na(eg2go$GOALL), ]
msg(length(unique(eg2go$ENSEMBL)), " of ", length(all_genes), " genes carry GO annotation")

go_by_ns <- lapply(c("BP", "CC", "MF"), function(ns) {
  sub_df <- eg2go[eg2go$ONTOLOGYALL == ns, ]
  terms  <- split(sub_df$ENSEMBL, sub_df$GOALL)
  terms  <- lapply(terms, function(g) intersect(unique(g), all_genes))
  terms[vapply(terms, length, 1L) >= MIN_TERM]
})
names(go_by_ns) <- c("BP", "CC", "MF")
msg("terms with >=", MIN_TERM, " genes in the universe: ",
    paste(sprintf("%s=%d", names(go_by_ns), vapply(go_by_ns, length, 1L)),
          collapse = " "))

n_all <- length(all_genes)

# Fisher's exact on a 2x2 of (in module, in term) is the hypergeometric tail,
# so phyper gives the identical one-sided p for every term at once instead of
# one fisher.test call per term per module - ~250k calls over 69 modules.
run_go <- function(module_genes, ns) {
  terms <- go_by_ns[[ns]]
  if (!length(terms)) return(NULL)
  n_mod <- length(module_genes)
  n_in  <- vapply(terms, function(g) sum(g %in% module_genes), 1L)
  n_uni <- vapply(terms, length, 1L)
  p <- phyper(n_in - 1L, n_uni, n_all - n_uni, n_mod, lower.tail = FALSE)
  res <- data.frame(go_id = names(terms), n_mod = unname(n_in),
                    n_universe = unname(n_uni), p_fisher = unname(p))
  res$fdr <- p.adjust(res$p_fisher, "BH")
  res$term_name <- suppressMessages(
    tryCatch(AnnotationDbi::Term(res$go_id), error = function(e) NA_character_))
  res[order(res$p_fisher), ]
}

# ---- kME on the 234 controls -----------------------------------------------
E_ctrl <- expr_for(d, group = "control")          # samples x genes, offset removed
msg("kME on ", nrow(E_ctrl), " controls x ", ncol(E_ctrl), " genes")
stopifnot(identical(colnames(E_ctrl), all_genes))
MEs <- WGCNA::moduleEigengenes(E_ctrl, colors = colors, excludeGrey = TRUE)$eigengenes

# ---- per module ------------------------------------------------------------
all_go <- list(); hubs <- list()
for (mod in modules) {
  mg <- all_genes[colors == mod]
  g  <- do.call(rbind, lapply(c("BP", "CC", "MF"), function(ns) {
    r <- run_go(mg, ns); if (is.null(r)) return(NULL)
    r$namespace <- ns; r
  }))
  g$module <- mod; g$n_module_genes <- length(mg)
  all_go[[mod]] <- g

  me <- MEs[[paste0("ME", mod)]]
  if (!is.null(me)) {
    k <- suppressWarnings(as.numeric(WGCNA::cor(E_ctrl[, mg, drop = FALSE], me,
                                                use = "p")))
    o <- order(k, decreasing = TRUE)[seq_len(min(TOP_HUBS, length(mg)))]
    hubs[[mod]] <- data.frame(module = mod, rank = seq_along(o),
                              gene = mg[o], symbol = unname(sym[mg[o]]),
                              kME = k[o])
  }
  msg(sprintf("  %-16s %4d genes  %4d terms FDR<%.2f", mod, length(mg),
              sum(g$fdr < FDR_CUT, na.rm = TRUE), FDR_CUT))
}
all_go <- do.call(rbind, all_go); hubs <- do.call(rbind, hubs)
rownames(all_go) <- NULL; rownames(hubs) <- NULL

write.csv(all_go, file.path(run_dir, "all_go_all_modules.csv"), row.names = FALSE)
write.csv(hubs,   file.path(run_dir, "hub_genes_top10.csv"),    row.names = FALSE)

sig <- all_go[!is.na(all_go$fdr) & all_go$fdr < FDR_CUT, ]
top <- do.call(rbind, lapply(split(sig, list(sig$module, sig$namespace), drop = TRUE),
                             function(x) head(x[order(x$p_fisher), ], TOP_PER_NS)))
rownames(top) <- NULL
write.csv(top, file.path(run_dir, "top3_go_per_module.csv"), row.names = FALSE)
msg(length(unique(sig$module)), " of ", length(modules),
    " modules have at least one term at FDR<", FDR_CUT)

# ---- what survives from 012 ------------------------------------------------
# Same-named modules are NOT the same gene sets: WGCNA names colours by size
# rank, so "darkturquoise" means whichever module ranked there in whichever
# decomposition. Comparing them by name is the question the spec asks, and the
# answer is mostly "nothing", which is the point.
if (file.exists(OLD_GO)) {
  old <- read.csv(OLD_GO)
  cmp <- do.call(rbind, lapply(sort(unique(old$module)), function(m) {
    o <- old[old$module == m, ]; n <- top[top$module == m, ]
    shared <- intersect(o$go_id, n$go_id)
    data.frame(module = m,
               old_size = if (nrow(o)) o$n_module_genes[1] else NA_integer_,
               new_size = if (m %in% modules) sum(colors == m) else 0L,
               old_top_terms = nrow(o), new_top_terms = nrow(n),
               shared_top_terms = length(shared),
               shared_term_names = paste(utils::head(
                 n$term_name[match(shared, n$go_id)], 3), collapse = "; "))
  }))
  cmp$note <- ifelse(cmp$new_size == 0, "module name absent from the 71-module reference", "")
  write.csv(cmp, file.path(run_dir, "old_vs_new_annotation.csv"), row.names = FALSE)
  msg("012 comparison: ", sum(cmp$shared_top_terms > 0), " of ", nrow(cmp),
      " old modules share a top GO term with the same-named new module")
  msg("  (term-level only - the old per-gene labels no longer exist, which is ",
      "why reference_gene_modules.csv is written here)")
} else {
  msg("012 results absent at ", OLD_GO, "; skipping the old-vs-new comparison")
}

msg("=== ", id, " COMPLETE ===")
