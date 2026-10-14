#!/usr/bin/env Rscript
# ===========================================================================
# 012-full-module-annotation
#
# Analyses 006 and 011 annotated two modules (darkturquoise, darkred) in
# detail. This analysis runs GO over-representation for all 33 non-grey
# modules, producing a summary table with top 3 GO terms per module in
# each namespace. Used to interpret base-001 preservation results when they
# arrive.
#
# Pre-specified: Fisher exact, one-sided, FDR<0.05 within module+namespace,
# min 10 annotated genes in universe, universe = all 12,368 network genes.
# Top 3 terms by p-value per module per namespace.
# kME hub genes (top 5 per module) computed on control samples (n=234).
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(WGCNA)
  library(AnnotationDbi)
  library(GO.db)
})
enableWGCNAThreads()
OUT <- "network/analyses/012-full-module-annotation/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)

msg("Loading project data ...")
source("network/00_load.R"); d <- load_project(cache = TRUE, verbose = FALSE)
all_genes <- rownames(d$logtpm)
msg(length(all_genes), " genes in network")

msg("Loading module labels ...")
lab    <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")
colors <- WGCNA::labels2colors(lab)
names(colors) <- all_genes
modules <- sort(unique(colors[colors != "grey"]))
msg(length(modules), " non-grey modules")

# Load group ANOVA results from analysis 009
anova09 <- read.csv("network/analyses/009-all-modules-across-groups/results/all_modules_group_anova.csv",
                    stringsAsFactors = FALSE)

# ---- GO enrichment setup ---------------------------------------------------
org_ok <- suppressPackageStartupMessages(
  tryCatch({ library(org.Hs.eg.db); TRUE }, error = function(e) FALSE)
)
msg("org.Hs.eg.db available: ", org_ok)
if (!org_ok) stop("org.Hs.eg.db required; install via BiocManager")

eg2go <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = all_genes,
                        columns = c("GOALL", "ONTOLOGYALL"),
                        keytype = "ENSEMBL")
)
eg2go <- eg2go[!is.na(eg2go$GOALL), ]
msg(length(unique(eg2go$ENSEMBL)), " genes with GO annotation")

# Pre-filter GO term lists (all namespaces) to save time in the inner loop
go_by_ns <- lapply(c("BP","CC","MF"), function(ns) {
  sub_df <- eg2go[eg2go$ONTOLOGYALL == ns, ]
  terms  <- split(sub_df$ENSEMBL, sub_df$GOALL)
  # Keep only terms with >=10 genes in the universe
  terms[sapply(terms, function(g) length(intersect(g, all_genes)) >= 10)]
})
names(go_by_ns) <- c("BP","CC","MF")

run_go <- function(module_genes, ns) {
  go_terms <- go_by_ns[[ns]]
  if (length(go_terms) == 0) return(NULL)
  n_all <- length(all_genes)
  n_mod <- length(module_genes)
  res <- do.call(rbind, lapply(names(go_terms), function(term) {
    in_term <- go_terms[[term]]
    a <- length(intersect(module_genes, in_term))
    b <- n_mod - a
    c <- length(intersect(all_genes, in_term)) - a
    dd <- n_all - n_mod - c
    pval <- fisher.test(matrix(c(a, b, c, dd), 2, 2), alternative = "greater")$p.value
    data.frame(go_id = term, n_mod = a,
               n_universe = length(intersect(all_genes, in_term)),
               p_fisher = pval, stringsAsFactors = FALSE)
  }))
  res$fdr <- p.adjust(res$p_fisher, "BH")
  res$term_name <- tryCatch(Term(res$go_id), error = function(e) NA_character_)
  res[order(res$p_fisher), ]
}

# ---- kME in control samples -----------------------------------------------
msg("Computing module eigengenes in control samples ...")
ctrl_idx <- !is.na(d$groups) & d$groups == "control"
E_ctrl <- t(d$logtpm[, ctrl_idx])
MEs_ctrl <- WGCNA::moduleEigengenes(E_ctrl, colors = colors,
                                     excludeGrey = TRUE)$eigengenes

# ---- Per-module annotation -------------------------------------------------
all_go_rows  <- list()
top_go_rows  <- list()
hub_rows     <- list()

for (mod in modules) {
  mod_genes <- all_genes[colors == mod]
  n_mod <- length(mod_genes)
  msg(mod, " (", n_mod, " genes) ...")

  # kME
  me_col <- paste0("ME", mod)
  if (me_col %in% names(MEs_ctrl)) {
    me_vec <- MEs_ctrl[[me_col]]
    E_mod  <- E_ctrl[, mod_genes, drop = FALSE]
    kme    <- apply(E_mod, 2, function(g) cor(g, me_vec, use = "pairwise.complete.obs"))
    kme_df <- data.frame(module = mod, ensembl_id = mod_genes,
                         kME = round(kme, 4), stringsAsFactors = FALSE)
    kme_df <- kme_df[order(-kme_df$kME), ]
    hub_rows[[mod]] <- head(kme_df, 5)
  }

  # GO per namespace
  for (ns in c("BP","CC","MF")) {
    r <- run_go(mod_genes, ns)
    if (is.null(r) || nrow(r) == 0) next
    r$module <- mod; r$namespace <- ns; r$n_module_genes <- n_mod
    all_go_rows[[paste(mod, ns)]] <- r

    sig <- r[!is.na(r$fdr) & r$fdr < 0.05, ]
    top <- head(r, 3)  # top 3 by p regardless of FDR
    top$n_sig_fdr05 <- nrow(sig)
    top$module <- mod; top$namespace <- ns; top$n_module_genes <- n_mod
    top_go_rows[[paste(mod, ns)]] <- top
  }
}

# ---- Save full results -----------------------------------------------------
msg("Saving full GO results ...")
all_go <- do.call(rbind, all_go_rows)
write.csv(all_go, file.path(OUT, "all_go_all_modules.csv"), row.names = FALSE)

top_go <- do.call(rbind, top_go_rows)
write.csv(top_go, file.path(OUT, "top3_go_per_module.csv"), row.names = FALSE)

hubs <- do.call(rbind, hub_rows)
write.csv(hubs, file.path(OUT, "hub_genes_top5.csv"), row.names = FALSE)

# ---- Summary table per module (one row per module) -------------------------
msg("Building per-module summary table ...")
summ <- do.call(rbind, lapply(modules, function(mod) {
  a09 <- anova09[anova09$module == mod, ]
  row_base <- data.frame(
    module = mod,
    n_genes = sum(colors == mod),
    ctrl_nolith_bonf = if (nrow(a09) > 0) round(a09$p_bonf_ctrl_nolith[1], 4) else NA,
    ctrl_lith_bonf   = if (nrow(a09) > 0) round(a09$p_bonf_ctrl_lith[1], 4) else NA,
    nolith_lith_bonf = if (nrow(a09) > 0) round(a09$p_bonf_nolith_lith[1], 4) else NA,
    r2_comp          = if (nrow(a09) > 0) round(a09$r2_comp[1], 4) else NA,
    stringsAsFactors = FALSE
  )
  for (ns in c("BP","CC","MF")) {
    key <- paste(mod, ns)
    top3 <- tryCatch(top_go_rows[[key]], error = function(e) NULL)
    if (is.null(top3) || nrow(top3) == 0) {
      row_base[[paste0("top_", ns)]] <- "none"
      row_base[[paste0("top_", ns, "_fdr")]] <- NA
    } else {
      row_base[[paste0("top_", ns)]] <- paste(top3$term_name[1], collapse = "")
      row_base[[paste0("top_", ns, "_fdr")]] <- round(top3$fdr[1], 4)
    }
  }
  row_base
}))

# Merge hub gene 1 symbol
# Resolve top hub gene Ensembl -> symbol for each module
hub_top <- hubs[!duplicated(hubs$module), c("module","ensembl_id","kME")]
syms <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
  keys = hub_top$ensembl_id, columns = "SYMBOL", keytype = "ENSEMBL"))
syms <- syms[!duplicated(syms$ENSEMBL), ]
hub_top <- merge(hub_top, syms, by.x = "ensembl_id", by.y = "ENSEMBL", all.x = TRUE)
hub_top$hub <- paste0(hub_top$SYMBOL, " (", hub_top$kME, ")")
summ <- merge(summ, hub_top[, c("module","hub")], by = "module", all.x = TRUE)
summ <- summ[order(summ$ctrl_nolith_bonf, na.last = TRUE), ]

write.csv(summ, file.path(OUT, "module_summary_table.csv"), row.names = FALSE)

cat("\n======= Module summary (top 20 by ctrl-nolith effect) =======\n")
print(summ[1:min(20, nrow(summ)),
           c("module","n_genes","ctrl_nolith_bonf","ctrl_lith_bonf","hub","top_BP")],
      row.names = FALSE)

cat("\nANALYSIS_012_DONE\n")
