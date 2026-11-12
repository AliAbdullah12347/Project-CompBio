#!/usr/bin/env Rscript
# ===========================================================================
# ann-058  Annotate the partition the BPD analysis actually used, and map it
#          onto the fixed reference
#
# WHY THIS IS NEEDED
#
# analyses/019 scored BPD and lithium module enrichment on grp-051's
# control_A partition (69 modules). ann-044 annotated the fixed verified
# reference (71 modules). Both use WGCNA colour names, and WGCNA assigns
# colours by size rank, so `salmon` in one is a different gene set from
# `salmon` in the other. Joining 019's results to ann-044's GO table by name
# produces a table that looks right and is meaningless -- the same trap the
# eosinophil finding documents, where base-001's mediumpurple2 and grp-051's
# orangered3 are the same 33 genes under different names.
#
# So: run the same GO enrichment on grp-051's own partition, and separately
# compute the gene-overlap map between the two decompositions so a module in
# one can be found in the other.
#
# WHAT IT WRITES
#   grp051_go_all.csv         every GO term x module, grp-051 partition
#   grp051_top_go.csv         best term per module
#   grp051_bpd_annotated.csv  019's BPD/LI enrichment + the correct GO label
#   partition_map.csv         grp-051 module <-> fixed reference, by Jaccard
#   grp051_hubs.csv           top kME genes per module, on the 234 controls
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "ann-058"
stopifnot(id == "ann-058")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))
msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
for (d_ in c(run_dir, cfg_dir)) dir.create(d_, showWarnings = FALSE, recursive = TRUE)

MIN_TERM <- 10L     # as ann-044
MIN_MOD  <- 20L     # as 019
N_HUB    <- 8L

suppressPackageStartupMessages({
  library(WGCNA); library(GO.db); library(org.Hs.eg.db); library(AnnotationDbi)
})
enableWGCNAThreads(); options(stringsAsFactors = FALSE)
source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)

## ---- the two partitions -----------------------------------------------------
E_ctrl <- expr_for(d, group = "control")
all_genes <- colnames(E_ctrl)

g51 <- read.csv("network/runs/grp-051/genes_control_A.csv")
col51 <- setNames(g51$module, g51$gene)[all_genes]
stopifnot(!anyNA(col51))

gm44 <- read.csv("network/runs/ann-044/reference_gene_modules.csv")
col44 <- setNames(gm44$module, gm44$gene)[all_genes]

mods51 <- setdiff(unique(col51), "grey")
mods44 <- setdiff(unique(col44), c("grey", NA))
msg("grp-051 control_A: ", length(mods51), " modules | fixed reference: ",
    length(mods44), " modules | ", length(all_genes), " genes")

## ---- GO universe, same construction as ann-044 ------------------------------
eg2go <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys = all_genes, columns = c("GOALL", "ONTOLOGYALL"),
  keytype = "ENSEMBL"))
eg2go <- eg2go[!is.na(eg2go$GOALL), ]
go_by_ns <- lapply(c("BP", "CC", "MF"), function(ns) {
  s <- eg2go[eg2go$ONTOLOGYALL == ns, ]
  t <- lapply(split(s$ENSEMBL, s$GOALL), function(g) intersect(unique(g), all_genes))
  t[vapply(t, length, 1L) >= MIN_TERM]
})
names(go_by_ns) <- c("BP", "CC", "MF")
n_all <- length(all_genes)
msg("GO terms in universe: ",
    paste(sprintf("%s=%d", names(go_by_ns), vapply(go_by_ns, length, 1L)), collapse = " "))

run_go <- function(module_genes, ns) {
  terms <- go_by_ns[[ns]]
  if (!length(terms)) return(NULL)
  n_in  <- vapply(terms, function(g) sum(g %in% module_genes), 1L)
  n_uni <- vapply(terms, length, 1L)
  p <- phyper(n_in - 1L, n_uni, n_all - n_uni, length(module_genes), lower.tail = FALSE)
  data.frame(go_id = names(terms), namespace = ns, n_mod = unname(n_in),
             n_universe = unname(n_uni), p_fisher = unname(p))
}

## ---- enrichment on grp-051's modules ----------------------------------------
msg("GO over ", length(mods51), " grp-051 modules ...")
GO <- do.call(rbind, lapply(mods51, function(m) {
  gs <- names(col51)[col51 == m]
  r  <- do.call(rbind, lapply(c("BP", "CC", "MF"), function(ns) run_go(gs, ns)))
  if (is.null(r)) return(NULL)
  r$module <- m; r$n_module_genes <- length(gs); r
}))
GO$fdr <- ave(GO$p_fisher, GO$module, FUN = function(x) p.adjust(x, "BH"))
GO$term_name <- suppressMessages(AnnotationDbi::Term(GO$go_id))
GO <- GO[order(GO$module, GO$p_fisher), ]
write.csv(GO, file.path(run_dir, "grp051_go_all.csv"), row.names = FALSE)

top <- GO[GO$fdr < 0.05, ]
top <- top[order(top$p_fisher), ]
top <- top[!duplicated(top$module), ]
write.csv(top, file.path(run_dir, "grp051_top_go.csv"), row.names = FALSE)
msg(nrow(top), " of ", length(mods51), " grp-051 modules have a GO term at FDR<0.05")

## ---- map the two decompositions onto each other -----------------------------
# Best-match Jaccard both ways. Without this no result stated in one
# partition's colour names can be located in the other's.
msg("mapping grp-051 modules onto the fixed reference by gene overlap ...")
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
map <- do.call(rbind, lapply(mods51, function(m) {
  gs <- names(col51)[col51 == m]
  js <- vapply(mods44, function(k) jac(gs, names(col44)[col44 == k]), numeric(1))
  i  <- which.max(js)
  data.frame(grp051_module = m, n_grp051 = length(gs),
             best_fixed_module = mods44[i], n_fixed = sum(col44 == mods44[i], na.rm = TRUE),
             jaccard = unname(js[i]),
             n_shared = length(intersect(gs, names(col44)[col44 == mods44[i]])))
}))
map <- map[order(-map$jaccard), ]
write.csv(map, file.path(run_dir, "partition_map.csv"), row.names = FALSE)
msg(sprintf("  median best-match Jaccard %.3f; %d of %d modules match above 0.5",
            median(map$jaccard), sum(map$jaccard > 0.5), nrow(map)))

## ---- hubs, on the 234 controls ----------------------------------------------
msg("kME and hubs on ", nrow(E_ctrl), " controls ...")
ME <- WGCNA::moduleEigengenes(E_ctrl, colors = col51, excludeGrey = TRUE)$eigengenes
kme <- WGCNA::signedKME(E_ctrl, ME)
sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = all_genes,
         column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
hubs <- do.call(rbind, lapply(mods51, function(m) {
  cn <- paste0("kME", m); if (!cn %in% colnames(kme)) return(NULL)
  gs <- names(col51)[col51 == m]
  # signedKME returns a data.frame, and df[rows, onecol] drops the names, so
  # names(v) would be NULL and every hub row would come back empty.
  v  <- sort(setNames(kme[gs, cn], gs), decreasing = TRUE)[seq_len(min(N_HUB, length(gs)))]
  data.frame(module = m, rank = seq_along(v), gene = names(v),
             symbol = unname(ifelse(is.na(sym[names(v)]), names(v), sym[names(v)])),
             kME = unname(v))
}))
write.csv(hubs, file.path(run_dir, "grp051_hubs.csv"), row.names = FALSE)

## ---- join 019's BPD/LI enrichment to the CORRECT annotation -----------------
E019 <- read.csv("network/analyses/019-bpd-module-connection/results/module_rank_enrichment.csv")
E019$go_term <- top$term_name[match(E019$module, top$module)]
E019$go_fdr  <- top$fdr[match(E019$module, top$module)]
E019$fixed_counterpart <- map$best_fixed_module[match(E019$module, map$grp051_module)]
E019$counterpart_jaccard <- map$jaccard[match(E019$module, map$grp051_module)]
E019$hub_genes <- vapply(E019$module, function(m)
  paste(head(hubs$symbol[hubs$module == m], 5), collapse = ", "), character(1))
write.csv(E019, file.path(run_dir, "grp051_bpd_annotated.csv"), row.names = FALSE)

bpd <- E019[E019$deg_list == "WB_BPD" & E019$padj < 0.05, ]
li  <- E019[E019$deg_list == "WB_LI"  & E019$padj < 0.05, ]
msg("")
msg("BPD-enriched: ", nrow(bpd), " | lithium-enriched: ", nrow(li),
    " | both: ", length(intersect(bpd$module, li$module)),
    " | BPD only: ", length(setdiff(bpd$module, li$module)),
    " | lithium only: ", length(setdiff(li$module, bpd$module)))
msg("")
msg("top BPD-enriched modules, correctly annotated:")
b <- bpd[order(bpd$p), ][seq_len(min(10, nrow(bpd))), ]
for (i in seq_len(nrow(b))) msg(sprintf("  %-15s n=%4d padj=%.1e  %-34s | hubs: %s",
  b$module[i], b$n[i], b$padj[i],
  if (is.na(b$go_term[i])) "(no GO at FDR<0.05)" else substr(b$go_term[i], 1, 34),
  substr(b$hub_genes[i], 1, 38)))

msg("=== ", id, " COMPLETE ===")
