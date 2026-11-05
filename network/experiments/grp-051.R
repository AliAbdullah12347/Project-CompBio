#!/usr/bin/env Rscript
# ===========================================================================
# grp-051 — group-specific networks, hub genes, module overlap, and the
#           modules' relationship to diagnosis. Built for the research
#           cluster (24 cores); runs anywhere R + WGCNA does.
#
# WHAT THIS ANSWERS, that nothing in the arm so far could.
#   Every earlier run builds ONE reference network and measures how well those
#   modules survive elsewhere. So there has never been a "lithium users'
#   network" to set beside a "non-users' network", and no run saved which gene
#   sits in which module. This rebuilds the network separately per group and
#   writes the gene-level tables out.
#
# TWO TRAPS IT IS BUILT AROUND.
#   1. Sample size. The existing references are n = 234 / 74 / 152, and module
#      count moves with n for reasons unrelated to biology. Every network here
#      is n = 74 with identical parameters.
#   2. Instability. Two networks built on different people from the SAME
#      population also disagree. So each draw builds two DISJOINT control
#      networks, and their agreement is the ceiling every cross-group number
#      is read against. Without it "lithium and control share 40% of modules"
#      is uninterpretable.
#
#   bp_nolith is exactly 74 people, so it cannot be subsampled: it is built
#   once and reused. Its comparisons therefore carry no draw-to-draw spread,
#   which is a real limit of this cohort, not of the method.
#
# Parameters follow ref-003.R so results sit on the same footing as the arm:
# bicor, signed, power 12, deepSplit 4, minModuleSize 30, mergeCutHeight 0,
# offset-removed expression.
# ===========================================================================
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

NTHREAD <- as.integer(Sys.getenv("GRP_THREADS", "4"))
N_DRAWS <- as.integer(Sys.getenv("GRP_DRAWS", "8"))
OUT     <- Sys.getenv("GRP_OUT", "network/runs/grp-051")
enableWGCNAThreads(NTHREAD)

source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)

CK <- file.path(OUT, "ck"); dir.create(CK, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                   paste0(..., collapse = ""))); flush.console() }
wr <- function(x, f) write.csv(x, file.path(OUT, f), row.names = FALSE)

POWER <- 12L; NET <- "signed"; DEEP <- 4L; MINMOD <- 30L; MERGE <- 0.0
N <- 74L; SEED <- 20261105L
sym <- d$sym
symof <- function(g) ifelse(g %in% names(sym) & nzchar(sym[g]), sym[g], g)

ids   <- colnames(d$logtpm)
ctrl  <- ids[which(d$groups == "control")]
nolit <- ids[which(d$groups == "bp_nolith")]
lith  <- ids[which(d$groups == "bp_lith")]
stopifnot(length(nolit) == N, length(ctrl) == 234, length(lith) == 152)
msg("groups: control ", length(ctrl), ", bp_nolith ", length(nolit), ", bp_lith ", length(lith))

# -- sample sets: per draw, two disjoint control halves plus a lithium draw --
set.seed(SEED)
SETS <- list(bp_nolith = nolit)
for (i in seq_len(N_DRAWS)) {
  cs <- sample(ctrl, 2 * N)
  SETS[[sprintf("control_A_%02d", i)]] <- cs[1:N]
  SETS[[sprintf("control_B_%02d", i)]] <- cs[(N + 1):(2 * N)]
  SETS[[sprintf("bp_lith_%02d",   i)]] <- sample(lith, N)
}
stopifnot(all(vapply(SETS, length, 1L) == N))

# -- build one network ------------------------------------------------------
build <- function(nm) {
  f <- file.path(CK, paste0(nm, ".rds"))
  if (file.exists(f)) { msg(nm, ": checkpoint hit"); return(readRDS(f)) }
  E <- expr_for(d, group = NULL)[SETS[[nm]], , drop = FALSE]
  stopifnot(nrow(E) == N, ncol(E) == nrow(d$logtpm))
  t0 <- Sys.time()
  bwm <- blockwiseModules(E, power = POWER, networkType = NET, TOMType = NET,
                          corType = "bicor", deepSplit = DEEP,
                          minModuleSize = MINMOD, mergeCutHeight = MERGE,
                          maxBlockSize = 15000L, numericLabels = FALSE,
                          saveTOMs = FALSE, verbose = 0)
  # kME = each gene's correlation with its own module's eigengene. High kME
  # means central to the module; the top few are its hub genes.
  kme <- signedKME(E, bwm$MEs, corFnc = "bicor")
  own <- paste0("kME", bwm$colors)
  k <- vapply(seq_along(bwm$colors), function(i)
    if (own[i] %in% colnames(kme)) kme[i, own[i]] else NA_real_, numeric(1))
  res <- list(colors = bwm$colors, MEs = bwm$MEs, genes = colnames(E), kME = k,
              secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  saveRDS(res, f); rm(E, bwm, kme); gc(verbose = FALSE)
  msg(nm, ": ", length(unique(res$colors)) - ("grey" %in% res$colors),
      " modules, ", sum(res$colors != "grey"), " genes assigned, ", round(res$secs), "s")
  res
}

NM <- names(SETS)
msg("building ", length(NM), " networks at n = ", N, " on ", NTHREAD, " threads")
NETS <- list()
for (nm in NM) NETS[[nm]] <- build(nm)

# ===========================================================================
# 1. per-network gene tables, module sizes, hub genes
# ===========================================================================
gene_table <- function(nm) {
  n <- NETS[[nm]]
  data.frame(network = nm, gene = n$genes, symbol = symof(n$genes),
             module = n$colors, kME = round(n$kME, 4))
}
# full gene-level tables for the four representative networks only; writing
# all 25 would be ~300 MB of near-duplicate rows.
REP <- c(control_A = "control_A_01", control_B = "control_B_01",
         bp_nolith = "bp_nolith",    bp_lith   = "bp_lith_01")
for (lab in names(REP)) {
  g <- gene_table(REP[[lab]])
  wr(g, sprintf("genes_%s.csv", lab))
  sz <- as.data.frame(table(g$module[g$module != "grey"]))
  names(sz) <- c("module", "n_genes")
  wr(sz[order(-sz$n_genes), ], sprintf("module_sizes_%s.csv", lab))
  gg <- g[g$module != "grey", ]
  hub <- do.call(rbind, lapply(split(gg, gg$module), function(x) head(x[order(-x$kME), ], 15)))
  wr(hub, sprintf("hubs_%s.csv", lab))
}

# module counts across every network, so the per-group spread is visible
cnt <- do.call(rbind, lapply(NM, function(nm) data.frame(
  network = nm,
  group = sub("_[0-9]+$", "", nm),
  n_modules = length(unique(NETS[[nm]]$colors)) - ("grey" %in% NETS[[nm]]$colors),
  n_assigned = sum(NETS[[nm]]$colors != "grey"),
  n_grey = sum(NETS[[nm]]$colors == "grey"),
  largest_module_pct = round(100 * max(table(NETS[[nm]]$colors[NETS[[nm]]$colors != "grey"])) /
                               length(NETS[[nm]]$colors), 2),
  secs = round(NETS[[nm]]$secs))))
wr(cnt, "module_counts.csv")

# ===========================================================================
# 2. partition agreement, every comparison read against the ceiling
# ===========================================================================
ari <- function(a, b) {
  tb <- table(a, b); n <- sum(tb); s <- function(x) sum(choose(x, 2))
  i <- s(tb); ea <- s(rowSums(tb)); eb <- s(colSums(tb)); tot <- choose(n, 2)
  e <- ea * eb / tot
  (i - e) / ((ea + eb) / 2 - e)
}
pair_tab <- function(A, B, write_as = NULL) {
  a <- NETS[[A]]; b <- NETS[[B]]; stopifnot(identical(a$genes, b$genes))
  keep <- a$colors != "grey" & b$colors != "grey"
  rand <- if (sum(keep) > 10) ari(a$colors[keep], b$colors[keep]) else NA_real_
  tb <- table(a$colors, b$colors)
  tb <- tb[rownames(tb) != "grey", colnames(tb) != "grey", drop = FALSE]
  sa <- table(a$colors)[rownames(tb)]; sb <- table(b$colors)[colnames(tb)]
  df <- do.call(rbind, lapply(rownames(tb), function(m) {
    ov <- tb[m, ]; j <- ov / (as.integer(sa[m]) + as.integer(sb) - ov); k <- which.max(j)
    data.frame(module_A = m, n_A = as.integer(sa[m]), best_match_B = colnames(tb)[k],
               n_B = as.integer(sb[k]), n_shared = as.integer(ov[k]),
               jaccard = round(max(j), 3))
  }))
  if (!is.null(write_as)) wr(df[order(-df$jaccard), ], write_as)
  data.frame(A = A, B = B, n_modules_A = nrow(tb), n_modules_B = ncol(tb),
             ARI = round(rand, 4), median_jaccard = round(median(df$jaccard), 3),
             frac_j50 = round(mean(df$jaccard >= 0.5), 3),
             frac_j25 = round(mean(df$jaccard >= 0.25), 3))
}

COMPS <- list()
for (i in seq_len(N_DRAWS)) {
  a <- sprintf("control_A_%02d", i); b <- sprintf("control_B_%02d", i)
  l <- sprintf("bp_lith_%02d", i)
  COMPS[[length(COMPS) + 1]] <- cbind(contrast = "ceiling: control vs control",   pair_tab(a, b))
  COMPS[[length(COMPS) + 1]] <- cbind(contrast = "bipolar: control vs bp_nolith", pair_tab(a, "bp_nolith"))
  COMPS[[length(COMPS) + 1]] <- cbind(contrast = "lithium: bp_nolith vs bp_lith", pair_tab("bp_nolith", l))
  COMPS[[length(COMPS) + 1]] <- cbind(contrast = "control vs bp_lith",            pair_tab(a, l))
}
CMP <- do.call(rbind, COMPS); wr(CMP, "pair_all_draws.csv")
# one worked example per contrast, with the per-module matches
pair_tab(REP[["control_A"]], REP[["control_B"]], "pair_example_ceiling.csv")
pair_tab(REP[["control_A"]], "bp_nolith",        "pair_example_bipolar.csv")
pair_tab("bp_nolith", REP[["bp_lith"]],          "pair_example_lithium.csv")

agg <- do.call(rbind, lapply(split(CMP, CMP$contrast), function(x) data.frame(
  contrast = x$contrast[1], n_comparisons = nrow(x),
  ARI_median = round(median(x$ARI, na.rm = TRUE), 4),
  ARI_min = round(min(x$ARI, na.rm = TRUE), 4), ARI_max = round(max(x$ARI, na.rm = TRUE), 4),
  median_jaccard = round(median(x$median_jaccard), 3),
  frac_j50 = round(median(x$frac_j50), 3))))
agg <- agg[order(-agg$ARI_median), ]; wr(agg, "agreement_summary.csv")

# ===========================================================================
# 3. hub genes: are the same genes central in each group?
# ===========================================================================
topk <- function(nm, k = 200) {
  n <- NETS[[nm]]; ok <- n$colors != "grey" & !is.na(n$kME)
  n$genes[ok][order(-n$kME[ok])][seq_len(min(k, sum(ok)))]
}
hub_rows <- list()
for (i in seq_len(N_DRAWS)) {
  a <- sprintf("control_A_%02d", i); b <- sprintf("control_B_%02d", i); l <- sprintf("bp_lith_%02d", i)
  f <- function(lab, x, y) data.frame(contrast = lab,
        shared = length(intersect(topk(x), topk(y))), pct = round(100 * length(intersect(topk(x), topk(y))) / 200, 1))
  hub_rows[[length(hub_rows) + 1]] <- f("ceiling: control vs control", a, b)
  hub_rows[[length(hub_rows) + 1]] <- f("bipolar: control vs bp_nolith", a, "bp_nolith")
  hub_rows[[length(hub_rows) + 1]] <- f("lithium: bp_nolith vs bp_lith", "bp_nolith", l)
}
HB <- do.call(rbind, hub_rows)
wr(do.call(rbind, lapply(split(HB, HB$contrast), function(x) data.frame(
  contrast = x$contrast[1], median_shared_of_200 = median(x$shared),
  median_pct = median(x$pct)))), "hub_overlap_summary.csv")

# genes that are hubs in every draw of a group — the stable centres
stable_hub <- function(prefix, k = 200) {
  nets <- grep(prefix, NM, value = TRUE)
  tb <- table(unlist(lapply(nets, topk, k = k)))
  data.frame(gene = names(tb), symbol = symof(names(tb)),
             n_draws = as.integer(tb), n_possible = length(nets))
}
SH <- rbind(cbind(group = "control",   stable_hub("^control_A_")),
            cbind(group = "bp_lith",   stable_hub("^bp_lith_")),
            cbind(group = "bp_nolith", stable_hub("^bp_nolith$")))
SH <- SH[order(SH$group, -SH$n_draws), ]; wr(SH, "stable_hub_genes.csv")

# ===========================================================================
# 4. what the modules have to do with diagnosis and with lithium
# ===========================================================================
# Eigengenes of the control reference projected onto ALL 474 people, then
# regressed on diagnosis and on lithium. This is the module-level version of
# the arm's central question, computed on a network built per group.
ref <- NETS[[REP[["control_A"]]]]
Eall <- expr_for(d, group = NULL)
ME <- moduleEigengenes(Eall, colors = ref$colors, excludeGrey = TRUE)$eigengenes
md <- d$meta
md$dx_bp  <- as.integer(d$groups %in% c("bp_nolith", "bp_lith"))
md$li     <- as.integer(d$groups == "bp_lith")
md$nolith <- as.integer(d$groups == "bp_nolith")
ilrv <- d$ilr[, c("b1", "b2", "b3", "b4")]

fit_one <- function(y, X, term) {
  f <- summary(lm(y ~ ., data = cbind(y = y, X)))$coefficients
  if (!(term %in% rownames(f))) return(c(NA, NA, NA))
  f[term, c(1, 3, 4)]
}
rows <- list()
for (mn in colnames(ME)) {
  y <- ME[[mn]]
  a <- fit_one(y, data.frame(li = md$li, age = md$age, sex = md$sex, rin = md$rin), "li")
  b <- fit_one(y, data.frame(dx = md$dx_bp, age = md$age, sex = md$sex, rin = md$rin), "dx")
  cc <- fit_one(y, cbind(data.frame(li = md$li, age = md$age, sex = md$sex, rin = md$rin), ilrv), "li")
  dd <- fit_one(y, data.frame(b1 = ilrv$b1, b2 = ilrv$b2, b3 = ilrv$b3, b4 = ilrv$b4), "b1")
  rows[[mn]] <- data.frame(module = sub("^ME", "", mn),
    li_beta = a[1], li_t = a[2], li_p = a[3],
    dx_beta = b[1], dx_t = b[2], dx_p = b[3],
    li_adj_beta = cc[1], li_adj_p = cc[3],
    b1_beta = dd[1], b1_p = dd[3])
}
MT <- do.call(rbind, rows)
for (p in c("li_p", "dx_p", "li_adj_p", "b1_p"))
  MT[[sub("_p$", "_padj", p)]] <- p.adjust(MT[[p]], "BH")
MT <- MT[order(MT$li_p), ]; wr(MT, "module_trait.csv")

wr(data.frame(
  n_modules        = nrow(MT),
  lithium_BH05     = sum(MT$li_padj  < 0.05, na.rm = TRUE),
  diagnosis_BH05   = sum(MT$dx_padj  < 0.05, na.rm = TRUE),
  composition_BH05 = sum(MT$b1_padj  < 0.05, na.rm = TRUE),
  lithium_BH05_after_composition_adj = sum(MT$li_adj_padj < 0.05, na.rm = TRUE)
), "module_trait_summary.csv")

# overlap of each module with the arm's DE gene lists
deg_files <- c(WB_LI = "de_final/results/DE_WB_LI.csv.gz",
               WB_BPD = "de_final/results/DE_WB_BPD.csv.gz")
ov <- list()
for (lab in names(deg_files)) {
  if (!file.exists(deg_files[[lab]])) next
  t <- read.csv(gzfile(deg_files[[lab]]))
  pc <- intersect(c("adj.P.Val", "padj", "FDR", "BH"), names(t))
  gc_ <- intersect(c("gene", "Gene", "ensembl", "ID"), names(t))
  if (!length(pc) || !length(gc_)) next
  degs <- intersect(t[[gc_[1]]][t[[pc[1]]] < 0.05], ref$genes)
  if (length(degs) < 5) next
  for (m in setdiff(unique(ref$colors), "grey")) {
    gin <- ref$genes[ref$colors == m]
    tt <- matrix(c(length(intersect(gin, degs)), length(setdiff(gin, degs)),
                   length(setdiff(degs, gin)),
                   length(ref$genes) - length(union(gin, degs))), 2, 2)
    ft <- fisher.test(tt, alternative = "greater")
    ov[[length(ov) + 1]] <- data.frame(deg_list = lab, module = m, n_module = length(gin),
      n_deg = length(degs), n_overlap = tt[1, 1],
      expected = round(length(gin) * length(degs) / length(ref$genes), 2),
      odds_ratio = round(unname(ft$estimate), 3), p = ft$p.value)
  }
}
if (length(ov)) {
  OV <- do.call(rbind, ov)
  OV$padj <- ave(OV$p, OV$deg_list, FUN = function(x) p.adjust(x, "BH"))
  wr(OV[order(OV$deg_list, OV$p), ], "module_deg_overlap.csv")
}

# ===========================================================================
# 5. GO enrichment — what the modules actually are
# ===========================================================================
if (requireNamespace("GO.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(GO.db); library(org.Hs.eg.db) })
  msg("GO enrichment on the control reference ...")
  uni <- ref$genes
  e <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = uni,
         keytype = "ENSEMBL", columns = c("ENTREZID")))
  e <- e[!is.na(e$ENTREZID) & !duplicated(e$ENSEMBL), ]
  g2e <- setNames(e$ENTREZID, e$ENSEMBL)
  go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(e$ENTREZID),
         keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
  go <- go[!is.na(go$GOALL), ]
  uni_e <- unique(e$ENTREZID)
  go <- go[go$ENTREZID %in% uni_e, ]
  bysets <- split(go$ENTREZID, paste(go$GOALL, go$ONTOLOGYALL, sep = "|"))
  bysets <- bysets[vapply(bysets, function(x) length(unique(x)) >= 10, TRUE)]
  msg("  ", length(bysets), " GO terms with >= 10 genes in the universe")
  res <- list()
  for (m in setdiff(unique(ref$colors), "grey")) {
    gin <- unique(g2e[ref$genes[ref$colors == m]]); gin <- gin[!is.na(gin)]
    if (length(gin) < 10) next
    pv <- vapply(bysets, function(s) {
      s <- unique(s); a <- length(intersect(gin, s))
      if (a < 3) return(1)
      fisher.test(matrix(c(a, length(gin) - a, length(s) - a,
                           length(uni_e) - length(union(gin, s))), 2, 2),
                  alternative = "greater")$p.value
    }, numeric(1))
    k <- order(pv)[1:min(5, length(pv))]
    nm2 <- names(bysets)[k]
    gid <- sub("\\|.*$", "", nm2); ont <- sub("^.*\\|", "", nm2)
    trm <- suppressMessages(AnnotationDbi::select(GO.db, keys = gid,
             keytype = "GOID", columns = "TERM"))
    res[[length(res) + 1]] <- data.frame(module = m, n_module = length(gin),
      go_id = gid, ontology = ont,
      term = trm$TERM[match(gid, trm$GOID)],
      p = pv[k], fdr = p.adjust(pv, "BH")[k])
  }
  if (length(res)) wr(do.call(rbind, res), "module_go_top5.csv")
  msg("GO done")
} else {
  msg("GO packages absent — skipping enrichment")
}

msg("=== module counts by group ===")
print(do.call(rbind, lapply(split(cnt, cnt$group), function(x) data.frame(
  group = x$group[1], n_networks = nrow(x),
  median_modules = median(x$n_modules), min = min(x$n_modules), max = max(x$n_modules),
  median_assigned = median(x$n_assigned)))), row.names = FALSE)
msg("=== agreement, every contrast against the ceiling ===")
print(agg, row.names = FALSE)
msg("=== hub overlap ==="); print(HB[!duplicated(HB$contrast), ], row.names = FALSE)
msg("done")
