#!/usr/bin/env Rscript
# Co-expression network analysis, reproducing Krebs et al. 2020 Psychological Medicine.
#
# Supplementary Methods, "Co-expression network analysis":
#   covariates = age + sex + tobacco + assessment group + RIN + plate + seqPC1-3
#   beta = 7, minModuleSize = 30, deepSplit = 2, mergeCutHeight = 0.25
#
# Everything else about the network construction is unstated, and the unstated part
# is what decides the answer. An earlier version of this script used
# blockwiseModules() with its defaults and produced 25 modules with 2,960 genes
# unassigned, against their 27 modules with 314 unassigned. src/03b_wgcna_tuning.R
# sweeps the undocumented parameters; this script is set to the winner. Three pieces
# of evidence, two of them from their own deposited file, say the authors used the
# step-by-step WGCNA tutorial workflow rather than blockwiseModules:
#
#   1. Their 27 module colours sit at scattered positions 2-52 of standardColors(),
#      with turquoise, yellow and green missing while red and brown are the two
#      largest modules. blockwiseModules always calls mergeCloseModules(relabel=TRUE),
#      which renames survivors turquoise, blue, brown ... in descending size, so 27
#      modules would have to carry standardColors()[1:27] with no gaps. The default,
#      relabel = FALSE, instead gives a merged branch the alphabetically first colour
#      among its members -- and every one of the 25 absorbed colours is alphabetically
#      after some surviving colour, with no exceptions.
#   2. The clincher is bisque4. It is standardColors()[47], so before merging it was
#      only the 47th largest module, yet it is their 6th largest at 622 genes. It is
#      also the alphabetically first colour in standardColors()[1:52], so under
#      relabel = FALSE it can never be absorbed and can absorb anything. Under
#      relabel = TRUE it could not exist.
#   3. blockwiseModules runs a kME quality-control pass that the tutorial does not
#      (minCoreKME = 0.5, minKMEtoStay = 0.3): it deletes weakly-cohesive modules and
#      moves individual low-kME genes to grey. It is the only step in the pipeline
#      that can take an already-clustered gene and unassign it, and switching it off
#      is what closes most of the 2,960-vs-314 gap. See results/wgcna_param_sweep.csv.
#
# The module-trait statistic is also mis-stated in their methods, which describe it as
# the correlation of module membership with gene significance (a correlation ACROSS
# GENES). That definition does not reproduce their published p-values. Back-solving
# r = 0.156 -> p = 9.4e-4 gives n = 444, i.e. the published numbers are eigengene-vs-
# trait correlations ACROSS SAMPLES. We compute both and treat the eigengene version
# as primary.
#
# The network is the expensive step (12,344 genes in one dense TOM on a machine with
# little free RAM), so the residual matrix and the fitted network are cached to
# data/processed/. Delete those .rds files to force a recompute.

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(WGCNA); library(dynamicTreeCut)
  library(fastcluster)
})
options(stringsAsFactors = FALSE)

# Resolve the project root from the script's own location so the pipeline can
# be run from any working directory.
script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- if (!is.na(script_path)) dirname(dirname(normalizePath(script_path))) else "."
if (!dir.exists(file.path(ROOT, "data"))) ROOT <- "."
PROC <- file.path(ROOT, "data", "processed")
REF  <- file.path(ROOT, "data", "reference")
RES  <- file.path(ROOT, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------- network configuration
# Fixed by the paper's supplementary methods.
BETA <- 7; MIN_MODULE_SIZE <- 30; DEEP_SPLIT <- 2; MERGE_CUT_HEIGHT <- 0.25
# Chosen by the sweep in src/03b_wgcna_tuning.R; see results/wgcna_param_sweep.csv.
NETWORK_TYPE        <- "unsigned"
TOM_TYPE            <- "signed"
COR_TYPE            <- "pearson"
DETECT_CUT_HEIGHT   <- NA      # NA = dynamicTreeCut's own default, as the tutorial
PAM_STAGE           <- TRUE
PAM_RESPECTS_DENDRO <- FALSE
KME_FILTER          <- FALSE   # blockwiseModules' kME purge; the tutorial has none
# mergeCloseModules' useAbs decides whether two modules with strongly ANTI-correlated
# eigengenes count as similar. FALSE is the default of both WGCNA workflows, so it is
# the assumption-free choice, but in an unsigned network -- where the adjacency has
# already discarded the sign -- TRUE is arguably the consistent one, and the methods
# do not say. The sweep scored both on this exact cut and they disagree about which
# is closer to the published network:
#   useAbs = FALSE -> 31 modules, grey 301, sizes 47-1697, mean 388.5, ARI 0.4482
#   useAbs = TRUE  -> 25 modules, grey 301, sizes 47-2107, mean 481.7, ARI 0.3766
# TRUE is nearer their 27 modules and their 2,760 largest; FALSE agrees better with
# their actual gene-by-gene assignment, which is the only one of those numbers that
# measures shared structure rather than a marginal summary. We keep FALSE, on that
# ground and because it needs no undocumented deviation, and report the trade-off
# rather than burying it. Both rows are in results/wgcna_param_sweep.csv.
MERGE_USE_ABS       <- FALSE

CACHE_RESID <- file.path(PROC, "wgcna_residuals.rds")
CACHE_NET   <- file.path(PROC, "wgcna_network_stepwise.rds")

# ---------------------------------------------------------------- metadata
cat("Loading metadata ...\n")
meta <- read.csv(file.path(PROC, "sample_metadata.csv"), check.names = FALSE)
meta <- meta[meta$qc_pass %in% c(TRUE, "True", "TRUE"), ]
stopifnot(!anyDuplicated(meta$title))

meta$case    <- factor(ifelse(meta$diagnosis == "Control", "control", "case"),
                       levels = c("control", "case"))
meta$sex     <- factor(meta$sex);   meta$tobacco <- factor(meta$tobacco)
meta$group   <- factor(meta$group); meta$plate   <- factor(meta$plate)
meta$lithium <- as.numeric(meta$lithium)

# The WGCNA covariate model: no diagnosis, no lithium -- those are the traits.
F_COV <- ~ age + sex + tobacco + group + rin + plate + seqpc1 + seqpc2 + seqpc3

# ------------------------------------------------- expression -> residuals
# Published WGCNA table defines the gene universe we must match.
ref_mod <- read.csv(file.path(REF, "File_S2_WGCNA_modules__Supplementary_File_2.csv"))
ref_key <- read.csv(file.path(REF, "File_S2_WGCNA_modules__module_name_color_key.csv"),
                    header = FALSE, col.names = c("module", "color"))
cat(sprintf("  published WGCNA table: %d genes, %d colours (%d modules + grey)\n",
            nrow(ref_mod), length(unique(ref_mod$moduleColor)),
            length(unique(ref_mod$moduleColor)) - 1L))

if (file.exists(CACHE_RESID)) {
  cat("Reusing cached residual matrix ...\n")
  cached  <- readRDS(CACHE_RESID)
  datExpr <- cached$datExpr; dropped9 <- cached$dropped9; gsg_flag <- cached$gsg_flag
  stopifnot(identical(rownames(datExpr), as.character(meta$title)))
} else {
  cat("Loading counts ...\n")
  counts <- as.matrix(read.delim(file.path(PROC, "counts_filtered.tsv"),
                                 row.names = 1, check.names = FALSE))
  # Join by key, never by position (titles all start with a digit).
  meta_ord <- meta[match(colnames(counts), meta$title), ]
  stopifnot(identical(as.character(meta_ord$title), colnames(counts)))
  cat(sprintf("  counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))

  design_cov <- model.matrix(F_COV, data = meta_ord)
  stopifnot(qr(design_cov)$rank == ncol(design_cov))

  cat("TMM + voom ...\n")
  y <- calcNormFactors(DGEList(counts = counts), method = "TMM")
  v <- voom(y, design_cov, plot = FALSE)
  rm(counts, y); invisible(gc())

  # Residualise on the covariates. One QR projection is the exact per-gene OLS
  # fit for all genes at once -- no loop, no approximation.
  cat("Residualising on covariates ...\n")
  Q     <- qr.Q(qr(design_cov))
  resid <- v$E - (v$E %*% Q) %*% t(Q)
  rm(v, Q); invisible(gc())

  # Which 9 genes did the authors drop between DE (12,353) and WGCNA (12,344)?
  dropped9 <- setdiff(rownames(resid), ref_mod$gene)
  stopifnot(length(setdiff(ref_mod$gene, rownames(resid))) == 0)

  # WGCNA's own gene screen is the obvious candidate explanation; test it.
  gsg <- goodSamplesGenes(t(resid), verbose = 0)
  gsg_flag <- rownames(resid)[!gsg$goodGenes]

  resid   <- resid[ref_mod$gene, , drop = FALSE]   # like-for-like 12,344
  datExpr <- t(resid)
  rm(resid); invisible(gc())

  # Reorder samples to the metadata order once, so every trait vector below can
  # be taken straight from `meta` without re-matching.
  datExpr <- datExpr[as.character(meta$title), , drop = FALSE]
  saveRDS(list(datExpr = datExpr, dropped9 = dropped9, gsg_flag = gsg_flag),
          CACHE_RESID)
}
cat(sprintf("  datExpr: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))

cat("\n=== The 9 genes in our DE set but not in their WGCNA set ===\n")
# Rank them by mean expression: if they were dropped for being lowly expressed,
# they should sit at the bottom of the 12,353.
de_ave <- read.csv(file.path(RES, "DE_lithium_all444.csv"))
de_ave$rank <- rank(de_ave$AveExpr)
d9 <- de_ave[match(dropped9, de_ave$gene), c("gene", "AveExpr", "rank")]
d9 <- d9[order(d9$AveExpr), ]
print(d9, row.names = FALSE, digits = 4)
cat(sprintf("  %d of 9 sit in the bottom 2%% of expression (rank <= %d of %d).\n",
            sum(d9$rank <= 0.02 * nrow(de_ave)), as.integer(0.02 * nrow(de_ave)),
            nrow(de_ave)))
cat(sprintf("  WGCNA goodSamplesGenes() flags %d of our 12,353 genes; %d of the 9.\n",
            length(gsg_flag), length(intersect(gsg_flag, dropped9))))
cat("  -> no zero-variance/missingness rule explains them; low expression is a\n")
cat("     partial pattern only (2 of the 9 are mid-range). The authors document\n")
cat("     no extra WGCNA filter, so the drop is unexplained by the methods.\n")

# ---------------------------------------------------------------- traits
lithium <- meta$lithium
bd      <- as.numeric(meta$case == "case")
nSamples <- nrow(datExpr)
cat(sprintf("  traits: lithium users %d / non-users %d ; cases %d / controls %d\n",
            sum(lithium == 1), sum(lithium == 0), sum(bd == 1), sum(bd == 0)))

GS_lith <- as.numeric(stats::cor(datExpr, lithium, use = "p"))
GS_bd   <- as.numeric(stats::cor(datExpr, bd,      use = "p"))
names(GS_lith) <- names(GS_bd) <- colnames(datExpr)

# Gene significance is computed straight off the residual matrix and does not depend
# on the clustering at all, so comparing it to their deposited GS.lithium column
# separates "is our expression matrix right" from "is our partition right". If this
# agrees and the partition does not, the disagreement is in the cut, not the data.
gs_pub <- setNames(ref_mod$GS.lithium, ref_mod$gene)[colnames(datExpr)]
cat(sprintf("\n=== Gene significance vs their deposited GS.lithium (%d genes) ===\n",
            length(gs_pub)))
cat(sprintf("  Pearson r = %.5f, Spearman = %.5f, sign concordance = %.4f\n",
            stats::cor(GS_lith, gs_pub), stats::cor(GS_lith, gs_pub, method = "spearman"),
            mean(sign(GS_lith) == sign(gs_pub))))

# ---------------------------------------------------------------- network
# The step-by-step WGCNA workflow, for the reasons given in the header. Every step
# here is also in blockwiseModules; what is absent is blockwiseModules' kME purge,
# its hard-coded detectCutHeight = 0.995, and its mergeCloseModules(relabel = TRUE).
if (file.exists(CACHE_NET)) {
  cat("\nReusing cached network ...\n")
  net <- readRDS(CACHE_NET)
  stopifnot(identical(net$genes, colnames(datExpr)))
} else {
  nthreads <- tryCatch({ enableWGCNAThreads(); WGCNAnThreads() },
                       error = function(e) { allowWGCNAThreads(); 1L })
  cat(sprintf("\nBuilding network (beta=%d, %s adjacency, %s TOM, %s correlation, %s threads)\n",
              BETA, NETWORK_TYPE, TOM_TYPE, COR_TYPE, nthreads))
  t0 <- Sys.time()
  # One dense 12,344 x 12,344 TOM; 1.2 GB, and the whole network in one block,
  # because a block-split network would not be comparable to theirs.
  dissTom <- TOMsimilarityFromExpr(datExpr, power = BETA, networkType = NETWORK_TYPE,
                                   TOMType = TOM_TYPE, corType = COR_TYPE,
                                   nThreads = nthreads, verbose = 0)
  dissTom <- 1 - dissTom
  invisible(gc())
  cat(sprintf("  TOM in %.1f min; clustering ...\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  geneTree <- fastcluster::hclust(as.dist(dissTom), method = "average")
  invisible(gc())

  dyn <- cutreeDynamic(dendro = geneTree, distM = dissTom, method = "hybrid",
                       deepSplit = DEEP_SPLIT,
                       cutHeight = if (is.na(DETECT_CUT_HEIGHT)) NULL else DETECT_CUT_HEIGHT,
                       minClusterSize = MIN_MODULE_SIZE,
                       pamStage = PAM_STAGE, pamRespectsDendro = PAM_RESPECTS_DENDRO,
                       verbose = 0)
  rm(dissTom); invisible(gc())
  cat(sprintf("  dynamic cut: %d modules, %d unassigned\n",
              length(setdiff(unique(dyn), 0)), sum(dyn == 0)))

  # labels2colors BEFORE merging, and relabel = FALSE, so that a merged branch keeps
  # the colour of one of its members instead of being renamed by rank. This is what
  # makes our colour names directly comparable to the published ones.
  merged <- mergeCloseModules(datExpr, labels2colors(dyn), cutHeight = MERGE_CUT_HEIGHT,
                              relabel = FALSE, useAbs = MERGE_USE_ABS, verbose = 0)
  net <- list(colors = as.character(merged$colors), dynamic = dyn,
              geneTree = geneTree, genes = colnames(datExpr))
  cat(sprintf("  done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(net, CACHE_NET)
}

ourColor <- net$colors
names(ourColor) <- colnames(datExpr)
modules  <- setdiff(sort(unique(ourColor)), "grey")
sizes    <- table(ourColor)[modules]

cat("\n=== Module structure ===\n")
cat(sprintf("  ours   : %d modules, sizes %d-%d, mean %.1f, unassigned(grey) %d\n",
            length(modules), min(sizes), max(sizes), mean(sizes),
            sum(ourColor == "grey")))
ref_sizes <- table(ref_mod$moduleColor)
ref_sizes_ng <- ref_sizes[names(ref_sizes) != "grey"]
cat(sprintf("  theirs : %d modules, sizes %d-%d, mean %.1f, unassigned(grey) %d\n",
            length(ref_sizes_ng), min(ref_sizes_ng), max(ref_sizes_ng),
            mean(ref_sizes_ng), ref_sizes[["grey"]]))
cat(sprintf("  paper  : 27 modules, sizes 48-2760, mean 441\n"))
cat(sprintf("           (441 is the mean over 28 groups incl. grey = %.1f;\n",
            mean(ref_sizes)))
cat(sprintf("            over the 27 real modules their own table gives %.1f)\n",
            mean(ref_sizes_ng)))
# mergeCloseModules never touches grey, so the unassigned count is decided entirely
# by the dynamic cut and is the single most diagnostic number in the whole run.
cat(sprintf("  dynamic cut before merging: %d modules, %d unassigned (their grey: %d)\n",
            length(setdiff(unique(net$dynamic), 0)), sum(net$dynamic == 0),
            ref_sizes[["grey"]]))

# -------------------------------------------------- eigengenes, GS, MM
MEs <- orderMEs(moduleEigengenes(datExpr, colors = ourColor)$eigengenes)
MEs <- MEs[, colnames(MEs) != "MEgrey", drop = FALSE]

MM   <- stats::cor(datExpr, MEs, use = "p")          # genes x modules

# Primary definition (reproduces their published r and p): eigengene vs trait.
r_lith <- as.numeric(stats::cor(MEs, lithium, use = "p"))
r_bd   <- as.numeric(stats::cor(MEs, bd,      use = "p"))
p_lith <- corPvalueStudent(r_lith, nSamples)
p_bd   <- corPvalueStudent(r_bd,   nSamples)
mod_of_ME <- sub("^ME", "", colnames(MEs))

# Methods' literal definition: correlation of module membership with gene
# significance, across the genes of the module.
mmgs <- function(color, gs) {
  idx <- which(ourColor == color)
  stats::cor(abs(MM[idx, paste0("ME", color)]), abs(gs[idx]), use = "p")
}
r_mmgs_lith <- vapply(mod_of_ME, mmgs, numeric(1), gs = GS_lith)

# Intramodular connectivity kIM, which the supplementary methods also report, is
# computed in src/06_wgcna_extended.R from this script's module assignment; it is not
# duplicated here. That script therefore has to be re-run after this one.

nMod  <- length(modules)
bonf  <- 0.05 / nMod
cat(sprintf("\n=== Module-trait association (Bonferroni 0.05/%d = %.3e) ===\n",
            nMod, bonf))

# ------------------------------------------------ DEG hypergeometric overlap
de <- read.csv(file.path(RES, "DE_lithium_all444.csv"))
de <- de[de$gene %in% colnames(datExpr), ]          # background = the 12,344
deg <- de$gene[de$adj.P.Val < 0.05]
N <- ncol(datExpr); K <- length(deg)
cat(sprintf("  lithium DEGs (FDR<0.05) in background: %d of %d genes\n", K, N))

ov  <- vapply(mod_of_ME, function(cl) sum(ourColor[deg] == cl), numeric(1))
msz <- vapply(mod_of_ME, function(cl) sum(ourColor == cl), numeric(1))
ov_p <- phyper(ov - 1, K, N - K, msz, lower.tail = FALSE)

# ------------------------------------------- correspondence with published
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sij <- sum(choose(tab, 2))
  sa  <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  expct <- sa * sb / choose(n, 2)
  (sij - expct) / ((sa + sb) / 2 - expct)
}
theirColor <- setNames(ref_mod$moduleColor, ref_mod$gene)[colnames(datExpr)]
ari_all <- ari(ourColor, theirColor)
keep_both <- ourColor != "grey" & theirColor != "grey"
ari_ng  <- ari(ourColor[keep_both], theirColor[keep_both])
cat(sprintf("\n=== Correspondence with their published moduleColor ===\n"))
cat(sprintf("  adjusted Rand index: %.4f (all genes) / %.4f (both non-grey, n=%d)\n",
            ari_all, ari_ng, sum(keep_both)))
# With relabel = FALSE the colour names are not arbitrary between the two runs, so a
# literal name match is meaningful rather than a coincidence of labelling.
cat(sprintf("  identical colour name: %d of %d genes (%.1f%%); shared colours %d of %d\n",
            sum(ourColor == theirColor), length(ourColor),
            100 * mean(ourColor == theirColor),
            length(intersect(modules, setdiff(unique(theirColor), "grey"))),
            length(ref_sizes_ng)))

xt <- table(ourColor, theirColor)
best <- apply(xt, 1, function(r) names(which.max(r)))
best_n <- apply(xt, 1, max)
their_n <- as.integer(ref_sizes[best])
jacc <- best_n / (as.integer(table(ourColor)[rownames(xt)]) + their_n - best_n)
names(best) <- names(best_n) <- rownames(xt); names(jacc) <- rownames(xt)
their_name <- setNames(ref_key$module, ref_key$color)

trait <- data.frame(
  module          = mod_of_ME,
  size            = as.integer(msz),
  r_lithium       = r_lith,
  p_lithium       = p_lith,
  r_bd            = r_bd,
  p_bd            = p_bd,
  deg_overlap     = as.integer(ov),
  deg_p           = ov_p,
  r_mm_gs_lithium = as.numeric(r_mmgs_lith),
  lithium_sig     = p_lith <= bonf,
  bd_sig          = p_bd   <= bonf,
  best_published  = ifelse(is.na(their_name[best[mod_of_ME]]), best[mod_of_ME],
                           paste0(their_name[best[mod_of_ME]], " (", best[mod_of_ME], ")")),
  match_overlap   = as.integer(best_n[mod_of_ME]),
  jaccard         = as.numeric(jacc[mod_of_ME]),
  row.names       = NULL)
trait <- trait[order(trait$p_lithium), ]

cat("\n  Modules passing Bonferroni for lithium use:\n")
sig <- trait[trait$lithium_sig, ]
if (nrow(sig) == 0) cat("    none\n") else
  print(sig[, c("module", "size", "r_lithium", "p_lithium", "deg_overlap",
                "deg_p", "best_published", "jaccard")],
        row.names = FALSE, digits = 3)
cat(sprintf("\n  Modules passing Bonferroni for BD: %d %s\n",
            sum(trait$bd_sig),
            if (sum(trait$bd_sig) == 0) "(paper: none -- agrees)" else
              paste0("(", paste(trait$module[trait$bd_sig], collapse = ", "), ")")))

# Their published modules are gene lists, so we can rebuild each one on OUR
# residual matrix. That separates "did the clustering reproduce" from "did the
# trait statistics reproduce" -- these numbers use their exact gene sets, so any
# disagreement here is in the statistic, not in the module detection.
cat("\n  Paper's five lithium modules, rebuilt on their exact gene sets:\n")
pub5 <- c(M1 = "brown", M7 = "lightsteelblue1", M9 = "plum2",
          M11 = "bisque4", M26 = "brown4")
pub_r   <- c(M1 = 0.156, M7 = -0.165, M9 = 0.153, M11 = 0.17, M26 = -0.175)
pub_p   <- c(M1 = 9.40e-4, M7 = 4.50e-4, M9 = 1.15e-3, M11 = 3.12e-4, M26 = 2.00e-4)
pub_ov  <- c(M1 = 431, M7 = 22, M9 = 17, M11 = 102, M26 = 17)
pub_ovp <- c(M1 = 2.03e-97, M7 = 1.00, M9 = 6.15e-7, M11 = 4.93e-13, M26 = 1.00)

cat(sprintf("    %-4s %-16s %5s | %-22s | %-22s | %s\n", "mod", "colour", "n",
            "r(eigengene,lithium)", "r(MM,GS) across genes", "DEG overlap"))
for (m in names(pub5)) {
  g  <- names(theirColor)[theirColor == pub5[m]]
  me <- moduleEigengenes(datExpr[, g, drop = FALSE],
                         colors = rep("x", length(g)))$eigengenes[[1]]
  rr <- stats::cor(me, lithium)
  # Sign of an eigengene is arbitrary; align it to the published direction.
  if (sign(rr) != sign(pub_r[m])) rr <- -rr
  pp <- corPvalueStudent(rr, nSamples)
  # The methods' literal definition, on the same gene set.
  mm2 <- stats::cor(datExpr[, g], me)
  r2  <- stats::cor(abs(as.numeric(mm2)), abs(GS_lith[g]))
  nov <- sum(g %in% deg)
  novp <- phyper(nov - 1, K, N - K, length(g), lower.tail = FALSE)
  cat(sprintf("    %-4s %-16s %5d | paper %+.3f %.1e        | %+.3f (n=%d genes)%s| paper %3d %.1e\n",
              m, pub5[m], length(g), pub_r[m], pub_p[m], r2, length(g),
              strrep(" ", max(1, 8 - nchar(sprintf("%d", length(g))))),
              pub_ov[m], pub_ovp[m]))
  cat(sprintf("    %-4s %-16s %5s | ours  %+.3f %.1e        |%s| ours  %3d %.1e\n",
              "", "", "", rr, pp, strrep(" ", 24), nov, novp))
}
cat("    -> the published r/p track the eigengene-vs-trait version (n=444),\n")
cat("       not the across-genes MM-vs-GS version their footnote describes.\n")

cat("\n  Best-matching module of ours for each published lithium module:\n")
for (m in names(pub5)) {
  ours_best <- rownames(xt)[which.max(xt[, pub5[m]])]
  cat(sprintf("    %-4s %-16s n=%4d -> our %-16s %d/%d genes shared\n",
              m, pub5[m], sum(theirColor == pub5[m]), ours_best,
              max(xt[, pub5[m]]), sum(theirColor == pub5[m])))
}
top_ov <- trait[order(-trait$deg_overlap), ][1, ]
cat(sprintf("    our largest DEG overlap: module %s, %d genes, p=%.2e\n",
            top_ov$module, top_ov$deg_overlap, top_ov$deg_p))

# ---------------------------------------------------------------- write
own_mm <- vapply(seq_len(ncol(datExpr)), function(i) {
  cl <- ourColor[i]
  if (cl == "grey") NA_real_ else MM[i, paste0("ME", cl)]
}, numeric(1))

write.csv(data.frame(gene = colnames(datExpr),
                     module = unname(ourColor),
                     membership = own_mm,
                     gs_lithium = unname(GS_lith),
                     gs_bd = unname(GS_bd),
                     published_module = unname(theirColor),
                     published_gs_lithium = unname(gs_pub),
                     row.names = NULL),
          file.path(RES, "wgcna_modules.csv"), row.names = FALSE)
write.csv(trait, file.path(RES, "wgcna_module_trait.csv"), row.names = FALSE)
cat("\nWrote results/wgcna_modules.csv and results/wgcna_module_trait.csv\n")
