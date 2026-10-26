#!/usr/bin/env Rscript
# ==============================================================================
# 11_diversity_prep.R -- build every diversity quantity once, for all 474
# subjects, so the test scripts only fit models. (PRESPEC.md section 8.)
#
# WHAT IS BUILT
#
#   gene universes   primary = filtered AND protein-coding AND not LM22
#                    sensitivity: + LM22 markers; all filtered genes (any
#                    biotype); primary minus mitochondrial genes
#   alpha            Shannon H (nats), Hill q=1 = exp(H), Hill q=2 = 1/sum p^2,
#                    on TPM proportions; Shannon on read-count proportions
#                    (sensitivity); Hill q=0 = richness on ALL protein-coding
#                    genes, because on the filtered set every gene is detected
#                    in >=90% of samples by construction, which would make
#                    richness constant
#   beta inputs      centred log-ratio (CLR) matrices for Aitchison distance;
#                    TPM proportions for Jensen-Shannon and Bray-Curtis
#   pathways         MSigDB Hallmark (primary) and KEGG legacy (sensitivity),
#                    frozen into data/ so a future msigdbr update cannot change
#                    them; abundances by fractional counting
#   cell-mix         Shannon diversity of the 5 lineage fractions (and of the
#                    proposal's 6 parts)
#   decomposition    H(bulk) = within-lineage + between-lineage (JSD) +
#                    residual, from the bMIND per-sample lineage profiles
#
# Every check that can be asserted, is.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(msigdbr)})
options(width = 160)
say  <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("=", nchar(t))))
dir.create("results/diversity", showWarnings = FALSE, recursive = TRUE)

prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
titles <- prep$meta$title
CNT <- prep$counts[, titles]                        # 12,368 filtered x 474
len <- readRDS(file.path(ROOT, "de_analysis/data/gene_lengths.rds"))
stopifnot(identical(names(len), rownames(CNT)))
e2s <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"), header = FALSE,
                  col.names = c("ens", "sym"), stringsAsFactors = FALSE)
sym <- setNames(toupper(e2s$sym), e2s$ens)
lm22 <- unique(toupper(trimws(read.delim(file.path(ROOT, "data/LM22.txt"),
                                         check.names = FALSE)[[1]])))

## ---- biotypes from the GENCODE v19 annotation already on disk ---------------
rule("1. gene universes")
# Streamed in chunks and filtered to "gene" lines as it goes: the annotation is
# ~1.2 GB uncompressed, and loading it whole would starve the resampling workers.
# Cached to a small TSV on first run.
bt_file <- "data/gencode_v19_gene_biotypes.tsv"
if (!file.exists(bt_file)) {
  con <- gzfile(file.path(ROOT, "cibersortx/data/gencode.v19.annotation.gtf.gz"), "r")
  keep <- character(0)
  repeat {
    ch <- readLines(con, n = 200000L)
    if (!length(ch)) break
    keep <- c(keep, ch[grepl("\tgene\t", ch, fixed = TRUE)])
  }
  close(con)
  att <- vapply(strsplit(keep, "\t", fixed = TRUE), `[`, character(1), 9)
  gid <- sub("[.].*$", "", sub('.*gene_id "([^"]+)".*', "\\1", att))
  gty <- sub('.*gene_type "([^"]+)".*', "\\1", att)
  write.table(data.frame(gene = gid, biotype = gty), bt_file, sep = "\t",
              quote = FALSE, row.names = FALSE)
  rm(keep, att)
}
bt <- read.delim(bt_file, stringsAsFactors = FALSE)
biotype <- setNames(bt$biotype, bt$gene)
# guard against exactly the failure that happened once: a mis-escaped
# backreference wrote a control character instead of every biotype
stopifnot(sum(biotype == "protein_coding") > 20000)
say("  GENCODE v19 genes with a biotype: %d | protein-coding: %d\n",
    length(biotype), sum(biotype == "protein_coding"))

g     <- rownames(CNT)
is_pc <- biotype[g] %in% "protein_coding"
is_lm <- sym[g] %in% lm22
is_mt <- grepl("^MT-", sym[g])
U <- list(
  primary      = g[is_pc & !is_lm],
  with_markers = g[is_pc],
  all_filtered = g,
  no_mito      = g[is_pc & !is_lm & !is_mt])
for (u in names(U)) say("  %-13s %6d genes\n", u, length(U[[u]]))
stopifnot(!any(is.na(is_pc)), length(U$primary) > 10000)

## ---- library size = all assigned reads, from the deposited metadata ---------
md <- read.csv(file.path(ROOT, "data/cohort_474/metadata_474_imputed.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
libsize <- setNames(as.numeric(md$library_size), md$title)[titles]
stopifnot(!any(is.na(libsize)))

## ---- alpha diversity ---------------------------------------------------------
rule("2. alpha diversity")
props <- function(M) sweep(M, 2, colSums(M), "/")
tpm_p <- function(genes) props(CNT[genes, , drop = FALSE] / len[genes])
shannon <- function(P) -colSums(ifelse(P > 0, P * log(P), 0))
hill2   <- function(P) 1 / colSums(P^2)

alpha <- data.frame(title = titles, libsize = libsize, log_libsize = log(libsize),
                    stringsAsFactors = FALSE)
P_primary <- tpm_p(U$primary)
alpha$H_primary      <- shannon(P_primary)
alpha$Hill1_primary  <- exp(alpha$H_primary)
alpha$Hill2_primary  <- hill2(P_primary)
alpha$H_with_markers <- shannon(tpm_p(U$with_markers))
alpha$H_all_filtered <- shannon(tpm_p(U$all_filtered))
alpha$H_no_mito      <- shannon(tpm_p(U$no_mito))
alpha$H_counts       <- shannon(props(CNT[U$primary, , drop = FALSE]))

# Hill q = 0 on ALL protein-coding non-LM22 genes in the raw matrix
raw <- read.delim(gzfile(file.path(ROOT, "data/raw_as_deposited/GSE124326_count_matrix.txt.gz")),
                  check.names = FALSE, row.names = 1)
raw <- raw[!grepl("^ENSGR", rownames(raw)), ]
rownames(raw) <- sub("[.].*$", "", rownames(raw))
colnames(raw) <- sub("[.]counts$", "", colnames(raw))
raw <- raw[, titles]
rg <- rownames(raw)
u0 <- rg[biotype[rg] %in% "protein_coding" & !(sym[rg] %in% lm22)]
alpha$Hill0_richness <- colSums(raw[u0, ] > 0)
say("  richness universe: %d protein-coding non-LM22 genes in the raw matrix\n", length(u0))
# library size check: the deposited value should be the column total
say("  deposited library size vs raw column totals: max relative diff %.2g\n",
    max(abs(libsize - colSums(raw)) / libsize))

# checks
stopifnot(max(abs(colSums(P_primary) - 1)) < 1e-12)
stopifnot(all(alpha$Hill1_primary >= alpha$Hill2_primary - 1e-9))   # Hill numbers fall with q
say("  Hill numbers decrease with q in every sample: PASS\n")
say("  Shannon H (primary): mean %.4f nats, range %.4f-%.4f | exp(H) = %.0f effective genes (range %.0f-%.0f)\n",
    mean(alpha$H_primary), min(alpha$H_primary), max(alpha$H_primary),
    mean(alpha$Hill1_primary), min(alpha$Hill1_primary), max(alpha$Hill1_primary))
cc <- cor(alpha[, c("H_primary", "H_with_markers", "H_all_filtered", "H_no_mito", "H_counts")])
say("  correlation of primary H with: markers-in %.3f | all-filtered %.3f | no-mito %.4f | counts %.3f\n",
    cc[1, 2], cc[1, 3], cc[1, 4], cc[1, 5])
# how much do the structural ncRNAs drive the all-genes version?
nc <- g[sym[g] %in% c("RN7SL1", "RN7SL2", "RN7SK", "RPPH1", "RMRP")]
alpha$share_struct_ncRNA <- colSums(CNT[nc, ]) / colSums(CNT)
say("  structural ncRNA share vs H(all filtered): r = %.3f | vs H(primary): r = %.3f\n",
    cor(alpha$share_struct_ncRNA, alpha$H_all_filtered), cor(alpha$share_struct_ncRNA, alpha$H_primary))

## ---- beta inputs -------------------------------------------------------------
rule("3. beta-diversity inputs")
clr <- function(genes) {          # samples x genes; gene lengths cancel in CLR distances
  L <- log(t(CNT[genes, , drop = FALSE]) + 0.5)
  L - rowMeans(L)
}
CLR <- list(primary = clr(U$primary), with_markers = clr(U$with_markers))
for (u in names(CLR)) say("  CLR %-13s %d samples x %d genes\n", u, nrow(CLR[[u]]), ncol(CLR[[u]]))

## ---- pathways ---------------------------------------------------------------
rule("4. pathway layer")
gs_file <- "data/genesets_msigdb_2026.1.Hs.rds"
if (!file.exists(gs_file)) {
  H <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "H"))
  K <- as.data.frame(msigdbr(species = "Homo sapiens", collection = "C2",
                             subcollection = "CP:KEGG_LEGACY"))
  saveRDS(list(hallmark = H[, c("gs_name", "gene_symbol", "db_version")],
               kegg_legacy = K[, c("gs_name", "gene_symbol", "db_version")],
               frozen = "2026-10-18", msigdbr_version = as.character(packageVersion("msigdbr"))),
          gs_file)
}
GS <- readRDS(gs_file)
say("  frozen gene sets: %s (MSigDB %s, msigdbr %s)\n", gs_file,
    unique(GS$hallmark$db_version), GS$msigdbr_version)

pathway_abund <- function(sets, genes) {
  # map symbols ONLY within this universe, so a symbol shared by two Ensembl IDs
  # cannot resolve to the one that is not in the universe
  symg <- sym[genes]
  ens <- genes[match(toupper(sets$gene_symbol), symg)]
  s2 <- sets[!is.na(ens), ]; ens <- ens[!is.na(ens)]
  memb <- unique(data.frame(ens = ens, set = s2$gs_name, stringsAsFactors = FALSE))
  k <- table(memb$ens)                                   # how many sets each gene is in
  W <- matrix(0, length(unique(memb$ens)), length(unique(memb$set)),
              dimnames = list(sort(unique(memb$ens)), sort(unique(memb$set))))
  W[cbind(memb$ens, memb$set)] <- 1 / as.numeric(k[memb$ens])   # fractional counting
  P <- tpm_p(genes)                                      # TPM proportions over the universe
  A <- t(W) %*% P[rownames(W), , drop = FALSE]           # sets x samples
  list(A = A, n_genes = nrow(W), n_sets = ncol(W),
       mass_covered = colSums(P[rownames(W), , drop = FALSE]))
}
PW <- list(
  hallmark         = pathway_abund(GS$hallmark, U$with_markers),
  hallmark_noLM22  = pathway_abund(GS$hallmark, U$primary),
  kegg             = pathway_abund(GS$kegg_legacy, U$with_markers))
for (nm in names(PW)) {
  pw <- PW[[nm]]
  say("  %-16s %3d sets | %4d genes mapped | median %.1f%% of TPM mass covered\n",
      nm, pw$n_sets, pw$n_genes, 100 * median(pw$mass_covered))
  Pp <- props(pw$A)
  alpha[[paste0("H_pw_", nm)]] <- shannon(Pp)
}
PCLR <- lapply(PW, function(pw) { L <- log(t(pw$A)); L - rowMeans(L) })
stopifnot(all(sapply(PW, function(pw) all(pw$A > 0))))
say("  every pathway abundance is strictly positive, so CLR needs no pseudocount: PASS\n")

## ---- cell-mix diversity (interpretation C) -----------------------------------
rule("5. diversity of the cell mix itself")
lz <- prep$lineage_z[titles, , drop = FALSE]
stopifnot(max(abs(rowSums(lz) - 1)) < 1e-9, all(lz > 0))
alpha$H_cellmix5 <- -rowSums(lz * log(lz))
dp <- readRDS("data/med_input_proposal.rds")
pz <- dp$parts_z[titles, , drop = FALSE]
alpha$H_cellmix6 <- -rowSums(pz * log(pz))
say("  Shannon of the 5 lineage fractions: mean %.3f (max possible log 5 = %.3f)\n",
    mean(alpha$H_cellmix5), log(5))
say("  Shannon of the proposal's 6 parts : mean %.3f (max possible log 6 = %.3f)\n",
    mean(alpha$H_cellmix6), log(6))

## ---- within / between-lineage decomposition (interpretation B) ---------------
rule("6. H(bulk) = within-lineage + between-lineage + residual")
bm <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
stopifnot(identical(dimnames(bm$A)[[2]], colnames(lz)))
say("  bMIND lineages: %s | scale log2(TPM+1), un-logged here\n", paste(colnames(lz), collapse = ", "))
gi <- match(U$primary, dimnames(bm$A)[[1]]); si <- match(titles, dimnames(bm$A)[[3]])
stopifnot(!any(is.na(gi)), !any(is.na(si)))
Hk <- matrix(NA_real_, length(titles), ncol(lz), dimnames = list(titles, colnames(lz)))
Hmix <- numeric(length(titles))
for (j in seq_along(titles)) {
  Pk <- pmax(2^bm$A[gi, , si[j]] - 1, 0)                 # genes x lineages, linear
  Pk <- sweep(Pk, 2, colSums(Pk), "/")
  Hk[j, ] <- -colSums(ifelse(Pk > 0, Pk * log(Pk), 0))
  Pm <- Pk %*% lz[j, ]                                   # the cell-mix-weighted mixture
  Hmix[j] <- -sum(ifelse(Pm > 0, Pm * log(Pm), 0))
}
alpha$H_within  <- rowSums(Hk * lz)
alpha$JSD_between <- Hmix - alpha$H_within
alpha$H_mixture <- Hmix
alpha$resid_decomp <- alpha$H_primary - Hmix
for (k in colnames(lz)) alpha[[paste0("H_lineage_", k)]] <- Hk[, k]
stopifnot(all(alpha$JSD_between >= -1e-10))
stopifnot(all(alpha$JSD_between <= alpha$H_cellmix5 + 1e-10))
say("  JSD >= 0 in every sample: PASS | JSD <= H(cell mix) in every sample: PASS\n")
say("  mean H(bulk) %.4f = within %.4f + between %.4f + residual %+.4f\n",
    mean(alpha$H_primary), mean(alpha$H_within), mean(alpha$JSD_between), mean(alpha$resid_decomp))
say("  residual: SD %.4f vs SD of H(bulk) %.4f | corr(H bulk, H mixture) = %.3f\n",
    sd(alpha$resid_decomp), sd(alpha$H_primary), cor(alpha$H_primary, Hmix))

## ---- save --------------------------------------------------------------------
saveRDS(list(alpha = alpha, CLR = CLR, PCLR = PCLR, P_primary = P_primary,
             universes = U, u0_richness = u0, pathways = lapply(PW, function(p) p[c("n_sets", "n_genes")]),
             titles = titles),
        "data/diversity_input.rds")
write.csv(alpha, "results/diversity/per_sample_diversity.csv", row.names = FALSE)
write.csv(data.frame(universe = names(U), n_genes = lengths(U)),
          "results/diversity/gene_universes.csv", row.names = FALSE)
say("\n  wrote data/diversity_input.rds and results/diversity/per_sample_diversity.csv (%d x %d)\n",
    nrow(alpha), ncol(alpha))
