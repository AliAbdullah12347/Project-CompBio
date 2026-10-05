#!/usr/bin/env Rscript
# ==============================================================================
# 01_prep.R -- build every object the analysis starts from. Runs once.
#
# Reads only from ../data/cohort_474/ (read-only) and the GENCODE v19 GTF.
# Writes only inside de_analysis/.
#
# DECISIONS MADE HERE, AND WHY
#
# 1. Gene lengths come from the GENCODE v19 GTF by UNION OF EXONS. Exons of
#    different transcripts overlap, so their lengths cannot be summed; the
#    intervals are merged first. Union is also the correct denominator for
#    HTSeq union mode, which is how these counts were generated.
#
# 2. bMIND receives log2(TPM + 1). This follows Boltz et al., who log-transformed
#    because their largest value exceeded 50 TPM, and it keeps our cell-type
#    estimates on the same footing as the only published bMIND analysis of this
#    kind of data.
#
# 3. Cell-type proportions are LINEAGE-level (5 parts), per config.yaml. bMIND
#    takes proportions as input, so unreliable fine-grained proportions would
#    propagate into every expression estimate.
#
# 4. The gene filter runs ONCE, on the full 474, and the resulting gene set is
#    used for every contrast. This differs from the exploratory work, where the
#    filter was re-run per subset. Fixing it here is deliberate: the four
#    hypotheses compare DEG counts ACROSS contrasts and levels, and a Jaccard
#    or percentage computed over different gene universes is not comparable.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
COH  <- file.path(dirname(HERE), "data", "cohort_474")
GTF  <- file.path(dirname(HERE), "cibersortx", "data", "gencode.v19.annotation.gtf.gz")
OUT  <- file.path(HERE, "data"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("== 01_prep ==\n")

## ---- counts ---------------------------------------------------------------
counts <- as.matrix(read.delim(file.path(COH, "counts_474.tsv.gz"),
                               row.names = 1, check.names = FALSE))
counts <- counts[!grepl("^ENSGR", rownames(counts)), ]
rownames(counts) <- sub("\\.\\d+$", "", rownames(counts))
stopifnot(!any(duplicated(rownames(counts))))
cat(sprintf("counts %d x %d (ENSGR dropped)\n", nrow(counts), ncol(counts)))

## ---- metadata -------------------------------------------------------------
md <- read.csv(file.path(COH, "metadata_474_imputed.csv"), check.names = FALSE,
               colClasses = c(title = "character"))
stopifnot(identical(md$title, colnames(counts)))
m <- data.frame(
  title = md$title,
  dx      = factor(md[["bipolar disorder diagnosis"]], c("Control", "BP1", "BP2")),
  lithium = as.integer(md[["lithium use (non-user=0, user = 1)"]]),
  age = as.numeric(md$age), sex = factor(md$Sex),
  group = factor(md[["assessment group"]]), rin = as.numeric(md$rin),
  plate = factor(md[["sequencing plate"]]),
  seqpc1 = as.numeric(md[["sequencing metric pc1"]]),
  seqpc2 = as.numeric(md[["sequencing metric pc2"]]),
  seqpc3 = as.numeric(md[["sequencing metric pc3"]]),
  tobacco = as.integer(md$tobacco_imp_01),
  depth_ok = as.logical(md$depth_ok), rin_ok = as.logical(md$rin_ok),
  stringsAsFactors = FALSE)
TOB <- as.matrix(md[, grep("^tobacco_imp_\\d+$", names(md))]); rownames(TOB) <- md$title
cat(sprintf("metadata %d | %s\n", nrow(m),
            paste(sprintf("%s=%d", levels(m$dx), table(m$dx)), collapse = " ")))

## ---- gene filter, ONCE, on all 474 ----------------------------------------
keep <- rowSums(counts > 10) >= 0.90 * ncol(counts)
counts <- counts[keep, ]
cat(sprintf("gene filter (>10 in >=90%% of 474): %d genes -- fixed for all contrasts\n",
            nrow(counts)))

## ---- gene lengths, union of exons -----------------------------------------
lenfile <- file.path(OUT, "gene_lengths.rds")
if (file.exists(lenfile)) {
  len <- readRDS(lenfile); cat("gene lengths: cached\n")
} else {
  # Accumulate exon records into flat vectors, then split ONCE. Growing a
  # per-gene matrix inside the read loop is quadratic and does not finish.
  # Only the genes that survived the filter are kept, which drops ~78% of the
  # work immediately.
  cat("parsing GTF for union exon lengths ... ")
  want <- rownames(counts)
  con <- gzfile(GTF, "rt"); G <- S <- E <- vector("list", 64); k <- 0
  repeat {
    l <- readLines(con, n = 200000); if (!length(l)) break
    l <- l[!startsWith(l, "#")]; if (!length(l)) next
    p <- strsplit(l, "\t", fixed = TRUE)
    f3 <- vapply(p, `[`, character(1), 3)
    p <- p[f3 == "exon"]; if (!length(p)) next
    gid <- sub("\\..*$", "", sub('.*gene_id "([^"]+)".*', "\\1",
                                 vapply(p, `[`, character(1), 9)))
    ok <- gid %in% want; if (!any(ok)) next
    k <- k + 1
    G[[k]] <- gid[ok]
    S[[k]] <- as.integer(vapply(p, `[`, character(1), 4))[ok]
    E[[k]] <- as.integer(vapply(p, `[`, character(1), 5))[ok]
  }
  close(con)
  gid <- unlist(G, use.names = FALSE)
  st  <- unlist(S, use.names = FALSE); en <- unlist(E, use.names = FALSE)
  cat(sprintf("%d exon records ... ", length(gid)))
  # Merge overlapping intervals per gene, then sum the merged widths.
  union_len <- function(i) {
    o <- order(st[i]); a <- st[i][o]; b <- en[i][o]
    # A new block starts wherever this exon begins after the running maximum
    # end of everything before it.
    runmax <- cummax(c(-Inf, b[-length(b)]))
    grp <- cumsum(a > runmax + 1)
    sum(tapply(b, grp, max) - tapply(a, grp, min) + 1)
  }
  len <- vapply(split(seq_along(gid), gid), union_len, numeric(1))
  saveRDS(len, lenfile); cat(sprintf("%d genes\n", length(len)))
}
shared <- intersect(rownames(counts), names(len))
counts <- counts[shared, ]; len <- len[shared]
cat(sprintf("after length match: %d genes\n", nrow(counts)))

## ---- TPM and log-CPM ------------------------------------------------------
rate <- counts / (len / 1000)
tpm  <- t(t(rate) / colSums(rate)) * 1e6
logtpm <- log2(tpm + 1)
dge <- calcNormFactors(DGEList(counts), method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 3)
cat(sprintf("TPM column sums: %.0f to %.0f\n", min(colSums(tpm)), max(colSums(tpm))))

## ---- cell composition, lineage level --------------------------------------
fr <- read.csv(file.path(COH, "fractions_bmode_474_CIBERSORTx.csv"),
               check.names = FALSE, colClasses = c(Mixture = "character"))
rownames(fr) <- fr$Mixture
fr <- as.matrix(fr[m$title, setdiff(names(fr), c("Mixture", "P-value", "Correlation", "RMSE"))])
LIN <- list(
  gran = c("Neutrophils", "Eosinophils", "Mast cells resting", "Mast cells activated"),
  mono = c("Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
           "Dendritic cells resting", "Dendritic cells activated"),
  T    = c("T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
           "T cells CD4 memory activated", "T cells follicular helper",
           "T cells regulatory (Tregs)", "T cells gamma delta"),
  NK   = c("NK cells resting", "NK cells activated"),
  B    = c("B cells naive", "B cells memory", "Plasma cells"))
stopifnot(setequal(unlist(LIN), colnames(fr)))
L <- sapply(LIN, function(cs) rowSums(fr[, cs, drop = FALSE]))
L <- L / rowSums(L)
rownames(L) <- m$title

zrepl <- function(X, delta = 0.65) {
  lim <- apply(X, 2, function(v) min(v[v > 0]) * delta)
  z <- X <= 0; ins <- rowSums(sweep(z, 2, lim, "*"))
  X[z] <- matrix(rep(lim, each = nrow(X)), nrow(X))[z]
  X[!z] <- (X * (1 - ins))[!z]; X / rowSums(X)
}
Lz <- zrepl(L)
sbp <- list(b1 = list(n = c("gran", "mono"), d = c("T", "NK", "B")),
            b2 = list(n = "gran", d = "mono"),
            b3 = list(n = "T", d = c("NK", "B")),
            b4 = list(n = "NK", d = "B"))
psi <- t(sapply(sbp, function(b) {
  v <- setNames(numeric(ncol(Lz)), colnames(Lz))
  r <- length(b$n); s <- length(b$d)
  v[b$n] <- sqrt(s / (r * (r + s))); v[b$d] <- -sqrt(r / (s * (r + s))); v }))
stopifnot(max(abs(unname(psi %*% t(psi)) - diag(nrow(psi)))) < 1e-12)
ILR <- log(Lz) %*% t(psi); colnames(ILR) <- names(sbp); rownames(ILR) <- m$title
cat(sprintf("lineages (mean): %s\n",
            paste(sprintf("%s=%.3f", colnames(L), colMeans(L)), collapse = " ")))

## ---- contrast membership --------------------------------------------------
sel <- list(
  LI  = list(keep = m$dx == "BP1",
             grp  = factor(ifelse(m$lithium[m$dx == "BP1"] == 1, "case", "ctrl"),
                           c("ctrl", "case")),
             cov  = c("age", "sex", "tobacco", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")),
  BPD = list(keep = (m$dx == "BP1" & m$lithium == 0) | m$dx == "Control",
             grp  = NULL,
             cov  = c("age", "sex", "tobacco", "group", "rin", "plate",
                      "seqpc1", "seqpc2", "seqpc3")))
sel$BPD$grp <- factor(ifelse(m$dx[sel$BPD$keep] == "Control", "ctrl", "case"), c("ctrl", "case"))
for (n in names(sel))
  cat(sprintf("contrast %-4s n=%3d (%d ctrl vs %d case)\n", n, sum(sel[[n]]$keep),
              sum(sel[[n]]$grp == "ctrl"), sum(sel[[n]]$grp == "case")))

saveRDS(list(counts = counts, logtpm = logtpm, logcpm = logcpm, len = len,
             meta = m, TOB = TOB, frac = fr, lineage = L, lineage_z = Lz,
             ILR = ILR, sel = sel, LIN = LIN),
        file.path(OUT, "prep.rds"))
cat(sprintf("\nwrote data/prep.rds (%.1f MB)\n", file.size(file.path(OUT, "prep.rds")) / 1e6))
