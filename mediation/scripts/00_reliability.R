#!/usr/bin/env Rscript
# ==============================================================================
# 00_reliability.R -- a reliability coefficient for each ILR balance, which
# regression calibration needs and which we did not have.
#
# THE PROBLEM
#
# Correcting a mediator for measurement error needs to know how reliable the
# mediator is. On record we had only two numbers: median ICC 0.693 across the
# five lineages, and 0.907 for the myeloid/lymphoid balance. Nothing for b2, b3
# or b4 individually, and the multi-method fractions those came from are no
# longer on disk.
#
# THE SOLUTION, USING DATA WE ALREADY HAVE
#
# There are TWO independent deconvolutions of these same subjects:
#
#   1. Krebs et al.'s CIBERSORT fractions, deposited in the GEO series matrix
#      and carried in our metadata (22 LM22 types per sample).
#   2. Our own CIBERSORTx run, B-mode batch correction, quantile normalisation
#      off, 100 permutations.
#
# Same signature matrix, different software generation, different settings,
# different operators. Running both through the identical lineage aggregation
# and ILR construction gives two measurements of each balance per subject, and
# the intraclass correlation between them is a reliability estimate.
#
# WHICH ICC, AND WHY
#
# ICC(A,1): two-way random effects, ABSOLUTE agreement, single measurement.
# Absolute rather than consistency, because calibration cares whether the two
# give the same value, not merely whether they rank subjects the same way. A
# systematic offset between methods is real measurement error for our purpose.
#
# THE DIRECTION OF THE REMAINING BIAS, STATED PLAINLY
#
# Both estimates use the same LM22 signature and the same algorithm family, so
# they share whatever error the signature itself introduces. The ICC therefore
# OVERSTATES reliability. An overstated reliability produces an understated
# calibration correction, which keeps the corrected mediated proportion a lower
# bound -- the same direction every other caveat in this project points. That
# is why this is usable despite not being two genuinely independent methods.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation"))
dir.create("results", showWarnings = FALSE, recursive = TRUE)
options(width = 150)

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
md <- read.csv(file.path(ROOT, "data/cohort_474/metadata_474_imputed.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(md$title, p$meta$title))

## ---- the 22 LM22 types, as named in each source -----------------------------
# config.yaml's lineage map, with the metadata's lowercase/dot spellings beside
# the CIBERSORTx column names.
LIN <- list(
  gran = list(cx = c("Neutrophils","Eosinophils","Mast cells resting","Mast cells activated"),
              md = c("neutrophils","eosinophils","mast.cells.resting","mast.cells.activated")),
  mono = list(cx = c("Monocytes","Macrophages M0","Macrophages M1","Macrophages M2",
                     "Dendritic cells resting","Dendritic cells activated"),
              md = c("monocytes","macrophages.m0","macrophages.m1","macrophages.m2",
                     "dendritic.cells.resting","dendritic.cells.activated")),
  T    = list(cx = c("T cells CD8","T cells CD4 naive","T cells CD4 memory resting",
                     "T cells CD4 memory activated","T cells follicular helper",
                     "T cells regulatory (Tregs)","T cells gamma delta"),
              md = c("t.cells.cd8","t.cells.cd4.naive","t.cells.cd4.memory.resting",
                     "t.cells.cd4.memory.activated","t.cells.follicular.helper",
                     "t.cells.regulatory..tregs.","t.cells.gamma.delta")),
  NK   = list(cx = c("NK cells resting","NK cells activated"),
              md = c("nk.cells.resting","nk.cells.activated")),
  B    = list(cx = c("B cells naive","B cells memory","Plasma cells"),
              md = c("b.cells.naive","b.cells.memory","plasma.cells")))

## ---- our CIBERSORTx fractions ----------------------------------------------
fx <- read.csv(file.path(ROOT, "data/cohort_474/fractions_bmode_474_CIBERSORTx.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
rownames(fx) <- fx[[1]]
fx <- fx[p$meta$title, , drop = FALSE]

## ---- Krebs' deposited fractions, from the metadata --------------------------
mdnames <- names(md)
norm <- function(x) tolower(gsub("[^a-z0-9]+", ".", tolower(x)))
findcol <- function(want) {
  hit <- mdnames[norm(mdnames) == norm(want)]
  if (length(hit) == 1) return(hit)
  hit <- mdnames[startsWith(norm(mdnames), norm(want))]
  if (length(hit) >= 1) return(hit[1])
  NA_character_ }

cat("== 00_reliability ==\nmatching the 22 LM22 types across both sources\n\n")
agg <- function(src) {
  out <- sapply(names(LIN), function(L) {
    cols <- if (src == "cx") LIN[[L]]$cx else LIN[[L]]$md
    if (src == "cx") {
      got <- intersect(cols, names(fx))
      if (!length(got)) return(rep(NA_real_, nrow(fx)))
      rowSums(fx[, got, drop = FALSE])
    } else {
      got <- na.omit(vapply(cols, findcol, character(1)))
      if (!length(got)) return(rep(NA_real_, nrow(md)))
      rowSums(sapply(got, function(g) suppressWarnings(as.numeric(md[[g]]))), na.rm = TRUE)
    }
  })
  out <- out / rowSums(out)          # re-close to sum 1
  out }

CX <- agg("cx"); KR <- agg("md")
for (L in names(LIN)) {
  gotc <- length(intersect(LIN[[L]]$cx, names(fx)))
  gotm <- sum(!is.na(vapply(LIN[[L]]$md, findcol, character(1))))
  cat(sprintf("  %-5s CIBERSORTx %d/%d types | deposited %d/%d types\n",
              L, gotc, length(LIN[[L]]$cx), gotm, length(LIN[[L]]$md)))
}
cat(sprintf("\nlineage means  ours: %s\n", paste(sprintf("%s=%.3f", names(LIN), colMeans(CX)), collapse = " ")))
cat(sprintf("lineage means Krebs: %s\n", paste(sprintf("%s=%.3f", names(LIN), colMeans(KR)), collapse = " ")))

## ---- identical ILR construction applied to both -----------------------------
zrepl <- function(X, delta = 0.65) {
  lim <- apply(X, 2, function(v) { nz <- v[v > 0]; if (!length(nz)) 1e-6 else min(nz) }) * delta
  z <- X <= 0; ins <- rowSums(sweep(z, 2, lim, "*"))
  X[z] <- matrix(rep(lim, each = nrow(X)), nrow(X))[z]
  X[!z] <- (X * (1 - ins))[!z]; X / rowSums(X) }

sbp <- list(b1 = list(n = c("gran","mono"), d = c("T","NK","B")),
            b2 = list(n = "gran",           d = "mono"),
            b3 = list(n = "T",              d = c("NK","B")),
            b4 = list(n = "NK",             d = "B"))
psi <- t(sapply(sbp, function(b) {
  v <- setNames(numeric(length(LIN)), names(LIN))
  r <- length(b$n); s <- length(b$d)
  v[b$n] <- sqrt(s / (r * (r + s))); v[b$d] <- -sqrt(r / (s * (r + s))); v }))
stopifnot(max(abs(unname(psi %*% t(psi)) - diag(nrow(psi)))) < 1e-12)

ilr <- function(F) log(zrepl(F)) %*% t(psi)
I_cx <- ilr(CX); I_kr <- ilr(KR)
colnames(I_cx) <- colnames(I_kr) <- names(sbp)

## ---- ICC(A,1): two-way random effects, absolute agreement, single rating ----
icc_A1 <- function(x, y) {
  M <- cbind(x, y); n <- nrow(M); k <- ncol(M)
  gm <- mean(M); rm <- rowMeans(M); cm <- colMeans(M)
  MSR <- k * sum((rm - gm)^2) / (n - 1)
  MSC <- n * sum((cm - gm)^2) / (k - 1)
  SST <- sum((M - gm)^2)
  MSE <- (SST - k * sum((rm - gm)^2) - n * sum((cm - gm)^2)) / ((n - 1) * (k - 1))
  (MSR - MSE) / (MSR + (k - 1) * MSE + k * (MSC - MSE) / n) }

cat("\n=== reliability of each ILR balance ===\n")
R <- do.call(rbind, lapply(names(sbp), function(b) {
  x <- I_cx[, b]; y <- I_kr[, b]
  data.frame(balance = b,
             contrast = sprintf("%s vs %s", paste(sbp[[b]]$n, collapse = "+"),
                                paste(sbp[[b]]$d, collapse = "+")),
             icc_absolute = round(icc_A1(x, y), 4),
             pearson = round(cor(x, y), 4),
             spearman = round(cor(x, y, method = "spearman"), 4),
             sd_ours = round(sd(x), 4), sd_krebs = round(sd(y), 4),
             stringsAsFactors = FALSE) }))
print(R, row.names = FALSE)

cat("\n=== also, lineage proportions themselves (for comparison with the 0.693 on record) ===\n")
P <- do.call(rbind, lapply(names(LIN), function(L)
  data.frame(lineage = L, icc_absolute = round(icc_A1(CX[, L], KR[, L]), 4),
             mean_ours = round(mean(CX[, L]), 4), mean_krebs = round(mean(KR[, L]), 4),
             stringsAsFactors = FALSE)))
print(P, row.names = FALSE)
cat(sprintf("\n  median lineage ICC: %.4f   (0.693 was the figure on record)\n", median(P$icc_absolute)))

write.csv(R, "results/balance_reliability.csv", row.names = FALSE)
write.csv(P, "results/lineage_reliability.csv", row.names = FALSE)
saveRDS(list(ILR_cx = I_cx, ILR_kr = I_kr, frac_cx = CX, frac_kr = KR, reliability = R),
        "results/reliability.rds")
cat("\nwrote mediation/results/balance_reliability.csv\n")
