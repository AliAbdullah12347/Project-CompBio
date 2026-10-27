#!/usr/bin/env Rscript
# ==============================================================================
# 16_celltype_diversity.R -- alpha and beta diversity of the CELL-TYPE-SPECIFIC
# expression estimates (bMIND), for lithium (primary) and bipolar vs control
# (exploratory). Protocol: PRESPEC.md section 9 (addendum C). Handoff:
# CELLTYPE_DIVERSITY_HANDOFF.md.
#
# Runs from a fresh clone: every input is a TRACKED file (no .rds caches).
#
# Environment variables (all optional):
#   LINEAGES  comma list from gran,mono,T,NK,B          default: all five
#   NBOOT     alpha bootstrap replicates                default: 10000
#   NPERM     alpha permutation replicates              default: 10000
#   NPERM_V   vegan permutations (beta)                 default: 9999
#   NBOOT_V   bootstrap for the beta shared fraction    default: 1000
#   NCORES    cores for vegan permutation tests         default: 1
#   OUT       output folder                             default: results/celltype_diversity
#   COMBINE_ONLY=1  skip analysis; merge per-lineage outputs and apply Holm
#
# The script STOPS unless it first reproduces the published whole-blood result,
# so a differently assembled dataset can never produce plausible wrong numbers.
# ==============================================================================

t_start <- Sys.time()
here <- normalizePath(getwd(), winslash = "/")
ROOT <- if (basename(here) == "mediation") dirname(here) else here
stopifnot(file.exists(file.path(ROOT, "network/00_load.R")))
setwd(ROOT); set.seed(481)
suppressPackageStartupMessages({library(vegan); library(permute)})
source("network/00_load.R"); source("mediation/scripts/med.R")
options(width = 170)
say  <- function(...) cat(sprintf(...)); flush.console()
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("=", nchar(t))))
envi <- function(k, d) as.integer(Sys.getenv(k, d))
LIN_ALL <- c("gran", "mono", "T", "NK", "B"); INTERP <- c("gran", "mono", "T")
LIN   <- strsplit(Sys.getenv("LINEAGES", paste(LIN_ALL, collapse = ",")), ",")[[1]]
NBOOT <- envi("NBOOT", "10000"); NPERM <- envi("NPERM", "10000")
NPV   <- envi("NPERM_V", "9999"); NBV <- envi("NBOOT_V", "1000"); NCORES <- envi("NCORES", "1")
OUT   <- Sys.getenv("OUT", "mediation/results/celltype_diversity")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stopifnot(all(LIN %in% LIN_ALL))

## =============================================================================
## 0. assemble every input from tracked files, and check it
## =============================================================================
if (Sys.getenv("COMBINE_ONLY") != "1") {
rule("0. inputs")
d <- load_project()
titles <- d$meta$title
stopifnot(nrow(d$logtpm) == 12368, ncol(d$logtpm) == 474,
          identical(as.vector(table(d$groups)), c(234L, 74L, 152L)))
am <- function(f) read.csv(file.path("data/analysis_matrices", f), check.names = FALSE, stringsAsFactors = FALSE)
fz <- am("fractions_5lineages_zeroreplaced.csv"); fz <- fz[match(titles, fz[[1]]), ]
FZ <- as.matrix(fz[, LIN_ALL]); rownames(FZ) <- titles
stopifnot(!any(is.na(FZ)), all(FZ > 0))
ILR <- as.matrix(d$ilr[, c("b1", "b2", "b3", "b4")]); rownames(ILR) <- titles
gl <- am("gene_lengths_union_exons.csv"); len <- setNames(gl[[2]], gl[[1]])[rownames(d$counts)]
stopifnot(!any(is.na(len)))
md <- read.csv("data/cohort_474/metadata_474_imputed.csv", check.names = FALSE, stringsAsFactors = FALSE)
loglib <- log(setNames(as.numeric(md$library_size), md$title)[titles])
SE <- read.csv("data/analysis_matrices/bmind/bmind_mean_posterior_se.csv", check.names = FALSE)
SE <- as.matrix(SE[match(titles, SE$title), LIN_ALL]); rownames(SE) <- titles
bt <- read.delim("mediation/data/gencode_v19_gene_biotypes.tsv", stringsAsFactors = FALSE)
biotype <- setNames(bt$biotype, bt$gene)
stopifnot(sum(biotype == "protein_coding") > 20000)
g  <- rownames(d$counts)
is_pc <- biotype[g] %in% "protein_coding"; is_lm <- toupper(d$sym[g]) %in% toupper(d$lm22)
U  <- list(primary = g[is_pc & !is_lm], with_markers = g[is_pc])
stopifnot(length(U$primary) == 10288, length(U$with_markers) == 10626)
GS <- read.csv("mediation/data/genesets_msigdb_2026.1.Hs.csv", stringsAsFactors = FALSE)
HM <- GS[GS$collection == "hallmark", ]
say("  474 samples | primary universe %d genes | Hallmark %d sets | lineages this run: %s\n",
    length(U$primary), length(unique(HM$gs_name)), paste(LIN, collapse = ","))

read_lineage <- function(k) {
  f <- sprintf("data/analysis_matrices/bmind/bmind_logtpm_%s.csv.gz", k)
  m <- as.matrix(read.csv(gzfile(f), row.names = 1, check.names = FALSE))
  stopifnot(identical(rownames(m), g), setequal(colnames(m), titles))
  m[, titles]
}

## ---- design for each exposure (identical to the whole-blood scripts) --------
cfg_for <- function(which) {
  if (which == "lithium") {
    keep <- d$groups %in% c("bp_nolith", "bp_lith"); x <- as.integer(d$groups[keep] == "bp_lith")
  } else {
    keep <- d$groups %in% c("control", "bp_nolith"); x <- as.integer(d$groups[keep] == "bp_nolith")
  }
  mm <- d$meta[keep, ]
  mm$sex <- factor(mm$sex); mm$plate <- factor(mm$plate); mm$group <- factor(mm$group)
  f <- if (which == "lithium") ~ age + sex + tobacco + rin + plate + seqpc1 + seqpc2 + seqpc3
       else ~ age + sex + tobacco + rin + plate + seqpc1 + seqpc2 + seqpc3 + group
  C <- cbind(model.matrix(f, mm)[, -1, drop = FALSE], log_libsize = loglib[mm$title])
  C <- C[, apply(C, 2, function(v) length(unique(v)) > 1), drop = FALSE]
  colnames(C) <- make.names(colnames(C))
  list(name = which, titles = mm$title, x = x, C = C, M = ILR[mm$title, , drop = FALSE])
}
CFG <- list(lithium = cfg_for("lithium"), bipolar = cfg_for("bipolar"))
stopifnot(length(CFG$lithium$x) == 226, sum(CFG$lithium$x) == 152,
          length(CFG$bipolar$x) == 308, sum(CFG$bipolar$x) == 74)

## ---- estimators ----------------------------------------------------------------
fast_eff <- function(x, M, C, Y) {          # NIE / NDE / TE for many outcomes at once
  cM <- sweep(M, 2, colMeans(M)); cC <- sweep(C, 2, colMeans(C)); K <- ncol(M)
  Bm <- qr.coef(qr(cbind(1, x, cC)), cM)
  Th <- qr.coef(qr(cbind(1, x, cM, cM * x, cC)), Y)
  T2 <- Th[2 + seq_len(K), , drop = FALSE]; T3 <- Th[2 + K + seq_len(K), , drop = FALSE]
  NIE <- colSums((T2 + T3) * Bm[2, ]); NDE <- Th[2, ] + colSums(T3 * Bm[1, ])
  rbind(NIE = NIE, NDE = NDE, TE = NIE + NDE)
}
shannon <- function(P) -colSums(ifelse(P > 0, P * log(P), 0))
props   <- function(M) sweep(M, 2, colSums(M), "/")
pc_scores <- function(Y) { Yc <- scale(Y, scale = FALSE)
  s <- svd(Yc, nu = min(dim(Yc)) - 1, nv = 0); s$u %*% diag(s$d[seq_len(ncol(s$u))]) }
r2adj <- function(S, X) { n <- nrow(S); Sc <- scale(S, scale = FALSE); Xc <- scale(X, scale = FALSE)
  r2 <- sum(qr.fitted(qr(cbind(1, Xc)), Sc)^2) / sum(Sc^2); p <- qr(Xc)$rank
  1 - (1 - r2) * (n - 1) / (n - p - 1) }
shared <- function(S, x, C, M) { a <- r2adj(S, cbind(x, C)) - r2adj(S, C)
  b <- r2adj(S, cbind(x, C, M)) - r2adj(S, cbind(C, M)); c(R2_x = a, shared = 1 - b / a) }
ptest <- function(mod) { a <- anova(mod, permutations = how(nperm = NPV), parallel = NCORES)
  c(F = a$F[1], p = a$`Pr(>F)`[1]) }

alpha_tests <- function(cfg, Y, shrink) {
  x <- cfg$x; C <- cfg$C; M <- cfg$M; n <- length(x); Cs <- cbind(C, shrink)
  a1  <- t(apply(Y, 2, function(y) summary(lm(y ~ x + C))$coef[2, c(1, 4)]))
  a1s <- t(apply(Y, 2, function(y) summary(lm(y ~ x + Cs))$coef[2, c(1, 4)]))
  pint <- apply(Y, 2, function(y) { cM <- scale(M, scale = FALSE)
    anova(lm(y ~ x + cM + C), lm(y ~ x + cM + I(cM * x) + C))$`Pr(>F)`[2] })
  fe <- fast_eff(x, M, C, Y)
  bt <- replicate(NBOOT, { i <- sample.int(n, n, replace = TRUE)
    fast_eff(x[i], M[i, , drop = FALSE], C[i, , drop = FALSE], Y[i, , drop = FALSE])["NIE", ] })
  bt <- matrix(bt, ncol = ncol(Y), byrow = TRUE)
  pn <- matrix(replicate(NPERM, fast_eff(sample(x), M, C, Y)["NIE", ]), ncol = ncol(Y), byrow = TRUE)
  data.frame(outcome = colnames(Y), sd = apply(Y, 2, sd),
    TE = a1[, 1], TE_p = a1[, 2], TE_adj_shrink = a1s[, 1], TE_adj_shrink_p = a1s[, 2],
    NIE = fe["NIE", ], NIE_lo = apply(bt, 2, quantile, .025), NIE_hi = apply(bt, 2, quantile, .975),
    NIE_p_perm = sapply(seq_len(ncol(Y)), function(j) (1 + sum(abs(pn[, j]) >= abs(fe["NIE", j]))) / (1 + NPERM)),
    p_interaction = pint, stringsAsFactors = FALSE)
}

beta_tests <- function(cfg, S, shrink, label) {
  x <- cfg$x; C <- cfg$C; M <- cfg$M; n <- length(x)
  df <- data.frame(x = x, C, shrink, M, check.names = TRUE)
  fC <- paste(colnames(C), collapse = " + "); fS <- paste(colnames(shrink), collapse = " + ")
  fM <- paste(colnames(M), collapse = " + ")
  mk <- function(f) rda(as.formula(f), data = df)
  B1  <- ptest(mk(paste("S ~ x + Condition(", fC, ")")))
  B1s <- ptest(mk(paste("S ~ x + Condition(", fC, "+", fS, ")")))
  B3  <- ptest(mk(paste("S ~ x + Condition(", fC, "+", fM, ")")))
  Pb  <- ptest(mk(paste("S ~", fM, "+ Condition(x +", fC, ")")))
  B2  <- permutest(betadisper(dist(S), factor(x)), permutations = how(nperm = NPV))
  Sr  <- qr.resid(qr(cbind(1, C, shrink)), S)
  B2r <- permutest(betadisper(dist(Sr), factor(x)), permutations = how(nperm = NPV))
  sf  <- shared(S, x, C, M)
  bs  <- replicate(NBV, { i <- sample.int(n, n, replace = TRUE)
    shared(S[i, , drop = FALSE], x[i], C[i, , drop = FALSE], M[i, , drop = FALSE])["shared"] })
  bd <- betadisper(dist(S), factor(x))
  data.frame(layer = label, B1_F = B1["F"], B1_p = B1["p"], B1_adj_shrink_p = B1s["p"],
    R2adj_exposure = sf["R2_x"], B2_p = B2$tab$`Pr(>F)`[1], B2_resid_p = B2r$tab$`Pr(>F)`[1],
    disp_unexposed = mean(bd$distances[x == 0]), disp_exposed = mean(bd$distances[x == 1]),
    B3_direct_p = B3["p"], path_b_p = Pb["p"], shared = sf["shared"],
    shared_lo = quantile(bs, .025, na.rm = TRUE), shared_hi = quantile(bs, .975, na.rm = TRUE),
    stringsAsFactors = FALSE, row.names = NULL)
}

## =============================================================================
## 1. GATE: reproduce the published whole-blood result before anything else
## =============================================================================
rule("1. gate -- reproduce the whole-blood alpha result")
cl <- CFG$lithium
Pb <- props(d$counts[U$primary, cl$titles] / len[U$primary])
Hb <- shannon(Pb)
gate_te  <- summary(lm(Hb ~ cl$x + cl$C))$coef[2, ]
gate_nie <- fast_eff(cl$x, cl$M, cl$C, matrix(Hb, ncol = 1))["NIE", 1]
say("  whole-blood Shannon, lithium: TE %.4f (p = %.4f) | NIE %.4f\n", gate_te[1], gate_te[4], gate_nie)
say("  published (results/diversity/alpha_decomposition_DB.csv): TE -0.0405 (p = 0.0797) | NIE -0.0341\n")
if (abs(gate_te[1] + 0.0405) > 2e-4 || abs(gate_nie + 0.0341) > 2e-4)
  stop("GATE FAILED: the assembled data do not reproduce the whole-blood result. Do not continue.")
say("  GATE PASSED\n")

## =============================================================================
## 2-4. per lineage
## =============================================================================
for (k in LIN) {
  rule(sprintf("lineage %s%s", k, if (k %in% INTERP) "" else "  (FLAGGED: precision too low to interpret)"))
  t0 <- Sys.time(); A <- read_lineage(k)
  lin_props <- function(genes, ti) props(pmax(2^A[genes, ti, drop = FALSE] - 1, 0))
  for (cn in names(CFG)) {
    cfg <- CFG[[cn]]; ti <- cfg$titles
    shrink <- cbind(log_own_fraction = log(FZ[ti, k]), mean_posterior_se = SE[ti, k])
    P1 <- lin_props(U$primary, ti)
    Y <- cbind(H_primary = shannon(P1), logHill2 = log(1 / colSums(P1^2)),
               H_with_markers = shannon(lin_props(U$with_markers, ti)))
    # pathway Shannon, fractional counting over Hallmark within the marker-inclusive universe
    symg <- toupper(d$sym[U$with_markers]); ens <- U$with_markers[match(toupper(HM$gene_symbol), symg)]
    mem <- unique(data.frame(ens = ens, set = HM$gs_name)[!is.na(ens), ])
    kk <- table(mem$ens); W <- matrix(0, length(unique(mem$ens)), length(unique(mem$set)),
      dimnames = list(sort(unique(mem$ens)), sort(unique(mem$set))))
    W[cbind(mem$ens, mem$set)] <- 1 / as.numeric(kk[mem$ens])
    Pw <- t(W) %*% lin_props(U$with_markers, ti)[rownames(W), , drop = FALSE]
    Y <- cbind(Y, H_pw_hallmark = shannon(props(Pw)))
    ar <- alpha_tests(cfg, Y, shrink); ar$lineage <- k; ar$exposure <- cn
    write.csv(ar, file.path(OUT, sprintf("alpha_%s_%s.csv", k, cn)), row.names = FALSE)
    say("  %-8s alpha  H_primary: TE %+.4f p=%.4f | adj.shrink p=%.4f | NIE %+.4f [%+.4f, %+.4f]\n", cn,
        ar$TE[1], ar$TE_p[1], ar$TE_adj_shrink_p[1], ar$NIE[1], ar$NIE_lo[1], ar$NIE_hi[1])
    # beta: centred log2(TPM+1) profile (Aitchison-type) primary; Jensen-Shannon secondary
    Lg <- t(A[U$primary, ti, drop = FALSE]); Lg <- Lg - rowMeans(Lg)
    br <- beta_tests(cfg, pc_scores(Lg), shrink, "aitchison")
    Pt <- t(P1); Hs <- -rowSums(ifelse(Pt > 0, Pt * log(Pt), 0)); m <- nrow(Pt); J <- matrix(0, m, m)
    for (a in seq_len(m - 1)) { b <- (a + 1):m
      Mx <- (matrix(Pt[a, ], length(b), ncol(Pt), byrow = TRUE) + Pt[b, , drop = FALSE]) / 2
      J[b, a] <- J[a, b] <- sqrt(pmax(-rowSums(ifelse(Mx > 0, Mx * log(Mx), 0)) - (Hs[a] + Hs[b]) / 2, 0)) }
    Sj <- cmdscale(as.dist(J), k = m - 1, eig = FALSE)       # Euclidean embedding of the JS distance
    br <- rbind(br, beta_tests(cfg, Sj, shrink, "jensen_shannon"))
    br$lineage <- k; br$exposure <- cn
    write.csv(br, file.path(OUT, sprintf("beta_%s_%s.csv", k, cn)), row.names = FALSE)
    say("  %-8s beta   aitchison: B1 p=%.4f | adj.shrink p=%.4f | B2 p=%.4f (resid %.4f) | direct p=%.4f | shared %.2f\n",
        cn, br$B1_p[1], br$B1_adj_shrink_p[1], br$B2_p[1], br$B2_resid_p[1], br$B3_direct_p[1], br$shared[1])
  }
  say("  lineage %s done in %.1f min\n", k, as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
}

## =============================================================================
## 5. combine and apply the pre-declared Holm correction
## =============================================================================
rule("5. combine")
af <- list.files(OUT, "^alpha_.*\\.csv$", full.names = TRUE); bf <- list.files(OUT, "^beta_.*\\.csv$", full.names = TRUE)
if (length(af) < 10 || length(bf) < 10) {
  say("  %d alpha and %d beta files present (10 each needed): run the remaining lineages, then COMBINE_ONLY=1\n",
      length(af), length(bf))
} else {
  # colClasses matters: read.csv would otherwise turn the lineage name "T" into TRUE
  rd <- function(f) read.csv(f, stringsAsFactors = FALSE, colClasses = c(lineage = "character", exposure = "character"))
  AR <- do.call(rbind, lapply(af, rd)); BR <- do.call(rbind, lapply(bf, rd))
  stopifnot(setequal(unique(AR$lineage), LIN_ALL), setequal(unique(BR$lineage), LIN_ALL))
  AR$interpretable <- AR$lineage %in% INTERP; BR$interpretable <- BR$lineage %in% INTERP
  AR$TE_p_holm <- NA; BR$B1_p_holm <- NA
  for (cn in unique(AR$exposure)) {
    i <- AR$exposure == cn & AR$outcome == "H_primary" & AR$interpretable
    AR$TE_p_holm[i] <- p.adjust(AR$TE_p[i], "holm")
    j <- BR$exposure == cn & BR$layer == "aitchison" & BR$interpretable
    BR$B1_p_holm[j] <- p.adjust(BR$B1_p[j], "holm")
  }
  write.csv(AR, file.path(OUT, "ALPHA_ALL.csv"), row.names = FALSE)
  write.csv(BR, file.path(OUT, "BETA_ALL.csv"), row.names = FALSE)
  pr <- AR[AR$outcome == "H_primary", c("exposure", "lineage", "TE", "TE_p", "TE_p_holm", "TE_adj_shrink_p", "NIE", "NIE_lo", "NIE_hi", "NIE_p_perm")]
  cat("\nPRIMARY alpha (Shannon, protein-coding non-LM22):\n"); print(pr[order(pr$exposure, pr$lineage), ], row.names = FALSE, digits = 3)
  pb <- BR[BR$layer == "aitchison", c("exposure", "lineage", "B1_p", "B1_p_holm", "B1_adj_shrink_p", "B2_p", "B2_resid_p", "B3_direct_p", "shared")]
  cat("\nPRIMARY beta (Aitchison):\n"); print(pb[order(pb$exposure, pb$lineage), ], row.names = FALSE, digits = 3)
  say("\n  wrote %s/ALPHA_ALL.csv and BETA_ALL.csv\n", OUT)
}
say("\ntotal %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
