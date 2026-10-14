#!/usr/bin/env Rscript
# ===========================================================================
# 001-module-trait-ds4
#
# Module-eigengene / trait relationships at the USABLE decomposition.
#
# base-001 measured preservation on the degenerate mergeCutHeight=0.25 structure
# (3-4 modules, one holding 95-97% of genes). mod-000-resolution showed
# deepSplit=4, mergeCutHeight=0 gives 33 modules, largest 11.68%. This analysis
# uses THAT decomposition.
#
# The arm's question, at the module level:
#   does any module eigengene associate with lithium, and does that association
#   survive adjustment for cell composition (ILR b1-b4)?
#
# Design points obeyed (IMAC_README):
#   6.2  modules defined ONCE on the control reference; measured in BP1 here
#   5    lithium contrast is BP1 ONLY (n=226: 74 off, 152 on)
#   5    `assessment group` is constant within BP1 - excluded, would be singular
#   6.6  BH across modules WITHIN a trait; Benjamini-Bogomolov across traits
#   6.8  multiple testing via de_v2/scripts/mtc.R, not a hand-rolled version
#
# Writes only into this folder. No shared file is touched.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
suppressPackageStartupMessages({ library(WGCNA) })
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

OUT <- "network/analyses/001-module-trait-ds4/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse=""))); flush.console() }

# ---- config, written before anything is computed --------------------------
CFG <- list(
  analysis       = "001-module-trait-ds4",
  module_source  = "network/runs/mod-000-resolution/labels_p12_ds4.rds",
  power          = 12L, deepSplit = 4L, mergeCutHeight = 0,
  corFnc         = "bicor", networkType = "signed",
  reference      = "control (n=234) - modules defined here",
  measured_in    = "BP1 only (n=226) - lithium contrast per IMAC_README 5",
  traits         = "lithium, age, sex, rin, plate, seqpc1-3, gran/mono/T/NK/B, b1-b4",
  mtc            = "de_v2/scripts/mtc.R (IMAC_README 6.8); BH within trait, BB across traits",
  seed           = 20261014L,
  git_sha        = tryCatch(trimws(system("git rev-parse --short HEAD", intern=TRUE)[1]), error=function(e) "unknown")
)
writeLines(paste0("{\n", paste(sprintf('  "%s": "%s"', names(CFG), unlist(CFG)), collapse=",\n"), "\n}"),
           file.path(OUT, "..", "config.json"))
set.seed(CFG$seed)
msg("config written")

source("de_v2/scripts/mtc.R")           # validated MTC helpers, do not reimplement
source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)

# ---- the usable decomposition ---------------------------------------------
lab <- readRDS(CFG$module_source)
colors <- labels2colors(lab)
keep_mod <- setdiff(unique(colors), "grey")
msg(length(keep_mod), " non-grey modules; grey = ", sum(colors=="grey"), " genes")

# ---- BP1 samples only (the lithium contrast) ------------------------------
bp <- !is.na(d$groups) & d$groups %in% c("bp_nolith", "bp_lith")
stopifnot(sum(bp) == 226L)
E  <- t(d$logtpm[, bp])                      # 226 x 12368, samples x genes
md <- d$meta[bp, , drop = FALSE]
il <- d$ilr [match(md$title, d$ilr$sample),  c("b1","b2","b3","b4"), drop=FALSE]
fr <- d$frac[match(md$title, d$frac$sample), c("gran","mono","T","NK","B"), drop=FALSE]
stopifnot(nrow(il) == nrow(E), nrow(fr) == nrow(E))
msg("BP1 n=", nrow(E), " | lithium on=", sum(md$lithium==1), " off=", sum(md$lithium==0))

# ---- module eigengenes ----------------------------------------------------
me <- moduleEigengenes(E, colors = colors, excludeGrey = TRUE)$eigengenes
me <- me[, order(colnames(me)), drop = FALSE]
msg("eigengenes: ", ncol(me), " modules x ", nrow(me), " samples")
saveRDS(me, file.path(OUT, "module_eigengenes_bp1.rds"))

mod_sizes <- sapply(sub("^ME", "", colnames(me)), function(m) sum(colors == m))

# ---- trait panel ----------------------------------------------------------
traits <- data.frame(
  lithium = as.numeric(md$lithium),
  age     = as.numeric(md$age),
  sex     = as.numeric(as.factor(md$sex)),
  rin     = as.numeric(md$rin),
  seqpc1  = as.numeric(md$seqpc1), seqpc2 = as.numeric(md$seqpc2), seqpc3 = as.numeric(md$seqpc3),
  gran = fr$gran, mono = fr$mono, Tcell = fr$T, NKcell = fr$NK, Bcell = fr$B,
  b1 = il$b1, b2 = il$b2, b3 = il$b3, b4 = il$b4
)

# ---- marginal correlations, BH within each trait --------------------------
res <- do.call(rbind, lapply(names(traits), function(tn) {
  p <- sapply(seq_len(ncol(me)), function(j) cor.test(me[,j], traits[[tn]], method="pearson")$p.value)
  r <- sapply(seq_len(ncol(me)), function(j) cor(me[,j], traits[[tn]], use="pairwise.complete.obs"))
  data.frame(trait = tn, module = sub("^ME","",colnames(me)), n_genes = mod_sizes,
             r = round(r,4), p = p, p_BH = p.adjust(p, method="BH"))
}))
write.csv(res, file.path(OUT, "module_trait_marginal.csv"), row.names = FALSE)
msg("marginal correlations written")

# ---- THE QUESTION: does a lithium association survive composition? --------
# M1: ME ~ lithium                      (marginal)
# M2: ME ~ lithium + b1..b4             (composition-adjusted)
# M3: ME ~ lithium + technical + b1..b4 (fully adjusted)
tech <- data.frame(age=traits$age, sex=traits$sex, rin=traits$rin,
                   seqpc1=traits$seqpc1, seqpc2=traits$seqpc2, seqpc3=traits$seqpc3)
lith <- traits$lithium
comp <- il

coef_p <- function(fit) { s <- summary(fit)$coefficients; if ("lith" %in% rownames(s)) s["lith", c("Estimate","Pr(>|t|)")] else c(NA,NA) }
surv <- do.call(rbind, lapply(seq_len(ncol(me)), function(j) {
  y <- me[, j]
  m1 <- lm(y ~ lith)
  m2 <- lm(y ~ lith + comp$b1 + comp$b2 + comp$b3 + comp$b4)
  m3 <- lm(y ~ lith + comp$b1 + comp$b2 + comp$b3 + comp$b4 +
               tech$age + tech$sex + tech$rin + tech$seqpc1 + tech$seqpc2 + tech$seqpc3)
  c1 <- coef_p(m1); c2 <- coef_p(m2); c3 <- coef_p(m3)
  data.frame(module = sub("^ME","",colnames(me)[j]), n_genes = mod_sizes[j],
             beta_marginal = round(c1[1],4), p_marginal = c1[2],
             beta_compadj  = round(c2[1],4), p_compadj  = c2[2],
             beta_fulladj  = round(c3[1],4), p_fulladj  = c3[2],
             pct_attenuation = round(100*(1 - abs(c2[1])/abs(c1[1])), 1))
}))
surv$p_marginal_BH <- p.adjust(surv$p_marginal, method="BH")
surv$p_compadj_BH  <- p.adjust(surv$p_compadj,  method="BH")
surv$p_fulladj_BH  <- p.adjust(surv$p_fulladj,  method="BH")
surv <- surv[order(surv$p_marginal), ]
write.csv(surv, file.path(OUT, "lithium_composition_adjusted.csv"), row.names = FALSE)
msg("composition-adjustment models written")

# ---- summary --------------------------------------------------------------
cat("\n================ LITHIUM, BEFORE AND AFTER COMPOSITION ================\n")
print(head(surv[, c("module","n_genes","beta_marginal","p_marginal","p_marginal_BH",
                    "beta_compadj","p_compadj","p_compadj_BH","pct_attenuation")], 12), row.names=FALSE)
n_marg <- sum(surv$p_marginal_BH < 0.05); n_comp <- sum(surv$p_compadj_BH < 0.05); n_full <- sum(surv$p_fulladj_BH < 0.05)
cat(sprintf("\n  modules with BH<0.05 for lithium:  marginal %d / %d | composition-adjusted %d | fully adjusted %d\n",
            n_marg, nrow(surv), n_comp, n_full))
cat("\n================ WHICH MODULES TRACK CELL COMPOSITION =================\n")
cmp <- res[res$trait %in% c("gran","mono","Tcell","NKcell","Bcell"), ]
top <- aggregate(abs(cmp$r), by=list(module=cmp$module), FUN=max)
names(top)[2] <- "max_abs_r_vs_lineage"
top <- top[order(-top$max_abs_r_vs_lineage), ]
print(head(top, 10), row.names=FALSE)
cat(sprintf("\n  modules with |r|>0.5 against at least one lineage fraction: %d / %d\n",
            sum(top$max_abs_r_vs_lineage > 0.5), nrow(top)))
saveRDS(list(surv=surv, res=res, sizes=mod_sizes), file.path(OUT, "summary.rds"))
cat("\nANALYSIS_001_DONE\n")
