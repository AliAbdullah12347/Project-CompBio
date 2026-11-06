#!/usr/bin/env Rscript
# ===========================================================================
# eos-052 — is orangered3's surviving lithium association a real
#           non-compositional effect, or composition the 5-lineage ILR
#           cannot see?
#
# grp-051b found that conditioning on the four ILR balances takes the 29
# lithium-associated modules down to 1: orangered3 (37 genes, padj 0.0059),
# which is also the most lithium-DEG-enriched module in the network
# (OR 11.4). Its control hub genes are IL5RA, IL1RL1, PTGDR2, ALOX15, CEBPE,
# LGALS12 — eosinophil and basophil markers.
#
# The ILR coordinates resolve five lineages, with eosinophils and basophils
# folded into "granulocytes" alongside neutrophils. So a shift *within* the
# granulocyte compartment is invisible to that adjustment. If orangered3
# tracks such a shift, the module is composition after all and the
# adjustment was simply too coarse.
#
# CIRCULARITY. orangered3's genes are themselves eosinophil markers, so an
# eosinophil score built from expression would be tautological. The tests
# below use the deposited CIBERSORTx fractions, which are estimated from the
# LM22 signature, and the overlap between orangered3 and LM22 is measured and
# reported so the remaining circularity is visible rather than assumed away.
# ===========================================================================
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
source("network/00_load.R"); d <- load_project(cache = TRUE, verbose = FALSE)
OUT <- "network/runs/eos-052"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, f) write.csv(x, file.path(OUT, f), row.names = FALSE)
say <- function(...) cat(paste0(..., collapse = ""), "\n")

g   <- read.csv("network/runs/grp-051/genes_control_A.csv")
stopifnot(identical(g$gene, rownames(d$logtpm)))
f22 <- read.csv("data/analysis_matrices/fractions_LM22_22types.csv", check.names = FALSE)
f22 <- f22[match(colnames(d$logtpm), f22$sample), ]
stopifnot(identical(f22$sample, colnames(d$logtpm)))

E  <- expr_for(d, group = NULL)
md <- d$meta
md$li <- as.integer(d$groups == "bp_lith")
md$dx <- as.integer(d$groups %in% c("bp_nolith", "bp_lith"))
keep <- !is.na(d$groups)
ilr <- d$ilr[, c("b1", "b2", "b3", "b4")]

# ---------------------------------------------------------------------------
# 0. what the fine fractions can and cannot support
# ---------------------------------------------------------------------------
fine <- c("Eosinophils", "Mast cells resting", "Neutrophils", "Monocytes")
qual <- do.call(rbind, lapply(fine, function(cn) data.frame(
  cell_type = cn, mean = mean(f22[[cn]]), median = median(f22[[cn]]),
  pct_nonzero = round(100 * mean(f22[[cn]] > 0), 1))))
wr(qual, "fraction_quality.csv"); say("== estimability of the fine fractions =="); print(qual, row.names = FALSE)

# ---------------------------------------------------------------------------
# 1. does lithium shift the granulocyte compartment's composition?
#    Eosinophils are zero in most samples, so they are tested as presence /
#    absence; a linear model on a mostly-zero variable is not meaningful.
# ---------------------------------------------------------------------------
cv <- data.frame(dx = md$dx, age = md$age, sex = md$sex, rin = md$rin,
                 plate = factor(md$plate))[keep, ]
li <- md$li[keep]
tests <- list()
for (cn in c("Neutrophils", "Mast cells resting", "Monocytes")) {
  y <- f22[[cn]][keep]
  s <- summary(lm(y ~ ., data = cbind(li = li, cv)))$coefficients
  tests[[cn]] <- data.frame(measure = cn, test = "lm, adj dx+age+sex+rin+plate",
                            beta = s["li", 1], p = s["li", 4])
}
eo <- as.integer(f22$Eosinophils[keep] > 0)
s  <- summary(glm(eo ~ ., data = cbind(li = li, cv), family = binomial))$coefficients
tests[["Eosinophils"]] <- data.frame(measure = "Eosinophils detected (yes/no)",
  test = "logistic, adj dx+age+sex+rin+plate", beta = s["li", 1], p = s["li", 4])
# neutrophil share OF the granulocyte compartment: the quantity the 5-lineage
# ILR is blind to, since it treats granulocytes as one part.
gran <- f22$Neutrophils + f22$Eosinophils + f22$`Mast cells resting`
nshare <- f22$Neutrophils / gran
s <- summary(lm(nshare[keep] ~ ., data = cbind(li = li, cv)))$coefficients
tests[["nshare"]] <- data.frame(measure = "neutrophil share of granulocytes",
  test = "lm, adj dx+age+sex+rin+plate", beta = s["li", 1], p = s["li", 4])
T1 <- do.call(rbind, tests); T1$padj <- p.adjust(T1$p, "BH")
wr(T1, "lithium_vs_fine_composition.csv")
say("\n== does lithium move the fine composition? =="); print(T1, row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. orangered3: is it a coherent module at all, and does it track eosinophils?
# ---------------------------------------------------------------------------
or3 <- g$gene[g$module == "orangered3"]
say("\n== orangered3: ", length(or3), " genes ==")
say("   symbols: ", paste(head(g$symbol[g$module == "orangered3"], 20), collapse = ", "))

# circularity: how much of the module is LM22 signature?
lm22 <- d$lm22
in22 <- sum(g$symbol[g$module == "orangered3"] %in% lm22)
say(sprintf("   LM22 signature genes inside the module: %d of %d (%.0f%%)",
            in22, length(or3), 100 * in22 / length(or3)))

# coherence per group: variance explained by the module's own first PC. A real
# module holds together in every group; an artefact of one clustering does not.
coh <- do.call(rbind, lapply(c("control", "bp_nolith", "bp_lith"), function(gr) {
  idx <- which(d$groups == gr); X <- E[idx, or3, drop = FALSE]
  pv <- prcomp(scale(X))$sdev^2
  data.frame(group = gr, n = length(idx), pc1_var_pct = round(100 * pv[1] / sum(pv), 1))
}))
wr(coh, "orangered3_coherence.csv")
say("\n== orangered3 coherence by group =="); print(coh, row.names = FALSE)

ME <- moduleEigengenes(E, colors = g$module, excludeGrey = TRUE)$eigengenes
me <- ME$MEorangered3
cr <- do.call(rbind, lapply(c("Eosinophils", "Mast cells resting", "Neutrophils", "Monocytes"), function(cn)
  data.frame(cell_type = cn,
             spearman_r = round(cor(me[keep], f22[[cn]][keep], method = "spearman"), 3),
             p = cor.test(me[keep], f22[[cn]][keep], method = "spearman",
                          exact = FALSE)$p.value)))
cr <- rbind(cr, data.frame(cell_type = "neutrophil share of granulocytes",
  spearman_r = round(cor(me[keep], nshare[keep], method = "spearman"), 3),
  p = cor.test(me[keep], nshare[keep], method = "spearman", exact = FALSE)$p.value))
wr(cr, "orangered3_vs_fractions.csv")
say("\n== orangered3 eigengene vs the fine fractions =="); print(cr, row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. the decisive test: does the lithium effect on orangered3 survive once the
#    granulocyte compartment's internal composition is in the model?
# ---------------------------------------------------------------------------
base <- cbind(li = li, cv, ilr[keep, ])
mods <- list(
  "lithium only"                      = data.frame(li = li, cv[, c("dx","age","sex","rin","plate")]),
  "+ 5-lineage ILR (grp-051b)"        = base,
  "+ mast/basophil fraction"          = cbind(base, mast = f22$`Mast cells resting`[keep]),
  "+ eosinophils detected"            = cbind(base, eo = eo),
  "+ neutrophil share of granulocytes"= cbind(base, nshare = nshare[keep]),
  "+ all three"                       = cbind(base, mast = f22$`Mast cells resting`[keep],
                                              eo = eo, nshare = nshare[keep]))
T3 <- do.call(rbind, lapply(names(mods), function(nm) {
  s <- summary(lm(me[keep] ~ ., data = mods[[nm]]))$coefficients
  data.frame(model = nm, li_beta = round(s["li", 1], 5), li_p = s["li", 4])
}))
T3$pct_attenuation <- round(100 * (1 - abs(T3$li_beta) / abs(T3$li_beta[1])), 1)
wr(T3, "orangered3_nested_models.csv")
say("\n== lithium effect on orangered3, as composition is resolved more finely ==")
print(T3, row.names = FALSE)

# same nested ladder across every lithium-significant module, so this is not a
# story about one module
mt <- read.csv("network/runs/grp-051/module_trait_joint.csv")
sig <- mt$module[!is.na(mt$li_padj) & mt$li_padj < 0.05]
T4 <- do.call(rbind, lapply(sig, function(m) {
  y <- ME[[paste0("ME", m)]][keep]
  p <- function(X) summary(lm(y ~ ., data = X))$coefficients["li", 4]
  data.frame(module = m, p_ilr = p(base),
             p_ilr_fine = p(cbind(base, mast = f22$`Mast cells resting`[keep],
                                  eo = eo, nshare = nshare[keep])))
}))
T4$padj_ilr <- p.adjust(T4$p_ilr, "BH"); T4$padj_fine <- p.adjust(T4$p_ilr_fine, "BH")
wr(T4[order(T4$p_ilr), ], "all_lithium_modules_fine_adjustment.csv")
say(sprintf("\n== of %d lithium modules: %d survive 5-lineage ILR, %d survive ILR + fine granulocyte terms",
            nrow(T4), sum(T4$padj_ilr < 0.05), sum(T4$padj_fine < 0.05)))
say("done")
