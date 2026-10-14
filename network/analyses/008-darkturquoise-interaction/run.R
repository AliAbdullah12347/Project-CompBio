#!/usr/bin/env Rscript
# ===========================================================================
# 008-darkturquoise-interaction
#
# Analysis 007 found that the composition-darkturquoise correlation differs
# across groups (control r_b2=-0.13 vs BP r_b2=-0.44). If the slopes differ,
# a simple additive composition adjustment is invalid — the control-lith
# difference surviving adjustment (analysis 007) could be an artefact of
# using the wrong model.
#
# Pre-specified:
#   (1) LM: dt_ME ~ group * b2 + b1 + b3 + b4 (group×b2 interaction)
#   (2) LRT: compare to additive model
#   (3) If interaction survives: repeat pairwise within interaction model
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(WGCNA) })
enableWGCNAThreads()
OUT <- "network/analyses/008-darkturquoise-interaction/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)

msg("Loading data from analysis 007 sample_data.csv ...")
df <- read.csv("network/analyses/007-darkturquoise-across-groups/results/sample_data.csv",
               stringsAsFactors = FALSE)
df$group <- factor(df$group, levels = c("control","bp_nolith","bp_lith"))
msg(nrow(df), " samples, groups: ", paste(table(df$group), collapse="/"))

# ------- Group slopes of b2 on dt_ME ----------------------------------------
msg("Per-group slopes of dt_ME ~ b2 ...")
slopes <- do.call(rbind, lapply(levels(df$group), function(g) {
  s <- df[df$group == g, ]
  m <- lm(dt_ME ~ b2, data = s)
  coefs <- coef(m)
  ci    <- confint(m, "b2")
  data.frame(group = g, n = nrow(s),
             intercept  = round(coefs[1], 4),
             slope_b2   = round(coefs[2], 4),
             ci_lo      = round(ci[1], 4),
             ci_hi      = round(ci[2], 4),
             r_b2       = round(cor(s$dt_ME, s$b2), 3))
}))
cat("\n======= Per-group b2 slopes =======\n")
print(slopes, row.names = FALSE)
write.csv(slopes, file.path(OUT, "per_group_slopes.csv"), row.names = FALSE)

# ------- Interaction model --------------------------------------------------
msg("Fitting interaction model: dt_ME ~ group * b2 + b1 + b3 + b4 ...")
m_add  <- lm(dt_ME ~ group + b1 + b2 + b3 + b4, data = df)
m_int  <- lm(dt_ME ~ group * b2 + b1 + b3 + b4, data = df)
lrt <- anova(m_add, m_int)
cat("\n======= LRT: additive vs group × b2 interaction =======\n")
print(lrt)
write.csv(as.data.frame(lrt), file.path(OUT, "lrt_interaction.csv"), row.names = TRUE)

cat(sprintf("\nInteraction F(%d,%d) = %.3f, p = %.5f\n",
    lrt$Df[2], lrt$Res.Df[2], lrt$F[2], lrt$"Pr(>F)"[2]))

# Interaction coefficients
cat("\n======= Interaction model coefficients =======\n")
sm <- summary(m_int)
coef_df <- as.data.frame(coef(sm))
coef_df$term <- rownames(coef_df)
print(coef_df[, c("term","Estimate","Std. Error","t value","Pr(>|t|)")],
      row.names = FALSE)
write.csv(coef_df, file.path(OUT, "interaction_coefs.csv"), row.names = FALSE)

# ------- Pairwise within interaction model ----------------------------------
msg("Pairwise contrasts within interaction model ...")
# Using the interaction model, test group contrasts holding b2 at its mean
b2_mean <- mean(df$b2)
pairs <- list(
  c("control","bp_nolith"),
  c("control","bp_lith"),
  c("bp_nolith","bp_lith")
)

# Predict at b2=mean, b1=mean, b3=mean, b4=mean for each group
b1m <- mean(df$b1); b3m <- mean(df$b3); b4m <- mean(df$b4)
newdf <- do.call(rbind, lapply(levels(df$group), function(g) {
  data.frame(group=g, b1=b1m, b2=b2_mean, b3=b3m, b4=b4m, stringsAsFactors=FALSE)
}))
newdf$group <- factor(newdf$group, levels = levels(df$group))
preds <- predict(m_int, newdata=newdf, se.fit=TRUE)
pred_df <- data.frame(group=newdf$group,
                      fit=round(preds$fit,5), se=round(preds$se.fit,5))

cat("\n======= Predicted dt_ME at mean composition, by group =======\n")
print(pred_df, row.names=FALSE)
write.csv(pred_df, file.path(OUT, "marginal_means.csv"), row.names = FALSE)

# Proper contrast tests using vcov() — contrast at b2=mean, b1=mean, b3=mean, b4=mean
# Model: dt_ME ~ group*b2 + b1 + b3 + b4
# Coef order: Intercept, groupbp_nolith, groupbp_lith, b2, b1, b3, b4, gnolith:b2, glith:b2
b2m <- mean(df$b2)
vc  <- vcov(m_int)
cf  <- coef(m_int)

# Predict marginal mean for each group at mean composition
# control: Intercept + b2*b2m + b1*b1m + b3*b3m + b4*b4m
# nolith:  Intercept + groupbp_nolith + (b2 + gnolith:b2)*b2m + ...
# lith:    Intercept + groupbp_lith   + (b2 + glith:b2)*b2m   + ...
contrast_ctrl  <- c(1, 0, 0, b2m, b1m, b3m, b4m, 0,   0  )  # control
contrast_nolith <- c(1, 1, 0, b2m, b1m, b3m, b4m, b2m, 0  )  # nolith
contrast_lith  <- c(1, 0, 1, b2m, b1m, b3m, b4m, 0,   b2m)  # lith

pred_ctrl   <- sum(cf * contrast_ctrl)
pred_nolith <- sum(cf * contrast_nolith)
pred_lith   <- sum(cf * contrast_lith)
se_ctrl     <- sqrt(as.numeric(t(contrast_ctrl)   %*% vc %*% contrast_ctrl))
se_nolith   <- sqrt(as.numeric(t(contrast_nolith) %*% vc %*% contrast_nolith))
se_lith     <- sqrt(as.numeric(t(contrast_lith)   %*% vc %*% contrast_lith))

cat(sprintf("\nMarginal means at mean composition:\n  control = %.5f (SE %.5f)\n  bp_nolith = %.5f (SE %.5f)\n  bp_lith = %.5f (SE %.5f)\n",
            pred_ctrl, se_ctrl, pred_nolith, se_nolith, pred_lith, se_lith))

# Contrast: lith vs control
contrast_diff_lith_ctrl <- contrast_lith - contrast_ctrl
diff_lc   <- as.numeric(cf %*% contrast_diff_lith_ctrl)
se_lc     <- sqrt(as.numeric(t(contrast_diff_lith_ctrl) %*% vc %*% contrast_diff_lith_ctrl))
t_lc      <- diff_lc / se_lc
# df from interaction model residuals
df_resid  <- m_int$df.residual
p_lc      <- 2 * pt(-abs(t_lc), df = df_resid)

# Contrast: nolith vs control
contrast_diff_nc  <- contrast_nolith - contrast_ctrl
diff_nc   <- as.numeric(cf %*% contrast_diff_nc)
se_nc     <- sqrt(as.numeric(t(contrast_diff_nc) %*% vc %*% contrast_diff_nc))
t_nc      <- diff_nc / se_nc
p_nc      <- 2 * pt(-abs(t_nc), df = df_resid)

# Contrast: lith vs nolith
contrast_diff_ln <- contrast_lith - contrast_nolith
diff_ln   <- as.numeric(cf %*% contrast_diff_ln)
se_ln     <- sqrt(as.numeric(t(contrast_diff_ln) %*% vc %*% contrast_diff_ln))
t_ln      <- diff_ln / se_ln
p_ln      <- 2 * pt(-abs(t_ln), df = df_resid)

pwise_int <- data.frame(
  contrast = c("control vs bp_nolith","control vs bp_lith","bp_nolith vs bp_lith"),
  diff     = round(c(diff_nc, diff_lc, diff_ln), 5),
  se       = round(c(se_nc, se_lc, se_ln), 5),
  t        = round(c(t_nc, t_lc, t_ln), 3),
  p_raw    = round(c(p_nc, p_lc, p_ln), 5),
  p_bonf   = round(pmin(1, c(p_nc, p_lc, p_ln) * 3), 5)
)
cat("\n======= Pairwise contrasts in interaction model (Bonferroni n=3) =======\n")
print(pwise_int, row.names = FALSE)
write.csv(pwise_int, file.path(OUT, "pairwise_interaction_model.csv"), row.names = FALSE)

cat("\nANALYSIS_008_DONE\n")
