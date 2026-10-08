#!/usr/bin/env Rscript
# ==============================================================================
# 09_ztest_pm.R -- the simplest possible version of the question "is our
# proportion mediated bigger than the permutation null?"
#
# Standardises the observed median proportion mediated against the distribution
# of the SAME statistic computed on the SAME genes with the lithium labels
# shuffled (04b_perm_pm.R). z = (observed - mean(null)) / sd(null).
#
# Reported alongside the exact empirical permutation p, because with 300 draws
# the empirical p is the more trustworthy of the two and the z-test is only a
# summary of it.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); options(width = 150)
say <- function(...) cat(sprintf(...))

P <- readRDS("results/perm/perm_pm_fixedset.rds")
B <- list(lineage = readRDS("results/boot/lineage_boot.rds"),
          proposal = readRDS("results/boot/proposal_boot.rds"))
L <- readRDS("results/sens/calibration.rds")
pc <- readRDS("results/sens/calibration_proposal.rds")
corrected <- c(lineage = L$cal$corrected[2], proposal = pc$summary$PM_corrected[2])

out <- NULL
for (vn in names(P)) {
  o <- P[[vn]]; nul <- o$null
  mu <- mean(nul); sdv <- sd(nul)
  z <- (o$obs - mu) / sdv
  p_z <- 2 * pnorm(-abs(z))
  # normality of the null is not guaranteed, so report the exact empirical p too
  out <- rbind(out, data.frame(
    tree = vn, n_genes = o$n_sig, n_perm = length(nul),
    observed = o$obs, null_mean = mu, null_sd = sdv,
    z = z, p_z = p_z, p_empirical = o$p,
    boot_se = B[[vn]]$pm_se, stringsAsFactors = FALSE))
}
say("== z-test: observed median proportion mediated vs the permutation null ==\n\n")
print(transform(out, observed = round(observed, 4), null_mean = round(null_mean, 4),
                null_sd = round(null_sd, 4), z = round(z, 3),
                p_z = round(p_z, 3), p_empirical = round(p_empirical, 3),
                boot_se = round(boot_se, 4)), row.names = FALSE)

say("\n  Shapiro-Wilk test of null normality (does the z-test's assumption hold?):\n")
for (vn in names(P))
  say("    %-9s W = %.4f, p = %.4f  -> %s\n", vn,
      shapiro.test(P[[vn]]$null)$statistic, shapiro.test(P[[vn]]$null)$p.value,
      if (shapiro.test(P[[vn]]$null)$p.value > 0.05) "normal enough" else
        "NOT normal; trust the empirical p, not the z")

say("\n  Does the measurement-error correction change the verdict?\n")
for (vn in names(P)) {
  f <- corrected[[vn]] / P[[vn]]$obs
  say("    %-9s correction multiplies PM by about %.3f. Applied to BOTH the observed\n", vn, f)
  say("              value and the null, z is unchanged at %.3f: a common rescaling\n",
      (P[[vn]]$obs - mean(P[[vn]]$null)) / sd(P[[vn]]$null))
  say("              cancels in the numerator and denominator.\n")
}
write.csv(out, "results/perm/pm_ztest.csv", row.names = FALSE)
say("\n  wrote results/perm/pm_ztest.csv\n")
