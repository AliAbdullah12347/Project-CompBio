# ==============================================================================
# boot_common.R -- the resampling engine shared by 03_bootstrap.R (subject
# resampling, for confidence intervals) and 04_permutation.R (label shuffling,
# for the null distribution).
#
# Both scripts need the SAME per-replicate computation, so it lives here once.
#
# WHAT IS REFITTED INSIDE EVERY REPLICATE, AND WHY
#
#   TMM normalisation  yes -- library-size factors are estimated from the data,
#                      so their uncertainty belongs inside the resampling.
#   voom weights       yes -- the mean-variance trend is estimated too. Holding
#                      the weights fixed at their full-sample values would treat
#                      an estimated quantity as known and understate the
#                      variance. It costs about 4.5 s of the 7 s per replicate;
#                      we have the time, so we pay it.
#   mediator models    yes -- refitted on the resampled subjects, unweighted.
#   centring           yes -- recomputed within the replicate, because the
#                      sample means are themselves estimates.
#
# The mediator models are UNWEIGHTED by design. voom weights describe count
# noise in gene expression; they say nothing about the precision of a subject's
# cell fractions, and the lithium-to-cell-mix effect is one gene-independent
# quantity. See METHODS.md for the consequence this has for the exact
# total = direct + indirect identity.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})

# One replicate. `idx` selects subjects (with replacement for the bootstrap, or
# 1:n for a permutation); `xx` is the exposure vector to use (resampled for the
# bootstrap, shuffled for the permutation).
#
# `with_p = TRUE` additionally returns delta-method p-values for NIE and TE, so
# the permutation null can count BH discoveries per replicate.
med_replicate <- function(idx, xx, CNT, M, C, lib = NULL, with_p = FALSE) {
  Mi <- M[idx, , drop = FALSE]
  Ci <- C[idx, , drop = FALSE]
  cM <- centre(Mi); cC <- centre(Ci)
  Xd <- med_design_full(xx, cM$X, cC$X)
  if (qr(Xd)$rank < ncol(Xd)) return(NULL)        # degenerate resample, skip

  dge <- normLibSizes(DGEList(CNT[, idx, drop = FALSE]), method = "TMM")
  v   <- voom(dge, Xd)

  mf <- mediator_fit(cM$X, xx, cC$X)
  of <- outcome_fit(v$E, v$weights, Xd, se = with_p, moderate = with_p)
  K  <- ncol(Mi); ix <- med_idx(K)
  th <- of$theta
  NDE <- as.vector(th[, ix$iX] + th[, ix$iI, drop = FALSE] %*% mf$b0)
  NIE <- as.vector((th[, ix$iM, drop = FALSE] + th[, ix$iI, drop = FALSE]) %*% mf$b1)
  out <- list(NIE = NIE, TE = NDE + NIE, b1 = mf$b1,
              b1_t = mf$b1 / mf$b1_se)
  if (with_p) {
    eff <- med_effects(of, mf, want_se = TRUE)
    out$NIE_p <- 2 * pnorm(-abs(eff$est$NIE / eff$se$NIE))
    out$TE_p  <- 2 * pnorm(-abs(eff$est$TE  / eff$se$TE))
  }
  out
}

# Build a PSOCK cluster (Windows has no fork) with everything the workers need.
# Four workers, not six: an earlier stage of this project thrashed a 7.7 GB
# machine at six and ran roughly five times slower than it should have.
make_med_cluster <- function(ncore = 4L, root, exports) {
  cl <- parallel::makeCluster(ncore)
  parallel::clusterExport(cl, "root", envir = list2env(list(root = root)))
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({library(edgeR); library(limma)})
    source(file.path(root, "mediation/scripts/med.R"))
    source(file.path(root, "mediation/scripts/boot_common.R"))
    NULL
  })
  parallel::clusterExport(cl, exports, envir = parent.frame())
  cl
}
