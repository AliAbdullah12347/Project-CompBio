#!/usr/bin/env Rscript
# netmet-065  netmet-064 on covariate-residualised expression. A shim; the
# variant table lives in netmet-064.R so the two cannot drift apart.
commandArgs <- function(trailingOnly = FALSE) "netmet-065"
source("network/experiments/netmet-064.R")
