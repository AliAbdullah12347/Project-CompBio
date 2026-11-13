#!/usr/bin/env Rscript
# netmet-062  netmet-060 with the expression residualised for age, sex,
# tobacco, RIN, plate and the three sequencing PCs before any network is
# built. A shim, not a copy: the variant table lives in netmet-060.R so the
# two can never drift apart. See that file for why this exists.
commandArgs <- function(trailingOnly = FALSE) "netmet-062"
source("network/experiments/netmet-060.R")
