#!/usr/bin/env Rscript
# netmet-063  netmet-061 run over netmet-062's draws. A shim; the source
# table lives in netmet-061.R.
commandArgs <- function(trailingOnly = FALSE) "netmet-063"
source("network/experiments/netmet-061.R")
