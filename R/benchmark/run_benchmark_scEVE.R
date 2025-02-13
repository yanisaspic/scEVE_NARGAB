"Run this script to benchmark the default instance of the scEVE framework.

	2025/02/05 @yanisaspic"

source("./config/R_LIBS.R")
.libPaths(R_LIBS$scEVE)
suppressPackageStartupMessages(library("sceve"))
.libPaths(R_LIBS$utils)
source("./R/src/load_data.R")

dataset <- commandArgs(trailingOnly=TRUE)[[1]]

params <- sceve::get_default_parameters()
benchmark <- get_benchmark_sceve(c(dataset), params, "scEVE")
write.csv(benchmark, glue::glue("./results/benchmarks/{dataset}/scEVE.csv"), row.names=FALSE)