"Run this script to benchmark the default instance of the scEVE framework.

	2025/02/14 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
suppressPackageStartupMessages(library(sceve))

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset)

params <- sceve::get_default_parameters()
benchmark <- get_benchmark_sceve.data(data, params, "scEVE")
write.csv(benchmark, glue::glue("./results/benchmarks/{dataset}/scEVE.csv"), row.names=FALSE)