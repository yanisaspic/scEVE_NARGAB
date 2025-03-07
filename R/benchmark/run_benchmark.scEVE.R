"Run this script to benchmark the default instance of the scEVE framework.

	2025/03/07 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
suppressPackageStartupMessages(library(feve))

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset)

params <- feve::get_parameters("scEVE")
benchmark <- feve::get_benchmark(data, params, "scEVE")
write.csv(benchmark, glue::glue("./results/benchmarks/{dataset}/scEVE.csv"), row.names=FALSE)