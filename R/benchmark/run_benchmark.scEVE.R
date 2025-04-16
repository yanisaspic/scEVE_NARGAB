"Run this script to benchmark the default instance of the scEVE framework.

	2025/04/11 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
suppressPackageStartupMessages(library(feve))

dataset_label <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset_label)

params <- feve::get_parameters("scEVE")
benchmark <- feve::get_benchmark(data, params, "scEVE")
write.csv(benchmark, glue::glue("./results/benchmarks/{dataset_label}/scEVE.csv"), row.names=FALSE)