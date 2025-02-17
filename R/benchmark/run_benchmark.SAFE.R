"Run this script to benchmark SAFE.

	2025/02/16 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
source("./R/src/get_benchmark_method.R")
source("./R/src/methods/ensemble_methods.R")

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset)

benchmark <- get_benchmark_method(data, use_SAFE.benchmark, "SAFE")
write.csv(benchmark, glue::glue("./results/benchmarks/{dataset}/SAFE.csv"), row.names=TRUE)