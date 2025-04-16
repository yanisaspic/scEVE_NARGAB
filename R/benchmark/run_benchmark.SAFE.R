"Run this script to benchmark SAFE.

	2025/04/11 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
source("./R/src/get_benchmark_method.R")
source("./R/src/methods/ensemble_methods.R")

dataset_label <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset_label)

benchmark <- get_benchmark_method(data, use_SAFE.benchmark, "SAFE")
write.csv(benchmark, glue::glue("./results/benchmarks/{dataset_label}/SAFE.csv"), row.names=FALSE)