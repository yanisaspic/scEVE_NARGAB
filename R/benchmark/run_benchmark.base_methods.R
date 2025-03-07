"Run this script to benchmark the four clustering methods integrated in scEVE, by default.

	2025/03/07 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
source("./R/src/get_benchmark_method.R")
source("./R/src/methods/base_methods.R")

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset)

base_methods <- c(
	"densityCut"=use_densityCut.benchmark,
	"monocle3"=use_monocle3.benchmark,
	"Seurat"=use_Seurat.benchmark,
	"SHARP"=use_SHARP.benchmark)

f <- function(method_label) {get_benchmark_method(data, base_methods[[method_label]], method_label)}
benchmarks <- lapply(X=names(base_methods), FUN=f)
benchmark <- do.call(rbind, benchmarks)
write.csv(benchmark, glue::glue("./results/benchmarks/{dataset}/base_methods.csv"), row.names=FALSE)