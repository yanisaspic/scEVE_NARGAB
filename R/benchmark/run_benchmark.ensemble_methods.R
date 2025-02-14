"Run this script to benchmark four ensemble clustering methods available in R packages.

	2025/02/14 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
source("./R/src/get_benchmark_method.R")
source("./R/src/methods/ensemble_methods.R")

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset)

ensemble_methods <- c(
	"RSEC"=use_RSEC.benchmark,
	"SAFE"=use_SAFE.benchmark,
	"SAME"=use_SAME.benchmark,
	"scEFSC"=use_scEFSC.benchmark)

# f <- function(method_label) {get_benchmark_method(data, base_methods[[method_label]], method_label)}
# benchmarks <- lapply(X=names(ensemble_methods), FUN=f)
# benchmark <- do.call(rbind, benchmarks)
# write.csv(benchmark, glue::glue("./results/benchmarks/{dataset}/ensemble_methods.csv"), row.names=TRUE)