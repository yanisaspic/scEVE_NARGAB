"Run this script to benchmark the scEFSC ensemble clustering method.

	2025/02/05 @yanisaspic"

source("./config/R_LIBS.R")
.libPaths(R_LIBS$scEFSC)

source("./R/src/load_any_data.R")
source("./R/src/get_benchmark_method.R")
source("./R/src/methods/scEFSC.R")

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
data <- load_any_data(dataset)
