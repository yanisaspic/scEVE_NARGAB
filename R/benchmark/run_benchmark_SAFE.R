"Run this script to benchmark the SAFE ensemble clustering method.

	2025/02/05 @yanisaspic"

source("./config/R_LIBS.R")
.libPaths(R_LIBS$SAFE)

source("./R/src/load_any_data.R")
source("./R/src/get_benchmark_method.R")
source("./R/src/methods/SAFE.R")

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
