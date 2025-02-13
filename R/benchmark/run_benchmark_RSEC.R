"Run this script to benchmark the RSEC ensemble clustering method.

	2025/02/05 @yanisaspic"

source("./config/R_LIBS.R")

.libPaths(R_LIBS$scEVE)
source("./R/src/load_any_data.R")
source("./R/src/get_benchmark_method.R")

.libPaths(R_LIBS$RSEC)
source("./R/src/methods/RSEC.R")

print(packageVersion("TMExplorer"))
print(packageVersion("clusterExperiment"))

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
data <- load_any_data(dataset)

print(data)