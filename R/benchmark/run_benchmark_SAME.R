"Run this script to benchmark the SAME ensemble clustering method.

	2025/02/05 @yanisaspic"

source("./config/R_LIBS.R")

.libPaths(R_LIBS$scEVE)
source("./R/src/load_any_data.R")
source("./R/src/get_benchmark_method.R")

.libPaths(R_LIBS$SAME)
source("./R/src/methods/SAME.R")

# dataset <- commandArgs(trailingOnly=TRUE)[[1]]
# data <- load_any_data(dataset)

data("data_SAME")

#test <- SC3::calculate_distance(data_SAME$Zheng.expr)
test <- use_SAME(data_SAME$Zheng.expr)
print(test)