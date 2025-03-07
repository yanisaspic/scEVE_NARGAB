"Run this script to get the R requirements and their versions.

	2025/03/07 @yanisaspic"

source("./R/src/methods/base_methods.R")
source("./R/src/methods/ensemble_methods.R")
source("./R/src/plots/clustering_analyses.R")
source("./R/src/plots/real_benchmarks.R")
source("./R/src/plots/synthetic_benchmarks.R")
source("./R/src/plots/tables.R")
source("./R/src/plots/utils.R")
source("./R/src/plots/real_benchmarks.R")
source("./R/src/get_benchmark_method.R")
source("./R/src/get_scrnaseq_data.R")

suppressPackageStartupMessages({
  library(sessioninfo)})

sessioninfo::session_info(to_file="./config/session.log")