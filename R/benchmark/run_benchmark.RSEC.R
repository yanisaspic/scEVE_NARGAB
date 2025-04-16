"Run this script to benchmark RSEC.

    2025/04/11 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
source("./R/src/get_benchmark_method.R")
source("./R/src/methods/ensemble_methods.R")

dataset_label <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset_label)

get_benchmark_rsec <- function(data, random_state=1) {
  #' Using computational as well as intrinsic and extrinsic clustering metrics, measure
  #' the performance of RSEC on a dataset.
  #'
  #' @param data a named list with two elements: `dataset` and `ground_truth`.
  #' `dataset` is a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' `ground_truth` is a named factor associating cells to their cluster annotations.
  #' @param random_state a numeric.
  #'
  #' @return a data.frame with eight columns: `method`, `time (s)`,
  #' `peak_memory_usage (Mb)`, `ARI`, `NMI`, `nPurity`, `SI` and `n_samples`.
  #'
  get_memory_usage <- function(memory) {memory[[11]] + memory[[12]]}
  memory_usage.init <- get_memory_usage(gc(reset=TRUE))
  time.init <- Sys.time()
  predictions <- use_RSEC.benchmark(data$dataset, random_state)
  time <- as.numeric(Sys.time() - time.init, units="secs")
  peak_memory_usage <- get_memory_usage(gc()) - memory_usage.init
  
  benchmark <- c("method"="RSEC", "time (s)"=time, "peak_memory_usage (Mb)"=peak_memory_usage,
                 get_clustering_metrics(data, predictions), "n_samples"=length(predictions))
  benchmark.star <- c("method"="RSEC*", "time (s)"=time, "peak_memory_usage (Mb)"=peak_memory_usage,
                      get_clustering_metrics(data, predictions[predictions != -1]),
                      "n_samples"=length(predictions[predictions != -1]))
  
  benchmark_rsec <- as.data.frame(rbind(benchmark, benchmark.star))
  return(benchmark_rsec)
}

benchmark <- get_benchmark_rsec(data)
write.csv(benchmark, glue::glue("./results/benchmarks/{dataset_label}/RSEC.csv"), row.names=FALSE)