"Functions called to benchmark clustering methods.

	2025/04/11 @yanisaspic"

suppressPackageStartupMessages({
  library(aricode)
  library(bluster)
  library(scater)
  library(scran)
  library(scuttle)
  library(SingleCellExperiment)})

get_data.bluster <- function(dataset) {
  #' Get a matrix usable for the functions of bluster.
  #'
  #' @param dataset a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @return an object associating selected genes and cells to their reduced dimensions.
  #'
  #' @import scater
  #' @import scran
  #' @import scuttle
  #' @import SingleCellExperiment
  #'
  data <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as.matrix(dataset)))
  data <- scuttle::logNormCounts(data)
  variable_genes <- scran::getTopHVGs(scran::modelGeneVar(data), n=1000)
  set.seed(1)
  data <- scater::runPCA(data, ncomponents=20, subset_row=variable_genes)
  data <- SingleCellExperiment::reducedDim(data)
  return(data)
}

get_clustering_metrics.intrinsic <- function(dataset, preds) {
  #' Using intrinsic metrics, measure a clustering performance.
  #'
  #' @param dataset a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param preds a named factor associating cells to their predicted clusters.
  #'
  #' @return a named vector of two: `nPurity` and `SI`.
  #'
  #' @import bluster
  #'
  n_clusters_predicted <- length(unique(preds))
  if (n_clusters_predicted < 2) {return(c("nPurity"=NA, "SI"=NA))}
  data <- get_data.bluster(dataset)
  neighborhood_purity <- bluster::neighborPurity(data, preds)
  silhouette_index <- bluster::approxSilhouette(data, preds)
  clustering_metrics.intrinsic <- c("nPurity"=mean(neighborhood_purity$purity),
                                    "SI"=mean(silhouette_index$width))
  return(clustering_metrics.intrinsic)
}

get_clustering_metrics <- function(data, preds) {
  #' Using both intrinsic and extrinsic clustering metrics, measure a clustering performance.
  #' Extrinsic metrics compare cluster predictions to the cell annotations of the dataset.
  #' Intrinsic metrics compare the gene expression of cells in and out of their clusters.
  #' The `ARI` and the `NMI` are extrinsic metrics.
  #' The `nPurity` and the `SI` are intrinsic metrics.
  #' For every metric, higher is better and the maximum value is 1.
  #'
  #' @param data a named list with two elements: `dataset` and `ground_truth`.
  #' `dataset` is a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' `ground_truth` is a named factor associating cells to their cluster annotations.
  #' @param preds a named factor associating cells to their predicted clusters.
  #'
  #' @return a named vector with four names: `ARI`, `NMI`, `nPurity` and `SI`.
  #'
  #' @import aricode
  #'
  if (length(preds) < 2) {return(c("ARI"=NA, "NMI"=NA, "nPurity"=NA, "SI"=NA))}
  ground_truth <- data$ground_truth[names(preds)]
  dataset <- data$dataset[, names(preds)]
  clustering_metrics <- c("ARI"=aricode::ARI(ground_truth, preds),
                          "NMI"=aricode::NMI(ground_truth, preds),
                          get_clustering_metrics.intrinsic(dataset, preds))
  return(clustering_metrics)
}

get_benchmark_method <- function(data, clustering_method, method_label, random_state=1) {
  #' Using computational as well as intrinsic and extrinsic clustering metrics, measure
  #' the performance of a clustering method on a dataset.
  #'
  #' @param data a named list with two elements: `dataset` and `ground_truth`.
  #' `dataset` is a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' `ground_truth` is a named factor associating cells to their cluster annotations.
  #' @param clustering_method a function taking `dataset` as input, and outputing
  #' a factor associating cells to their predicted clusters.
  #' @param method_label a character.
  #' @param random_state a numeric.
  #'
  #' @return a data.frame with seven columns: `method`, `time (s)`,
  #' `peak_memory_usage (Mb)`, `ARI`, `NMI`, `nPurity` and `SI`.
  #'
  get_memory_usage <- function(memory) {memory[[11]] + memory[[12]]}
  memory_usage.init <- get_memory_usage(gc(reset=TRUE))
  time.init <- Sys.time()
  preds <- clustering_method(data$dataset, random_state)
  time <- as.numeric(Sys.time() - time.init, units="secs")
  peak_memory_usage <- get_memory_usage(gc()) - memory_usage.init
  benchmark <- c("method"=method_label, "time (s)"=time,
                 "peak_memory_usage (Mb)"=peak_memory_usage,
                 get_clustering_metrics(data, preds))
  benchmark <- as.data.frame(t(benchmark))
  return(benchmark)
}