"Functions called by multiple scripts related to plots.

	2025/03/07 @yanisaspic"

suppressPackageStartupMessages(library(feve))

get_analyses <- function(path) {
  #' Load every clustering analysis conducted in a single list.
  #' 
  #' @param path a path where benchmark results are stored.
  #'
  #' @return a named list. Each name corresponds to a dataset analysed.
  #'
  get_dataset <- function(file) {
    tmp <- strsplit(file, split="/", fixed=TRUE)[[1]]
    dataset <- tmp[length(tmp)-1]
    return(dataset)}
  has_clusters <- function(analysis) {nrow(analysis$meta) > 1}
  
  files <- list.files(path, pattern=".xlsx", recursive=TRUE, full.names=TRUE)
  analyses <- lapply(X=files, FUN=feve::get_records)
  names(analyses) <- sapply(X=files, FUN=get_dataset)
  analyses <- Filter(f=has_clusters, x=analyses)
  return(analyses)
}

get_benchmarks <- function(path) {
  #' Load every benchmark conducted in a single data.frame.
  #' 
  #' @param path a path where benchmark results are stored.
  #'
  #' @return a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #'
  get_benchmark.dataset <- function(file) {
    tmp <- strsplit(file, split="/", fixed=TRUE)[[1]]
    dataset <- tmp[length(tmp)-1]
    benchmark.dataset <- read.csv(file)
    benchmark.dataset[, "dataset"] <- dataset
    return(benchmark.dataset)}
  
  files <- list.files(path, pattern=".csv", recursive=TRUE, full.names=TRUE)
  benchmarks <- lapply(X=files, FUN=get_benchmark.dataset)
  benchmark <- dplyr::bind_rows(benchmarks)
  
  for (col in c("time..s.", "peak_memory_usage..Mb.")) {benchmark[, col] <- log10(benchmark[, col])}
  colnames(benchmark)[colnames(benchmark)=="time..s."] <- "log10(s)"
  colnames(benchmark)[colnames(benchmark)=="peak_memory_usage..Mb."] <- "log10(Mbytes)"
  colnames(benchmark)[colnames(benchmark)=="Purity"] <- "nPurity"
  
  aes_params <- get_aesthetic_parameters()
  benchmark$method <- factor(benchmark$method, levels=names(aes_params$method_types))
  is_real <- function(dataset) {dataset %in% aes_params$sorted_datasets}
  benchmark[, "is_real"] <- sapply(X=benchmark$dataset, FUN=is_real)
  return(benchmark)
}

get_aesthetic_parameters <- function() {
  #' Get parameters required to draw the figures effectively.
  #' 
  #' @return a named list of parameters, with three names: `method_types`, `colormap` and `sorted_datasets`.
  #' 
  method_types <- c("densityCut"="base_method", "monocle3"="base_method", "Seurat"="base_method", "SHARP"="base_method",
                    "SAFE"="ensemble_method", "SAME"="ensemble_method",
                    "RSEC"="main", "RSEC*"="main", "scEVE"="main", "scEVE*"="main",
                    "ground_truth"="secondary")
  colormap <- list("densityCut"="#56B4E9", "monocle3"="#009E73",
                   "Seurat"="#D55E00", "SHARP"="#CC79A7",
                   "SAFE"="#F0E442", "SAME"="#E69F00",
                   "RSEC"="#666666", "RSEC*"="#666666",
                   "scEVE"="#0072B2", "scEVE*"="#0072B2",
                   "ground_truth"="#FFFFFF")
  sorted_datasets <- c("Li_HumCRC_b", "Li_HumCRC_a", "Baron_MouPan_1", "Baron_MouPan_2",
                       "Baron_HumPan_4", "Tasic_MouBra", "Baron_HumPan_2", "Baron_HumPan_1",
                       "Darmanis_HumGBM", "Baron_HumPan_3", "JerbyArnon_HumMLM", "Gillen_HumEPN",
                       "VanGalen_HumAML", "Lambrechts_HumNSCLC", "Peng_HumPDAC")
  aesthetic_parameters <- list(method_types=method_types, colormap=colormap, sorted_datasets=sorted_datasets)
  return(aesthetic_parameters)
}