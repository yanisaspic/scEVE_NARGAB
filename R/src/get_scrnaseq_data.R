"Function called to load scRNA-seq datasets.

	2025/02/13 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(SingleCellExperiment)
  library(SPARSim)
  library(SummarizedExperiment)
  library(TMExplorer)})

# _____________________________________________________________ real datasets

get_real_datasets_info <- function() {
  #' Get information regarding the real datasets available.
  #'
  #' @return a data.frame associating each dataset to its sequencing protocol,
  #' its number of cells, clusters and genes, as well as its year of publication,
  #' its accession number and its associated doi.
  #'
  characteristics <- c("sequencing_protocol", "n_cells", "n_clusters", "n_genes",
                       "year", "accession_number", "doi")
  metadata <- c("SMARTer (Fluidigm C1)", 364, 7, 57241, 2017, "GSE81861", "10.1038/ng.3818",
                "SMARTer (Fluidigm C1)", 561, 9, 55186, 2017, "GSE81861", "10.1038/ng.3818",
                "inDrop", 822, 13, 14878, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
                "inDrop", 1064, 13, 14878, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
                "inDrop", 1303, 14, 20125, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
                "SMARTer", 1679, 18, 24057, 2016, "GSE71585", "10.1038/nn.4216",
                "inDrop", 1724, 14, 20125, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
                "inDrop", 1937, 14, 20125, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
                "SMART-Seq2", 3589, 7, 23460, 2017, "GSE84465", "10.1016/j.celrep.2017.10.030",
                "inDrop", 3605, 14, 20125, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
                "SMART-Seq2", 6879, 9, 23686, 2018, "GSE115978", "10.1016/j.cell.2018.09.006",
                "10x Genomics", 18456, 18, 23580, 2020, "GSE125969", "10.1016/j.celrep.2020.108023",
                "Seq-Well", 22600, 17, 27899, 2018, "GSE116256", "10.1016/j.cell.2019.01.031",
                "10x Genomics", 51775, 17, 22180, 2018, "E-MTAB-6149,   E-MTAB-6653,", "10.1038/s41591-018-0096-5",
                "10x Genomics", 57530, 10, 24005, 2019, "CRA001160", "10.1038/s41422-019-0195-y")
  datasets <- c("Li_HumCRC_b", "Li_HumCRC_a", "Baron_MouPan_1", "Baron_MouPan_2", "Baron_HumPan_4",
                "Tasic_MouBra", "Baron_HumPan_2", "Baron_HumPan_1", "Darmanis_HumGBM", "Baron_HumPan_3",
                "JerbyArnon_HumMLM", "Gillen_HumEPN", "VanGalen_HumAML", "Lambrechts_HumNSCLC", "Peng_HumPDAC")
  output <- matrix(metadata, nrow=length(datasets), ncol=length(characteristics), byrow=TRUE)
  output <- as.data.frame.matrix(output, row.names=datasets)
  colnames(output) <- characteristics
  return(output)
}

get_cell_ids <- function(labels) {
  #' Get a vector of unique cell ids according to the label of each cell.
  #'
  #' @param labels a vector of cell types.
  #'
  #' @return a vector of characters.
  #'
  cell_ids <- c()
  counter <- list()
  for (cell_type in unique(labels)) {counter[[cell_type]] <- 1}
  for (n in 1:length(labels)) {
    cell_type <- labels[[n]]
    cell_ids[n] <- glue::glue("{cell_type}_{counter[[cell_type]]}")
    counter[[cell_type]] <- counter[[cell_type]] + 1}
  return(cell_ids)
}

get_real_data.TMExplorer <- function(dataset) {
  #' Get a data.frame and the ground truth of a dataset loaded with TMExplorer.
  #'
  #' @param dataset one of `Li_HumCRC_b`, `Darmanis_HumGBM`, `JerbyArnon_HumMLM`,
  #' `Gillen_HumEPN`, `VanGalen_HumAML`, `Lambrechts_HumNSCLC` or `Peng_HumPDAC`.
  #' (cf. `get_real_datasets_info()`)
  #'
  #' @return a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  real_datasets_info <- get_real_datasets_info()
  accession_number <- real_datasets_info[dataset, "accession_number"]
  data.TMExplorer <- TMExplorer::queryTME(accession_number)[[1]]
  
  if (accession_number == "GSE84465") {data.TMExplorer <- head(data.TMExplorer, -6)}
  # the 6 last rows are metadata in this dataset (Darmanis_HumGBM).

  throwaway_labels <- c("?", "NA", "", "Doublets")
  predictions <- SummarizedExperiment::colData(data.TMExplorer)
  ground_truth <- stats::setNames(predictions$label, rownames(predictions))
  ground_truth <- ground_truth[!ground_truth %in% throwaway_labels]
  expression.init <- SingleCellExperiment::counts(data.TMExplorer)
  expression.init <- expression.init[, names(ground_truth)]

  data <- list(expression.init=expression.init, ground_truth=ground_truth)
  return(data)
}

get_real_data.local <- function(dataset) {
  #' Get a data.frame and the ground truth of a dataset loaded locally.
  #'
  #' Run ./datasets/download_datasets.sh and ./datasets/setup_datasets.py first.
  #'
  #' @param dataset one of `Li_HumCRC_a`, `Baron_MouPan_1`, `Baron_MouPan_2`, `Baron_HumPan_4`,
  #' `Tasic_MouBra`, `Baron_HumPan_2`, `Baron_HumPan_1` or `Baron_HumPan_3`.
  #' (cf. `get_real_datasets_info()`)
  #'
  #' @return a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  expression.init <- read.csv(glue::glue("./datasets/scrna/{dataset}.csv"), row.names=1)
  get_label <- function(cell_id) {
    tmp <- strsplit(cell_id, split="_")[[1]]
    label <- paste(head(tmp, -1), collapse="_")
    return(label)}
  labels <- sapply(X=colnames(expression.init), FUN=get_label)
  ground_truth <- stats::setNames(labels, colnames(expression.init))
  data <- list(expression.init=expression.init, ground_truth=ground_truth)
  return(data)
}

# _____________________________________________________________ synthetic datasets

get_synthetic_datasets_info <- function() {
  #' Get information regarding the synthetic datasets available.
  #'
  #' @return a data.frame associating each dataset to its number of clusters, the balance of its
  #' cluster sizes, the independence of its cluster transcriptomes and its random state.
  #' 
  params <- list(k=1:10, balanced_clusters=c(TRUE, FALSE), independent_clusters=c(TRUE, FALSE), random_state=1:30)
  params <- expand.grid(params)
  return(params)
}

get_cluster_sizes <- function(k, balanced_clusters, dataset_size) {
  #' Get cluster sizes for a synthetic dataset.
  #' 
  #' @param k the number of clusters.
  #' @param balanced_clusters a boolean indicating if the cluster sizes are balanced.
  #' @param dataset_size the number of cells in the dataset.
  #' 
  #' @return a vector of sizes.
  #' 
  if (balanced_clusters) {cluster_sizes <- rep(dataset_size / k, k)}
  else {
    get_size <- function(n) {
      a <- 1/(2^n) * dataset_size
      b <- 2^(k-n) / (2^k - 1)
      c <- 1/(2^k) * dataset_size
      size <- a + b * c
      return(size)}
    n <- 1:k
    cluster_sizes <- sapply(X=n, FUN=get_size)}
  cluster_sizes <- as.integer(cluster_sizes)
  return(cluster_sizes)
}

get_cluster_foldchanges <- function(cluster_params) {
  #' Get the foldchanges values of a cluster.
  #' 
  #' @param cluster_params a named list, with three names: `intensity`, `variability` and `lib_size`.
  #' 
  #' @return a vector of foldchanges.
  #' 
  n_genes <- length(cluster_params$variability)
  negative_foldchanges <- runif(n=250, min=0.0001, max=0.25)
  positive_foldchanges <- runif(n=250, min=4, max=100)
  null_foldchanges <- rep(1, n_genes-500)
  foldchanges <- c(negative_foldchanges, positive_foldchanges, null_foldchanges)
  foldchanges <- sample(foldchanges)
  return(foldchanges)
}

get_cluster_params <- function(cluster_label, cluster_size, init_params) {
  #' Get a list of parameters specific to a cluster.
  #' 
  #' @param cluster_label a character.
  #' @param cluster_size an integer.
  #' @param init_params a named list, with three names: `intensity`, `variability` and `lib_size`.
  #' 
  #' @return a named list, with three names: `intensity`, `variability` and `lib_size`.
  #' 
  foldchanges <- get_cluster_foldchanges(init_params)
  cluster_params <- SPARSim::SPARSim_create_DE_genes_parameter(init_params, foldchanges, cluster_size,
                                                                   condition_name=cluster_label)
  cluster_params$variability <- cluster_params$variability * runif(n=length(cluster_params$variability),
                                                                   min=1, max=1.5)
  return(cluster_params)
}

get_related_clusters_params <- function(cluster_sizes, init_params) {
  #' Get lists of parameters specific to every cluster in a related scRNA-seq dataset.
  #' Related scRNA-seq clusters are generated by adding DEGs to the last cluster generated, recursively.
  #' 
  #' @param cluster_sizes a vector of integers.
  #' @param init_params a named list, with three names: `intensity`, `variability` and `lib_size`.
  #' 
  #' @return a named list of elements. Every element is a named list, with three names: 
  #' `intensity`, `variability` and `lib_size`.
  #' 
  related_cluster_params <- list()
  cluster_params.previous <- init_params
  for (i in 1:length(cluster_sizes)) {
    related_cluster_params[[i]] <- get_cluster_params(as.character(i), cluster_sizes[i], cluster_params.previous)
    cluster_params.previous <- related_cluster_params[[i]]}
  return(related_cluster_params)
}

get_independent_clusters_params <- function(cluster_sizes, init_params) {
  #' Get lists of parameters specific to every cluster in an independent scRNA-seq dataset.
  #' Independent scRNA-seq clusters are generated by adding DEGs to a base cluster, parallely.
  #' 
  #' @param cluster_sizes a vector of integers.
  #' @param init_params a named list, with three names: `intensity`, `variability` and `lib_size`.
  #' 
  #' @return a named list of elements. Every element is a named list, with three names: 
  #' `intensity`, `variability` and `lib_size`.
  #' 
  f <- function(i) {get_cluster_params(as.character(i), cluster_sizes[i], init_params)}
  independent_cluster_params <- lapply(X=1:length(cluster_sizes), FUN=f)
  return(independent_cluster_params)
}

simulate_scrnaseq_data <- function(k, balanced_clusters, independent_clusters, random_state,
                                   init_params, dataset_size) {
  #' Simulate a scRNA-seq dataset using `SPARSim`.
  #' 
  #' @param k the number of clusters.
  #' @param balanced_clusters a boolean indicating if the cluster sizes are balanced.
  #' @param independent_clusters a boolean indicating if the cluster transcriptomes are independent.
  #' @param random_state a numeric.
  #' @param init_params a named list, with three names: `intensity`, `variability` and `lib_size`.
  #' @param dataset_size the number of cells in the dataset.
  #' 
  #' @return a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  set.seed(random_state)
  cluster_sizes <- get_cluster_sizes(k, balanced_clusters, dataset_size)
  if (independent_clusters) {cluster_params <- get_independent_clusters_params(cluster_sizes, init_params)}
  else {cluster_params <- get_related_clusters_params(cluster_sizes, init_params)}
  
  tmp <- SPARSim::SPARSim_simulation(cluster_params)
  expression.init <- as.data.frame(tmp$count_matrix)
  expression.init <- expression.init[, sample(colnames(expression.init))]
  return(expression.init)
}

get_synthetic_data <- function(dataset) {
  #' Load a synthetic scRNA-seq dataset of raw count expression.
  #'
  #' @param dataset a numeric ranging from 1 to 360. Each numeric is associated to a synthetic
  #' dataset with specific properties.
  #'
  #' @return a named list with two elements: `expression.init` and `ground_truth`.
  #' `expression.init` is a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' `ground_truth` is a named factor associating cells to their cluster annotations.
  #'
  data("PBMC_10X_param_preset")
  init_params <- PBMC_10X_param_preset[[1]]
  synthetic_datasets_info <- get_synthetic_datasets_info()
  tmp <- synthetic_datasets_info[dataset, ]
  
  expression.init <- simulate_scrnaseq_data(tmp$k, tmp$balanced_clusters, tmp$independent_clusters, tmp$random_state,
                                            init_params, dataset_size=10000)
  get_cell_type <- function(cell_id) {strsplit(cell_id, split="_")[[1]][1]}
  ground_truth <- sapply(X=colnames(expression.init), FUN=get_cell_type)
  data <- list(expression.init=expression.init, ground_truth=ground_truth)
  return(data)
}

# ___________________________________________________________________________________

get_scrnaseq_data <- function(dataset) {
  #' Load a scRNA-seq dataset of raw count expression.
  #'
  #' @param dataset either the name of a real scRNA-seq dataset, or a numeric ranging from 1 to 600 (synthetic data).
  #'
  #' @return a named list with two elements: `expression.init` and `ground_truth`.
  #' `expression.init` is a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' `ground_truth` is a named factor associating cells to their cluster annotations.
  #'
  if (dataset %in% c("Li_HumCRC_b", "Darmanis_HumGBM", "JerbyArnon_HumMLM", "Gillen_HumEPN",
                     "VanGalen_HumAML", "Lambrechts_HumNSCLC", "Peng_HumPDAC")) {data <- get_real_data.TMExplorer(dataset)}
  else if (dataset %in% c("Li_HumCRC_a", "Baron_MouPan_1", "Baron_MouPan_2", "Baron_HumPan_4", "Tasic_MouBra",
                     "Baron_HumPan_2", "Baron_HumPan_1", "Baron_HumPan_3")) {data <- get_real_data.local(dataset)}
  else {data <- get_synthetic_data(dataset)}
  
  rownames(data$expression.init) <- gsub("_", "+", rownames(data$expression.init))
  cell_ids <- get_cell_ids(data$ground_truth)
  names(data$ground_truth) <- cell_ids
  colnames(data$expression.init) <- cell_ids
  return(data)
}