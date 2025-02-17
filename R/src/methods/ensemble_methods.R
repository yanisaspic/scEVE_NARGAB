suppressPackageStartupMessages({
  library(clusterExperiment)
  library(SAFEclustering)
  library(SAMEclustering)})

do_SC3.SAxE <- function(expression.init, random_state) {
  #' Predict clusters with the SC3 method, for the SAFE and SAME algorithms.
  #'
  #' This function is directly copied from the repositories of the SAFE and SAME packages.
  #' c.f. https://github.com/yycunc/SAFEclustering/blob/master/R/individual_clustering.R
  #' https://github.com/yycunc/SAMEclustering/blob/master/R/SAMEclustering.R
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #' 
  #' @return a named list (`base_clusters`), with a named vector associating cells to their predicted cluster.
  #'
  
  # default hyperparameters
  gene_filter <- FALSE
  svm_num_cells <- 5000
  
  # data pre-processing
  data <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = expression.init))
  SingleCellExperiment::normcounts(data) <- t(t(expression.init)/colSums(expression.init)) * 1000000
  SingleCellExperiment::logcounts(data) <- log2(SingleCellExperiment::normcounts(data) + 1)
  SummarizedExperiment::rowData(data)$feature_symbol <- rownames(data)
  data <- data[!duplicated(SummarizedExperiment::rowData(data)$feature_symbol), ]
  
  # prediction of k clusters
  data <- SC3::sc3_estimate_k(data)
  optimal_K <- S4Vectors::metadata(data)$sc3$k_estimation
  if (ncol(expression.init) < svm_num_cells) {
    data <- SC3::sc3(data, ks=optimal_K, biology=FALSE, gene_filter=gene_filter, n_cores=1,
                     rand_seed=random_state)}
  else {
    data <- SC3::sc3(data, ks=optimal_K, biology=FALSE, gene_filter=gene_filter, svm_max = svm_num_cells,
                     svm_num_cells=svm_num_cells, n_cores=1, rand_seed=random_state)
    data <- SC3::sc3_run_svm(exp_cell_exprs, ks=optimal_K)}
  
  # report results
  tmp <- SummarizedExperiment::colData(data)
  column_name <- paste("sc3_", optimal_K, "_clusters", sep = '')
  base_clusters <- stats::setNames(as.numeric(tmp[, column_name]), rownames(tmp))
  output <- list(base_clusters=base_clusters)
  return(output)
}

do_CIDR.SAxE <- function(expression.init, random_state) {
  #' Predict clusters with the CIDR method, for the SAFE and SAME algorithms.
  #'
  #' This function is directly copied from the repositories of the SAFE and SAME packages.
  #' c.f. https://github.com/yycunc/SAFEclustering/blob/master/R/individual_clustering.R
  #' https://github.com/yycunc/SAMEclustering/blob/master/R/SAMEclustering.R
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #' 
  #' @return a named list (`base_clusters`, `nPC`), with a named vector associating cells to their predicted cluster.
  #' 
  set.seed(random_state)
  
  # default hyperparameters
  percent_dropout <- 10
  
  # data pre-processing
  dropouts <- rowSums(expression.init == 0) / ncol(expression.init) * 100
  expression.init.cidr <- expression.init[-c(which(dropouts <= percent_dropout), which(dropouts >= 100-percent_dropout)),]
  data <- cidr::scDataConstructor(expression.init.cidr, tagType = "raw")
  data <- cidr::determineDropoutCandidates(data)
  data <- cidr::wThreshold(data)
  data <- cidr::scDissim(data)
  data <- cidr::scPCA(data, plotPC=FALSE)
  data <- cidr::nPC(data)
  
  # prediction of clusters
  nPC.cidr <- data@nPC
  data <- cidr::scCluster(data, nPC = nPC.cidr)
  
  # report results
  base_clusters <- stats::setNames(data@clusters, colnames(expression.init))
  output <- list(base_clusters=base_clusters, nPC=nPC.cidr)
  return(output)
}

do_Seurat.SAxE <- function(expression.init, random_state, nPC) {
  #' Predict clusters with the Seurat method, for the SAFE and SAME algorithms.
  #'
  #' This function is directly copied from the repositories of the SAFE and SAME packages.
  #' c.f. https://github.com/yycunc/SAFEclustering/blob/master/R/individual_clustering.R
  #' https://github.com/yycunc/SAMEclustering/blob/master/R/SAMEclustering.R
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #' @param nPC a numeric.
  #' 
  #' @return a named list (`base_clusters`), with a named vector associating cells to their predicted cluster.
  #' 
  
  # default hyperparameters
  high.genes <- 8000
  low.genes <- 200
  nGene_filter <- TRUE
  resolution <- 0.7
  
  # data pre-processing
  data <- Seurat::CreateSeuratObject(counts=expression.init, min.cells=0, min.features=0, project="single-cell clustering")
  if (nGene_filter == TRUE) {data <- subset(data, subset = nFeature_RNA > low.genes & nFeature_RNA < high.genes)}
  data <- Seurat::NormalizeData(data, normalization.method="LogNormalize", scale.factor=10000)
  data <- Seurat::FindVariableFeatures(data, selection.method="vst", nfeatures=2000)
  all.genes <- rownames(data)
  data <- Seurat::ScaleData(data, features=all.genes)
  data <- Seurat::RunPCA(data, features=Seurat::VariableFeatures(data), npcs=max(nPC, 20), seed.use=random_state)
  data <- Seurat::FindNeighbors(data, dims = 1:max(nPC, 20))
  
  # prediction of clusters
  data <- Seurat::FindClusters(data, resolution=resolution, random.seed=random_state)
  if (length(data@active.ident) < ncol(expression.init)) {
    base_clusters <- matrix(NA, ncol=ncol(expression.init), byrow=T)
    colnames(base_clusters) <- colnames(expression.init)
    tmp <- t(as.matrix(as.numeric(data@active.ident)))
    colnames(tmp) <- names(data@active.ident)
    for (i in 1:ncol(tmp)) {
      base_clusters[1, colnames(tmp)[i]] <- tmp[1, colnames(tmp)[i]]}}
  else {base_clusters <- t(as.matrix(as.numeric(data@active.ident)))}
  
  # report results
  base_clusters <- stats::setNames(as.vector(base_clusters), colnames(expression.init))
  output <- list(base_clusters=base_clusters)
  return(output)
}

do_tSNE_kMeans.SAxE <- function(expression.init, random_state, k.max) {
  #' Predict clusters with the t-SNE + k-means method, for the SAFE and SAME algorithms.
  #'
  #' This function is directly copied from the repositories of the SAFE and SAME packages.
  #' c.f. https://github.com/yycunc/SAFEclustering/blob/master/R/individual_clustering.R
  #' https://github.com/yycunc/SAMEclustering/blob/master/R/SAMEclustering.R
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #' @param k.max a numeric.
  #' 
  #' @return a named list (`base_clusters`), with a named vector associating cells to their predicted cluster.
  #' 
  set.seed(random_state)
  
  # default hyperparameters
  dimensions <- 3
  k.min <- 2
  percent_dropout <- 10
  tsne_min_cells <- 200
  tsne_min_perplexity <- 10
  perplexity <- ifelse(ncol(expression.init) < tsne_min_cells, tsne_min_perplexity, 30)
  
  # data pre-processing
  dropouts <- rowSums(expression.init == 0) / ncol(expression.init) * 100
  data <- expression.init[-c(which(dropouts <= percent_dropout), which(dropouts >= 100 - percent_dropout)),]
  data <- log2(t(t(data)/colSums(data)) * 1000000 + 1)
  
  # prediction of k clusters
  output.tsne <- Rtsne::Rtsne(t(data), dims=dimensions, perplexity=perplexity, check_duplicates=FALSE)
  tmp <- ADPclust::adpclust(output.tsne$Y, htype="amise", centroids="auto", nclust=k.min:k.max)
  output.kmeans <- stats::kmeans(output.tsne$Y, output.tsne$Y[tmp$centers, ], tmp$nclust)
  
  # report results
  base_clusters <- stats::setNames(as.vector(output.kmeans$cluster), colnames(expression.init))
  output <- list(base_clusters=base_clusters)
  return(output)
}

do_SIMLR.SAME <- function(expression.init, random_state, k.max) {
  #' Predict clusters with the SIMLR method, for the SAME algorithm.
  #'
  #' This function is directly copied from the repositories of the SAFE and SAME packages.
  #' c.f. https://github.com/yycunc/SAMEclustering/blob/master/R/SAMEclustering.R
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #' @param k.max a numeric.
  #' 
  #' @return a named list (`base_clusters`), with a named vector associating cells to their predicted cluster.
  #' 
  set.seed(random_state)
  
  # default hyperparameters
  k.min <- 2
  percent_dropout <- 10
  
  # data pre-processing
  dropouts <- rowSums(expression.init == 0) / ncol(expression.init) * 100
  data <- expression.init[-c(which(dropouts <= percent_dropout), which(dropouts >= 100 - percent_dropout)),]
  tmp <- data <- log2(t(t(data)/colSums(data)) * 1000000 + 1)
  data <- log10(data + 1)
  k.range <- k.min:k.max
  
  # prediction of k clusters
  output.k <- SIMLR::SIMLR_Estimate_Number_of_Clusters(tmp, NUMC=k.range, cores.ratio=1)
  k <- which.min(output.k$K1) + k.range[1] - 1
  if (ncol(expression.init) < 1000) {results <- SIMLR::SIMLR(data, c=k, cores.ratio=0)}
  else {results <- SIMLR::SIMLR_Large_Scale(data, c=k)}
  
  # report the results
  base_clusters <- stats::setNames(results$y$cluster, colnames(expression.init))
  output <- list(base_clusters=base_clusters)
  return(output)
}

get_base_clusters.SAxE <- function(expression.init, random_state, clustering_methods, diverse) {
  #' Predict a set of 4 base clusters with multiple methods.
  #'
  #' This function is directly copied from the repositories of the SAFE and SAME packages.
  #' c.f. https://github.com/yycunc/SAFEclustering/blob/master/R/individual_clustering.R
  #' https://github.com/yycunc/SAMEclustering/blob/master/R/SAMEclustering.R
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #' @param clustering_methods a vector of valid clustering names: `CIDR`, `SC3`, `Seurat` `SIMR` or `tSNE_kMeans`.
  #' @param diverse a boolean that indicates if the most divergent base clusters should be selected (`SAME` algorithm).
  #' 
  #' @return a named list (`base_clusters`), with a named vector associating cells to their predicted cluster.
  #' 
  k.max <- 2
  base_clusters <- list()
  
  # sequential predictions of base clusters_____________________________________
  if ("SC3" %in% clustering_methods) {
    base_clusters[["SC3"]] <- do_SC3.SAxE(expression.init, random_state)$base_clusters
    k.max <- max(k.max, max(base_clusters$SC3))}
  if ("CIDR" %in% clustering_methods) {
    output <- do_CIDR.SAxE(expression.init, random_state)
    nPC <- output$nPC
    base_clusters[["CIDR"]] <- output$base_clusters
    k.max <- max(k.max, max(base_clusters$CIDR))}
  if ("Seurat" %in% clustering_methods) {
    base_clusters[["Seurat"]] <- do_Seurat.SAxE(expression.init, random_state, nPC)$base_clusters
    k.max <- max(k.max, max(base_clusters$Seurat, na.rm=TRUE))}
  if ("tSNE_kMeans" %in% clustering_methods) {
    output <- do_tSNE_kMeans.SAxE(expression.init, random_state, k.max)
    base_clusters[["tSNE_kMeans"]] <- output$base_clusters
    k.max <- max(k.max, max(base_clusters$tSNE_kMeans))}
  if ("SIMLR" %in% clustering_methods) {
    base_clusters[["SIMLR"]] <- do_SIMLR.SAME(expression.init, random_state, k.max)$base_clusters}
  
  output <- do.call(rbind, base_clusters)
  if (diverse) {
    n <- nrow(output)
    similarity <- matrix(0, nrow=n, ncol=n, dimnames=list(clustering_methods, clustering_methods))
    
    for (row in clustering_methods) {
      for (col in clustering_methods) {
        similarity[row, col] <- cidr::adjustedRandIndex(output[row, ], output[col, ])}}
    
    least_diverse_method <- which.min(apply(X=similarity, MARGIN=1, FUN=var))
    output <- output[-least_diverse_method, ]}
  return(output)
}

get_formatted_predictions <- function(cells, predictions) {
  #' Format the cluster predictions of a method to obtain a standard output.
  #'
  #' @param cells a vector of cells names.
  #' @param predictions a vector of cluster predictions.
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  predictions <- stats::setNames(predictions, cells)
  predictions <- factor(predictions)
  return(predictions)
}

use_SAME.benchmark <- function(expression.init, random_state) {
  #' Predict clusters with the SAME algorithm.
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #' 
  #' @return a named factor that associates each cell to their cluster prediction.
  #' 
  expression.init <- as.matrix(expression.init)
  base_clusters <- get_base_clusters.SAxE(expression.init, random_state,
                                          clustering_methods=c("CIDR", "SC3", "Seurat", "SIMLR", "tSNE_kMeans"),
                                          diverse=TRUE)
  ensemble_clusters <- SAMEclustering::SAMEclustering(Y=as.matrix(t(base_clusters)), rep=3, SEED=random_state)
  predictions <- get_formatted_predictions(colnames(expression.init), ensemble_clusters$BICcluster)
  return(predictions)
}

use_SAFE.benchmark <- function(expression.init, random_state, program.dir="./config/dependencies") {
  #' Predict clusters with the SAFE algorithm.
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #' @param program.dir a path where executable gpmetis and shmetis programs are stored.
  #' (cf. https://github.com/yycunc/SAFEclustering)
  #' 
  #' @return a named factor that associates each cell to their cluster prediction.
  #' 
  expression.init <- as.matrix(expression.init)
  base_clusters <- get_base_clusters.SAxE(expression.init, random_state,
                                          clustering_methods=c("CIDR", "SC3", "Seurat", "tSNE_kMeans"),
                                          diverse=FALSE)
  ensemble_clusters <- SAFEclustering::SAFE(as.matrix(base_clusters), program.dir=program.dir,
                                            MCLA=TRUE, CSPA=TRUE, HGPA=TRUE, SEED=random_state)
  predictions <- get_formatted_predictions(colnames(expression.init), ensemble_clusters$optimal_clustering)
  return(predictions)
}

use_RSEC.benchmark <- function(expression.init, random_state) {
  #' Predict clusters with the RSEC algorithm.
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #' 
  #' @return a named factor that associates each cell to their cluster prediction.
  #' 
  expression.init <- as.matrix(expression.init)
  
  # data pre-processing
  data <- SummarizedExperiment::SummarizedExperiment(expression.init)
  dropout_genes <- which(rowSums(SummarizedExperiment::assay(data))==0)
  tmp <- apply(SummarizedExperiment::assay(data), 1, function(x) {length(x[x >= 10]) >= 10})
  data <- data[tmp, ]
  
  # prediction of k clusters
  output <- clusterExperiment::RSEC(data, isCount=TRUE, random.seed=random_state)
  predictions <- get_formatted_predictions(colnames(expression.init), primaryCluster(output))
  return(predictions)
}

use_scEFSC.benchmark <- function(expression.init, random_state) {
  #' Predict clusters with the scEFSC algorithm.
  #'
  #' @param expression.init a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #' 
  #' @return a named factor that associates each cell to their cluster prediction.
  #' 
}