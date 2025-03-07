"Base clustering methods used in the scEVE framework.

	2025/03/07 @yanisaspic"

suppressPackageStartupMessages({
  library(densitycut)
  library(monocle3)
  library(scater)
  library(Seurat)
  library(SeuratWrappers)
  library(SHARP)})

get_variable_genes <- function(expression, n_genes=5000) {
  #' Get the n most variable genes in a scRNA-seq dataset.
  #' The variable genes are identified by calling the function FindVariableFeatures() of Seurat.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param n_genes the number of highly variable genes to sample.
  #'
  #' @return a vector of genes.
  #'
  SeurObj <- Seurat::CreateSeuratObject(expression)
  SeurObj <- Seurat::FindVariableFeatures(SeurObj, nfeatures = n_genes)
  return(Seurat::VariableFeatures(SeurObj))
}

get_expression.selected <- function(expression) {
  #' Get a scRNA-seq dataset of raw count expression, with selected genes.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @return a scRNA-seq dataset of raw count expression, with selected genes.
  #'
  selected_genes <- get_variable_genes(expression)
  expression.selected <- expression[selected_genes,]
  return(expression.selected)
}

get_SeuratObject.selected <- function(expression.selected) {
  #' Get a SeuratObject from a scRNA-seq dataset of raw count expression, with selected genes.
  #' This function is used as a pre-processing step for the Seurat and monocle3 clustering algorithms.
  #'
  #' @param expression.selected a scRNA-seq dataset of raw count expression, with selected genes.
  #' Its rows are genes and its columns are cells.
  #'
  #' @return a SeuratObject, on which the function ScaleData() of Seurat has been applied already.
  #'
  SeurObj <- Seurat::CreateSeuratObject(expression.selected)
  Seurat::VariableFeatures(SeurObj) <- rownames(expression.selected)
  SeurObj <- Seurat::NormalizeData(SeurObj)
  SeurObj <- Seurat::ScaleData(SeurObj, features=Seurat::VariableFeatures(SeurObj))
  return(SeurObj)
}

use_Seurat.benchmark <- function(expression, random_state) {
  #' Predict clusters with the Seurat method.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  expression.selected <- get_expression.selected(expression)
  SeurObj <- get_SeuratObject.selected(expression.selected)
  SeurObj <- Seurat::RunPCA(SeurObj, features=Seurat::VariableFeatures(SeurObj), seed.use=random_state)
  SeurObj <- Seurat::FindNeighbors(SeurObj, features=Seurat::VariableFeatures(SeurObj))
  SeurObj <- Seurat::FindClusters(SeurObj, random.seed=random_state)
  predictions <- Seurat::Idents(SeurObj)
  return(predictions)
}

use_monocle3.benchmark <- function(expression, random_state) {
  #' Predict clusters with the monocle3 method.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  expression.selected <- get_expression.selected(expression)
  SeurObj <- get_SeuratObject.selected(expression.selected)
  SeurObj <- Seurat::RunUMAP(SeurObj, features=Seurat::VariableFeatures(SeurObj), seed.use=random_state)
  CDSObject <- SeuratWrappers::as.cell_data_set(SeurObj)
  CDSObject <- monocle3::cluster_cells(CDSObject, random_seed=random_state)
  predictions <- CDSObject@clusters@listData$UMAP$clusters
  return(predictions)
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

use_SHARP.benchmark <- function(expression, random_state) {
  #' Predict clusters with the SHARP method.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  expression.selected <- get_expression.selected(expression)
  results <- SHARP::SHARP(scExp=expression.selected, exp.type="count",
                          n.cores = 1, rN.seed=random_state)
  predictions <- get_formatted_predictions(cells=colnames(expression.selected),
                                           predictions=results$pred_clusters)
  return(predictions)
}

use_densityCut.benchmark <- function(expression, random_state) {
  #' Predict clusters with the densityCut method.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  expression.selected <- get_expression.selected(expression)
  logtpm.population <- log2(scater::calculateTPM(expression.selected) + 1)
  set.seed(random_state)
  data <- t(logtpm.population) # densityCut expects cells as rows and genes as columns.
  results <- densitycut::DensityCut(t(logtpm.population), show.plot = FALSE)
  predictions <- get_formatted_predictions(cells=colnames(logtpm.population), predictions=results$cluster)
  return(predictions)
}