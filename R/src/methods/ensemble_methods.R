"Ensemble clustering methods used in the scEVE paper.

	2025/02/14 @yanisaspic"

suppressPackageStartupMessages({
    # library(clusterExperiment)
    # library(SAFEclustering)
    library(SAMEclustering)})
    # library(scEFSC)})

use_RSEC.benchmark <- function(expression, random_state) {
  #' Predict clusters with the RSEC method.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  base_clusters <- SAMEclustering::individual_clustering(inputTags = expression, mt_filter = TRUE,
  percent_dropout = 10, SC3 = TRUE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nGene_filter = FALSE,
  nPC.seurat = NULL, resolution = 0.7, tSNE = TRUE, dimensions = 2, perplexity = 30, SIMLR = TRUE, diverse = TRUE, 
  save.results = FALSE, SEED = random_state)
  clusters <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED = random_state)
  return(clusters)
}

use_SAFE.benchmark <- function(expression, random_state) {
  #' Predict clusters with the SAFE method.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  base_clusters <- SAMEclustering::individual_clustering(inputTags = expression, mt_filter = TRUE,
  percent_dropout = 10, SC3 = TRUE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nGene_filter = FALSE,
  nPC.seurat = NULL, resolution = 0.7, tSNE = TRUE, dimensions = 2, perplexity = 30, SIMLR = TRUE, diverse = TRUE, 
  save.results = FALSE, SEED = random_state)
  clusters <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED = random_state)
  return(clusters)
}

use_SAME.benchmark <- function(expression, random_state) {
  #' Predict clusters with the SAME method.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  base_clusters <- SAMEclustering::individual_clustering(inputTags = expression, mt_filter = TRUE,
  percent_dropout = 10, SC3 = TRUE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nGene_filter = FALSE,
  nPC.seurat = NULL, resolution = 0.7, tSNE = TRUE, dimensions = 2, perplexity = 30, SIMLR = TRUE, diverse = TRUE, 
  save.results = FALSE, SEED = random_state)
  clusters <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED = random_state)
  return(clusters)
}

use_scEFSC.benchmark <- function(expression, random_state) {
  #' Predict clusters with the scEFSC method.
  #'
  #' @param expression a scRNA-seq dataset of raw count expression, without selected genes.
  #' Its rows are genes and its columns are cells.
  #' @param random_state a numeric.
  #'
  #' @return a named factor that associates each cell to their cluster prediction.
  #'
  base_clusters <- SAMEclustering::individual_clustering(inputTags = expression, mt_filter = TRUE,
  percent_dropout = 10, SC3 = TRUE, CIDR = TRUE, nPC.cidr = NULL, Seurat = TRUE, nGene_filter = FALSE,
  nPC.seurat = NULL, resolution = 0.7, tSNE = TRUE, dimensions = 2, perplexity = 30, SIMLR = TRUE, diverse = TRUE, 
  save.results = FALSE, SEED = random_state)
  clusters <- SAMEclustering(Y = t(cluster.result), rep = 3, SEED = random_state)
  return(clusters)
}