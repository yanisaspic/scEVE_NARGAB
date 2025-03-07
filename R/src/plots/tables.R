"Functions called to get tables of results.

	2025/03/07 @yanisaspic"

suppressPackageStartupMessages({
  library(cancersea)
  library(dplyr)
  library(glue)
  library(tibble)
  library(tidyr)})

get_cancer_markers <- function() {
  #' Get a data.frame associating cancer signatures to their marker genes.
  #' 
  #' @return a data.frame with two columns: `cancer_signature` and `symbol`.
  #' 
  f <- function(data, signature) {
    n_markers <- nrow(data)
    data[, "cancer_signature"] <- glue::glue("{signature} ({n_markers})")
    return(data)}
  markers <- list(
    f(cancersea::angiogenesis, "Angiogenesis"),
    f(cancersea::apoptosis, "Apoptosis"),
    f(cancersea::cell_cycle, "Cell cycle"),
    f(cancersea::differentiation, "Differentiation"),
    f(cancersea::dna_damage, "DNA damage"),
    f(cancersea::dna_repair, "DNA repair"),
    f(cancersea::emt, "EMT"),
    f(cancersea::hypoxia, "Hypoxia"),
    f(cancersea::inflammation, "Inflammation"),
    f(cancersea::invasion, "Invasion"),
    f(cancersea::metastasis, "Metastasis"),
    f(cancersea::proliferation, "Proliferation"),
    f(cancersea::quiescence, "Quiescence"),
    f(cancersea::stemness, "Stemness"))
  markers <- do.call(rbind, markers)
  markers <- markers[, -1]
  
  tmp <- unique(markers$symbol)
  label <- glue::glue("Any ({length(tmp)})")
  any <- data.frame(symbol=tmp, cancer_signature=as.character(label))
  markers <- rbind(markers, any)
  return(markers)
}

get_cancer_data <- function(records) {
  #' Get a data.frame associating predicted populations to their cancer marker genes.
  #' 
  #' @param records data.frame associating predicted populations to their marker genes.
  #' Its rows are genes, its columns are predicted populations,
  #' and the strength of the characterization (e.g. log2-transformed pvalues) are reported in the table.
  #'
  #' @return a data.frame with three columns: `cluster`, `cancer_signature` and `n_markers`.
  #'
  clusters <- rownames(records$meta)
  is_leaf <- function(cluster) {!(cluster %in% records$meta$parent)}
  clusters <- clusters[sapply(X=clusters, FUN=is_leaf)]
  get_cluster_markers <- function(cluster) {rownames(records$features[records$features[, cluster] > 0,])}
  cluster_markers <- sapply(X=clusters, FUN=get_cluster_markers)
  cluster_markers <- tibble::enframe(cluster_markers, name="cluster", value="symbol") %>%
    tidyr::unnest(symbol)

  cancer_markers <- get_cancer_markers()
  cancer_data <- cluster_markers %>%
    dplyr::inner_join(cancer_markers, by="symbol") %>%
    dplyr::count(cluster, cancer_signature, name="n_markers") %>%
    tidyr::complete(cluster, cancer_signature, fill=list(n_markers=0))
  average <- cancer_data %>%
    dplyr::group_by(cancer_signature) %>%
    dplyr::summarise(n_markers=mean(n_markers), .groups="drop") %>%
    dplyr::mutate(cluster="mu")
  
  cancer_data <- dplyr::bind_rows(cancer_data, average) %>%
    tidyr::pivot_wider(names_from=cluster, values_from=n_markers) %>%
    tibble::column_to_rownames("cancer_signature") %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits=2)))
  return(cancer_data)
}

get_contributions_data <- function(analyses) {
  #' Get a data.frame associating clustering methods to their contributions in
  #' predicting robust clusters.
  #' 
  #' @param analyses a named list. Each name corresponds to a dataset analysed.
  #' 
  #' @return a data.frame with three columns: `method`, `dataset` and `contribution`.
  #' 
  get_contributions <- function(dataset) {
    records <- analyses[[dataset]]
    contributions <- rowMeans(records$methods)
    contributions <- tibble::enframe(contributions, name="method", value="contribution") %>%
      tidyr::unnest(contribution)
    contributions[, "dataset"] <- dataset
    return(contributions)}
  
  contributions_data <- lapply(X=names(analyses), FUN=get_contributions)
  contributions_data <- do.call(rbind, contributions_data) %>%
    tidyr::complete(method, dataset, fill=list(contribution=0)) %>%
    tidyr::pivot_wider(names_from=method, values_from=contribution) %>%
    tibble::column_to_rownames("dataset")
  
  averages <- contributions_data %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), \(x) mean(x, na.rm=TRUE)))
  ses <- contributions_data %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ stats::sd(.) / sqrt(length(.))))
  contributions_statistics <- rbind(averages, ses)
  rownames(contributions_statistics) <- c("average", "standard_error")
  
  contributions_data <- rbind(contributions_data, contributions_statistics) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) 100*x)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits=2)))
  return(contributions_data)
}

get_leftout_data <- function(benchmarks) {
  #' Get a data.frame associating RSEC and scEVE to their missing cells.
  #' 
  #' @param benchmarks a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #'
  #' @return a data.frame with three columns: `method`, `dataset` and `missing_cells`.
  #' 
  benchmarks <- benchmarks[benchmarks$method %in% c("ground_truth", "RSEC*", "scEVE*"), ]
  leftout_data <- benchmarks %>%
    dplyr::group_by(dataset) %>%
    dplyr::mutate(ref=n_samples[method=="ground_truth"]) %>%
    dplyr::mutate(percent=100*n_samples/ref)
  leftout_data <- leftout_data[leftout_data$method %in% c("RSEC*", "scEVE*"),
                               c("method", "dataset", "percent")]
  
  leftout_data <- leftout_data %>%
    tidyr::pivot_wider(names_from=method, values_from=percent) %>%
    tibble::column_to_rownames("dataset") %>%
    dplyr::mutate_all(~ifelse(. == 0, NA, .))
  
  averages <- leftout_data %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), \(x) mean(x, na.rm=TRUE)))
  ses <- leftout_data %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ stats::sd(., na.rm=TRUE) / sqrt(length(.))))
  leftout_statistics <- rbind(averages, ses)
  rownames(leftout_statistics) <- c("average", "standard_error")

  leftout_data <- rbind(leftout_data, leftout_statistics) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), \(x) round(x, digits=2)))
  return(leftout_data)
}