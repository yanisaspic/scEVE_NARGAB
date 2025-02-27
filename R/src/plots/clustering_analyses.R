"Functions called to draw the summary figure of a clustering analysis.

	2025/02/24 @yanisaspic"

source("./R/src/plots/utils.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggtree)
  library(patchwork)
  library(scales)
  library(tidytree)})

get_tree_data <- function(records.meta) {
  #' Get a data structure representing the hierarchy of clusters predicted.
  #' The robustness of each cluster is also included.
  #' The resulting data structure is suitable for `ggtree::ggtree()`.
  #' 
  #' @param records.meta associating predicted populations to generic information, including:
  #' their `size`, their `robustness`, their `parent` and their `clustering_status`.
  #' 
  #' @return a data structure suitable for `ggtree::ggtree()`.
  #'
  relationships <- data.frame(parent=records.meta[-1,]$parent, node=rownames(records.meta[-1,]))
  tmp <- tidytree::as.phylo(relationships)
  tree_data <- dplyr::as_tibble(tmp)
  coordinates <- ggplot2::fortify(tmp)
  
  # align the node clusters to their specific rows______________________________  
  alignment <- rownames(records.meta)[order(records.meta$robustness, -records.meta$size)]
  records.meta <- records.meta[alignment, ]
  tree_data <- tree_data[match(alignment, tree_data$label), ]
  coordinates <- coordinates[match(alignment, coordinates$label), ]
  # ____________________________________________________________________________
  
  tree_data[, "robustness"] <- records.meta$robustness
  for (col in c("isTip", "x", "y", "branch", "angle")) {tree_data[, col] <- coordinates[, col]}
  return(tree_data)
}

get_tree_plot <- function(tree_data) {
  #' Get a plot representing the hierarchy of clusters predicted.
  #' 
  #' @param tree_data a data structure suitable for `ggtree::ggtree()`.
  #'
  #' @return a plot.
  #' 
  is_robust <- function(robustness) {ifelse(robustness > 0, "yes", "no")}
  tree_data[, "is_robust"] <- sapply(X=tree_data$robustness, FUN=is_robust)
  
  plot <- ggtree::ggtree(tree_data) +
    ggplot2::geom_label(ggplot2::aes(label=label, fill=is_robust, color=is_robust),
                        hjust=ifelse(tree_data$isTip, "right", "middle")) +
    ggplot2::scale_color_manual(values=c("white", "black")) +
    ggplot2::scale_fill_manual(values=list("yes"="grey90", "no"="grey10")) +
    ggplot2::geom_label(data=tree_data[-1,], ggplot2::aes(x=branch, label=round(robustness, 2)),
                        label.size=NA, hjust="right") +
    ggplot2::theme(legend.position="none")
  
  return(plot)
}

get_cluster_compositions <- function(records.cells) {
  #' Get a data.frame associating each cluster to its cell composition,
  #' according to the ground truth.
  #' 
  #' @param records.cells a data.frame associating cells to their predicted populations.
  #' Its rows are cells and and its columns are population. The cell values range from 0 to 1.
  #' 
  #' @return a data.frame with three columns: `cluster`, `cell_type` and `n`.
  #' 
  get_ground_truth <- function(cell_id) {
    tmp <- strsplit(cell_id, split="_")[[1]]
    ground_truth <- paste(head(tmp, -1), collapse="_")
    return(ground_truth)}
  
  get_composition.cluster <- function(cluster) {
    cells_of_cluster <- sceve::get_cells_of_population(cluster, records.cells)
    labels <- sapply(X=cells_of_cluster, FUN=get_ground_truth)
    tmp <- table(labels)
    composition.cluster <- data.frame(cluster=cluster, cell_type=names(tmp), n=as.numeric(tmp))
    return(composition.cluster)}
  
  cluster_compositions <- lapply(X=colnames(records.cells), FUN=get_composition.cluster)
  cluster_compositions <- do.call(rbind, cluster_compositions)
  return(cluster_compositions)
}

get_barplot.composition <- function(tree_data, records.cells) {
  #' Get a barplot associating each leaf cluster to its cell type composition.
  #' 
  #' @param tree_data a data structure suitable for `ggtree::ggtree()`.
  #' @param records.cells a data.frame associating cells to their predicted populations.
  #' Its rows are cells and and its columns are population. The cell values range from 0 to 1.
  #' 
  #' @return a plot.
  #' 
  leaves <- tree_data[tree_data$isTip, ]
  alignment <- leaves$label[order(leaves$y)]
  cluster_compositions <- get_cluster_compositions(records.cells)
  cluster_compositions <- cluster_compositions[cluster_compositions$cluster %in% leaves$label, ]
  n_cell_types <- length(unique(cluster_compositions$cell_type))
  
  plot <- ggplot2::ggplot(data=cluster_compositions) +
    ggplot2::geom_bar(ggplot2::aes(x=factor(cluster, levels=alignment), y=n, fill=cell_type),
                      position="fill", stat="identity", color="black") +
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=0), labels=scales::percent) +
    ggplot2::scale_fill_manual(values=pals::kelly(n_cell_types)) +
    ggplot2::coord_flip()
  
  n_rows_legend <- ceiling(n_cell_types / 7)

  plot <- plot +
    ggplot2::theme_classic() +
    ggplot2::ylab("cell types") +
    ggplot2::guides(fill=ggplot2::guide_legend(nrow=n_rows_legend))
  return(plot)
}

get_cluster_sizes <- function(records.cells) {
  #' Get a data.frame associating each cluster to its size.
  #' 
  #' @param records.cells a data.frame associating cells to their predicted populations.
  #' Its rows are cells and and its columns are population. The cell values range from 0 to 1.
  #' 
  #' @return a data.frame with two columns: `cluster` and `n`.
  #' 
  get_size.cluster <- function(cluster) {
    cells_of_cluster <- sceve::get_cells_of_population(cluster, records.cells)
    size.cluster <- data.frame(cluster=cluster, n=length(cells_of_cluster))
    return(size.cluster)}
  
  cluster_sizes <- lapply(X=colnames(records.cells), FUN=get_size.cluster)
  cluster_sizes <- do.call(rbind, cluster_sizes)
  return(cluster_sizes)
}

get_barplot.sizes <- function(tree_data, records.cells) {
  #' Get a barplot associating each leaf cluster to its size.
  #' 
  #' @param tree_data a data structure suitable for `ggtree::ggtree()`.
  #' @param records.cells a data.frame associating cells to their predicted populations.
  #' Its rows are cells and and its columns are population. The cell values range from 0 to 1.
  #' 
  #' @return a plot.
  #' 
  leaves <- tree_data[tree_data$isTip, ]
  alignment <- leaves$label[order(leaves$y)]
  cluster_sizes <- get_cluster_sizes(records.cells)
  cluster_sizes <- cluster_sizes[cluster_sizes$cluster %in% leaves$label, ]
  
  plot <- ggplot2::ggplot(data=cluster_sizes) +
    ggplot2::geom_bar(ggplot2::aes(x=factor(cluster, levels=alignment), y=n),
                      stat="identity", fill="black") +
    ggplot2::geom_text(aes(x=factor(cluster, levels=alignment), y=n, label=n), color="white",
                       hjust=1.5, fontface="bold") +
    ggplot2::scale_y_log10(expand=ggplot2::expansion(mult=0), guide="axis_logticks") +
    ggplot2::coord_flip()

  plot <- plot +
    ggplot2::theme_classic() +
    ggplot2::ylab("# cells")
  return(plot)
}

get_summary_plot <- function(records) {
  #' Get a composite plot summarizing a scEVE clustering analysis on a dataset with a ground truth.
  #' The composite reports the relationships, the cell composition and the size of every leaf cluster
  #' predicted by scEVE.
  #' 
  #' @param records a named list, with four data.frames: `cells`, `markers`, `meta` and `methods`.
  #'
  #' @return a composite plot.
  #' 
  tree_data <- get_tree_data(records$meta)
  tree_plot <- get_tree_plot(tree_data)
  tree_plot <- tree_plot +
    ggplot2::theme(plot.margin=ggplot2::unit(c(0, 0, 0, 0), "null")) +
    ggplot2::guides(fill="none", color="none")
  
  remove_y_axis <- function(barplot) {
    barplot <- barplot +
      ggplot2::theme(axis.line.y=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank())
    return(barplot)}
  
  barplot.composition <- get_barplot.composition(tree_data, records$cells)
  barplot.sizes <- get_barplot.sizes(tree_data, records$cells)
  
  composite_plot <- tree_plot +
    remove_y_axis(barplot.composition) +
    remove_y_axis(barplot.sizes) +
    patchwork::plot_layout(widths=c(12, 3.25, 3.25), guides="collect") &
    ggplot2::theme(legend.position="bottom")
  
  return(composite_plot)
}