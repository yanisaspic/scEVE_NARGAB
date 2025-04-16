"Functions called to draw the results of a benchmark conducted on real datasets.

	2025/04/11 @yanisaspic"

source("./R/src/plots/utils.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(egg)
  library(ggplot2)
  library(ggpubr)
  library(RColorBrewer)})

# _____________________________________________________________________________ real benchmarks

get_fontfaces <- function(benchmark, metric, higher_is_better) {
  #' Get a vector associating the best values on every dataset to a bold fontface.
  #'
  #' @param benchmark a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #' @param metric a character.
  #' @param higher_is_better a boolean indicating if higher values are better for a given metric.
  #'
  #' @return a vector of fontfaces.
  #'
  benchmark <- benchmark[!is.na(benchmark[, metric]), ]
  if (!higher_is_better) {benchmark[, metric] <- -1 * benchmark[, metric]}
  benchmark <- benchmark %>%
    group_by(dataset) %>%
    mutate(fontfaces = ifelse(.data[[metric]] == max(.data[[metric]]), "bold", "plain")) %>%
    ungroup()
  fontfaces <- benchmark$fontfaces
  return(fontfaces)
}

get_gradient <- function(benchmark, metric, higher_is_better) {
  #' Get a gradient associating three colors to the best, the worst and the median values of a heatmap.
  #' 
  #' @param benchmark a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #' @param metric a character.
  #' @param higher_is_better a boolean indicating if higher values are better for a given metric.
  #' 
  #' @return a ScaleContinuous object.
  #' 
  minimum <- min(benchmark[, metric], na.rm=TRUE)
  maximum <- max(benchmark[, metric], na.rm=TRUE)
  midpoint <- mean(benchmark[, metric], na.rm=TRUE)
  colormap <- c("#FB9A99", "#FFFFFF", "#B2DF8A")
  if (!higher_is_better) {colormap <- rev(colormap)}
  
  gradient <- ggplot2::scale_fill_gradient2(low=colormap[1], mid=colormap[2], high=colormap[3],
                                            midpoint=midpoint, limits=c(minimum, maximum),
                                            breaks=c(minimum, midpoint, maximum),
                                            na.value="grey90", labels=function(x) {sprintf("%.2f", x)})
  return(gradient)
}

get_heatmap <- function(benchmark, metric, higher_is_better) {
  #' Get a heatmap associating scRNA-seq datasets (in rows) and clustering methods (in columns).
  #' The value of an input metric is reported in the cells of the heatmap.
  #'
  #' @param benchmark a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #' @param metric a character.
  #' @param higher_is_better a boolean indicating if higher values are better for a given metric.
  #' 
  #' @return a plot.
  #' 
  aes_params <- get_aesthetic_parameters()
  benchmark <- benchmark %>%
    dplyr::mutate(method_type=aes_params$method_types[method])
  fontfaces <- get_fontfaces(benchmark, metric, higher_is_better)
  
  plot <- ggplot2::ggplot(benchmark, ggplot2::aes(x=method, y=dataset, fill=.data[[metric]])) +
    ggplot2::geom_tile(color="white", lwd=.5, linetype=1) +
    ggplot2::geom_text(data=benchmark[!is.na(benchmark[, metric]), ], aes(label=round(.data[[metric]], 2)),
                       fontface=fontfaces, size=4) +
    ggplot2::geom_text(data=benchmark[is.na(benchmark[, metric]), ], aes(label="NA")) +
    ggplot2::facet_grid(~method_type, scales="free_x", space="free_x") +
    get_gradient(benchmark, metric, higher_is_better)
  
  plot$labels$fill <- paste(metric, ifelse(higher_is_better, "[\u2b08]", "[\u2b0a]"), sep=" ")
  
  plot <- plot +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=30, vjust=0.8, hjust=0.8),
                   legend.key.height=ggplot2::unit(0.5, "lines"), legend.title.position="top",
                   legend.position="bottom", axis.title.x=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(), legend.title=ggplot2::element_text(hjust=0.6),
                   strip.text=ggplot2::element_blank())
  return(plot)
}

get_boxplots <- function(benchmark, metric, draw_legend) {
  #' Get boxplots associating clustering methods (on the x-axis) to clustering performances (on the y-axis).
  #'
  #' @param benchmark a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #' @param metric a character.
  #' @param draw_legend a boolean indicating if the plot legend should be drawn.
  #' 
  #' @return a plot.
  #' 
  aes_params <- get_aesthetic_parameters()
  benchmark <- benchmark %>%
    dplyr::mutate(method_type=aes_params$method_types[method])
  
  plot <- ggplot2::ggplot(benchmark[!is.na(benchmark[, metric]), ], ggplot2::aes(x=method, y=.data[[metric]])) +
    ggplot2::geom_boxplot(ggplot2::aes(fill=method)) +
    ggplot2::geom_point() +
    ggplot2::scale_fill_manual(values=aes_params$colormap) +
    ggplot2::facet_grid(~method_type, scales="free_x", space="free_x")
  
  plot <- plot +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(),
                   legend.position="top", panel.grid.major=ggplot2::element_line(linewidth=0.5),
                   legend.text=ggplot2::element_text(size=10), legend.margin=ggplot2::margin(0, 0, 0, 100),
                   axis.line.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(),
                   legend.title=ggplot2::element_blank(), strip.text=ggplot2::element_blank()) +
    ggplot2::guides(fill=ggplot2::guide_legend(nrow=1))
  
  if (!draw_legend) {plot <- plot + ggplot2::theme(legend.position="none")}
  return(plot)
}

get_boxplot_legend <- function(benchmark) {
  #' Get a boxplot legend stored in a plot.
  #' 
  #' @param benchmark a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #' 
  #' @return a plot.
  #' 
  #' 
  benchmark <- benchmark[!benchmark$method %in% c("RSEC*", "scEVE*"), ]
  tmp <- get_boxplots(benchmark, "SI", draw_legend=TRUE)
  legend <- ggpubr::get_legend(tmp)
  legend <- ggpubr::as_ggplot(legend)
  return(legend)
}

get_subplot <- function(benchmark, metric, higher_is_better, draw_legend) {
  #' Get a composite plot associating a heatmap and a boxplot.
  #'
  #' @param benchmark a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #' @param metric a character.
  #' @param higher_is_better a boolean indicating if higher values are better for a given metric.
  #' @param draw_legend a boolean indicating if the plot legend should be drawn.
  #' 
  #' @return a plot.
  #' 
  heatmap <- get_heatmap(benchmark, metric, higher_is_better)
  boxplots <- get_boxplots(benchmark, metric, draw_legend)
  plot <- egg::ggarrange(boxplots, heatmap, nrow=2, ncol=1, widths=1, heights=c(3,7), byrow=TRUE)
  return(plot)
}

get_plot.clustering_performance.real_benchmark <- function(benchmark, metrics) {
  #' Get a composite plot reporting an extrinsic and an intrinsic clustering metric, respectively.
  #' 
  #' @param benchmark a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #' @param metrics a vector of two metrics. 
  #' 
  #' @return a plot.
  #' 
  left_plot <- get_subplot(benchmark[benchmark$method != "ground_truth", ], metrics[1], 
                                higher_is_better=TRUE, draw_legend=FALSE)
  right_plot <- get_subplot(benchmark, metrics[2], higher_is_better=TRUE, draw_legend=FALSE)
  plot <- ggpubr::ggarrange(left_plot, right_plot, nrow=1, ncol=2, labels=c("a", "b"))
  plot <- ggpubr::ggarrange(get_boxplot_legend(benchmark), plot, nrow=2, ncol=1, heights=c(1,10))
  return(plot)
}

get_plot.computational_performance.real_benchmark <- function(benchmark) {
  #' Get a composite plot reporting the time and peak memory usage, respectively.
  #' 
  #' @param benchmark a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #' 
  #' @return a plot.
  #' 
  benchmark <- benchmark[!benchmark$method %in% c("RSEC*", "scEVE*", "ground_truth"), ]
  left_plot <- get_subplot(benchmark, "log10(s)", higher_is_better=FALSE, draw_legend=FALSE)
  right_plot <- get_subplot(benchmark, "log10(Mbytes)", higher_is_better=FALSE, draw_legend=FALSE)
  plot <- ggpubr::ggarrange(left_plot, right_plot, nrow=1, ncol=2, labels=c("a", "b"))
  plot <- ggpubr::ggarrange(get_boxplot_legend(benchmark), plot, nrow=2, ncol=1, heights=c(1,10))
  return(plot)
}

get_plots.real_benchmark <- function(benchmark) {
  #' Get plots summarizing the benchmark conducted on real datasets.
  #'
  #' @param benchmark a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `nPurity`, `SI`, `n_samples`, `dataset` and `is_real`.
  #' 
  #' @return a named list, with three names: `main_clustering_performance`, `supplementary_clustering_performance`
  #' and `computational_performance`.
  #' 
  aes_params <- get_aesthetic_parameters()
  benchmark$dataset <- factor(benchmark$dataset, levels=aes_params$sorted_datasets)
  plots.real_benchmark <- list(
    main_clustering_performance=get_plot.clustering_performance.real_benchmark(benchmark, c("NMI", "SI")),
    supplementary_clustering_performance=get_plot.clustering_performance.real_benchmark(benchmark, c("ARI", "nPurity")),
    computational_performance=get_plot.computational_performance.real_benchmark(benchmark))
  return(plots.real_benchmark)
}