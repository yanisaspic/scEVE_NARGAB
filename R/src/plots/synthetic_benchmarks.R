"Functions called to draw the results of a benchmark conducted on synthetic datasets.

	2025/02/21 @yanisaspic"

source("./R/src/plots/utils.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)})

# ________________________________________________________________________ synthetic benchmarks

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

get_lineplots <- function(benchmark, metric, draw_legend) {
  #' Get lineplots associating clustering methods to their performance (y-axis) on multiple datasets (x-axis)
  #' 
  #' @param benchmark a data.frame with fourteen columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `Purity`, `SI`, `n_cells`, `dataset`, `is_real`, `k`, `balanced_clusters`,
  #' `independent_clusters` and `random_state`.
  #' @param metrics a character
  #' @param draw_legend a boolean indicating if the plot legend should be drawn.
  #' 
  #' @return a plot.
  #' 
  aes_params <- get_aesthetic_parameters()
  tmp <- benchmark %>%
    dplyr::group_by(k, balanced_clusters, independent_clusters, method) %>%
    dplyr::summarise(y = mean(.data[[metric]], na.rm=TRUE),
                     se = stats::sd(.data[[metric]] / sqrt(dplyr::n()), na.rm=TRUE),
                     .groups = "drop")
  
  # this section is dedicated to aesthetics ____________________________________
  linetypes <- ifelse(tmp$method=="scEVE*", "dashed", "solid")
  linetypemaps <- stats::setNames(linetypes, tmp$method)
  labelmaps.balanced <- stats::setNames(c("balanced cluster sizes", "unbalanced cluster sizes"),
                                        c(TRUE, FALSE))
  labelmaps.independent <- stats::setNames(c("unrelated cluster transcriptomes",
                                             "related cluster transcriptomes"), c(TRUE, FALSE))
  tmp$balanced_clusters <- factor(tmp$balanced_clusters, levels=c(TRUE, FALSE))
  tmp$independent_clusters <- factor(tmp$independent_clusters, levels=c(TRUE, FALSE))
  # ____________________________________________________________________________

  plot <- ggplot2::ggplot(tmp, ggplot2::aes(x=k, y=y)) +
    ggplot2::geom_line(ggplot2::aes(color=method, linetype=method)) +
    ggplot2::geom_point(ggplot2::aes(color=method)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=y-se, ymax=y+se, color=method), width=0.2) +
    ggplot2::facet_grid(balanced_clusters~independent_clusters, switch="y",
                        labeller = ggplot2::labeller(balanced_clusters=labelmaps.balanced,
                                                     independent_clusters=labelmaps.independent)) +
    ggplot2::scale_linetype_manual(values=linetypemaps) +
    ggplot2::scale_x_continuous(breaks=1:10) +
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=0.05), position="right") +
    ggplot2::scale_color_manual(values=unlist(aes_params$colormap))
     
  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background=ggplot2::element_rect(fill="#333333"),
                   strip.text=ggplot2::element_text(color="white"),
                   panel.background=ggplot2::element_rect(fill="#ebebeb"),
                   panel.grid.major=ggplot2::element_line(colour="white"),
                   legend.text=ggplot2::element_text(size=10), legend.title=ggplot2::element_blank()) +
    ggplot2::ylab(metric) +
    ggplot2::xlab("# clusters") +
    ggplot2::guides(color=ggplot2::guide_legend(nrow=1), linetype="none")
  
  if (!draw_legend) {plot <- plot + ggplot2::theme(legend.position="none")}
  return(plot)
}

get_lineplot_legend <- function(benchmark) {
  #' Get a lineplot legend stored in a plot.
  #' 
  #' @param benchmark a data.frame with fourteen columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `Purity`, `SI`, `n_cells`, `dataset`, `is_real`, `k`, `balanced_clusters`,
  #' `independent_clusters` and `random_state`.
  #' 
  #' @return a plot.
  #' 
  #' 
  benchmark <- benchmark[!benchmark$method=="scEVE*", ]
  tmp <- get_lineplots(benchmark, "SI", draw_legend=TRUE)
  legend <- ggpubr::get_legend(tmp)
  legend <- ggpubr::as_ggplot(legend)
  return(legend)
}

get_lineplot.synthetic_benchmark <- function(benchmark, metrics) {
  #' Get a composite plot reporting an extrinsic and an intrinsic clustering metric, respectively.
  #' 
  #' @param benchmark a data.frame with fourteen columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `Purity`, `SI`, `n_cells`, `dataset`, `is_real`, `k`, `balanced_clusters`,
  #' `independent_clusters` and `random_state`.
  #' @param metrics a vector of two metrics. 
  #' 
  #' @return a plot.
  #' 
  left_plot <- get_lineplots(benchmark, metrics[1], draw_legend=FALSE)
  right_plot <- get_lineplots(benchmark, metrics[2], draw_legend=FALSE)
  plot <- ggpubr::ggarrange(left_plot, right_plot, nrow=1, ncol=2, labels=c("a", "b"))
  plot <- ggpubr::ggarrange(plot, get_lineplot_legend(benchmark), nrow=2, ncol=1, heights=c(10,1))
  return(plot)
}

get_plots.synthetic_benchmark <- function(benchmark) {
  #' Get plots summarizing the benchmark conducted on synthetic datasets.
  #'
  #' @param benchmark a data.frame with ten columns: `method`, `time (s)`, `peak_memory_usage (Mb)`,
  #' `ARI`, `NMI`, `Purity`, `SI`, `n_cells`, `dataset` and `is_real`.
  #' 
  #' @return a named list, with three names: `main_clustering_performance`, `supplementary_clustering_performance`
  #' and `computational_performance`.
  #' 
  synthetic_datasets_info <- get_synthetic_datasets_info()
  synthetic_datasets_info <- synthetic_datasets_info %>% dplyr::mutate(dataset=rownames(.))
  benchmark <- benchmark %>% dplyr::left_join(synthetic_datasets_info, by=c("dataset"="dataset"))
  benchmark <- benchmark[benchmark$method!="ground_truth", ]
  
  plots.synthetic_benchmark <- list(
    main_clustering_performance=get_lineplot.synthetic_benchmark(benchmark, c("NMI", "SI")),
    supplementary_clustering_performance=get_lineplot.synthetic_benchmark(benchmark, c("ARI", "Purity")),
    computational_performance=get_lineplot.synthetic_benchmark(benchmark, c("log10(s)", "log10(Mbytes)")))
  return(plots.synthetic_benchmark)
}