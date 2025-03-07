"Run this script to draw result figures.

	2025/03/07 @yanisaspic"

suppressPackageStartupMessages(library(glue))
source("./R/src/plots/clustering_analyses.R")
source("./R/src/plots/real_benchmarks.R")
source("./R/src/plots/synthetic_benchmarks.R")
source("./R/src/plots/utils.R")

analyses <- get_analyses("./results/analyses")
benchmarks <- get_benchmarks("./results/benchmarks")

# Figure 5 and Supplementary Data : summary of a clustering analysis
for (dataset in names(analyses)) {
    widths <- c(8, 4, 4)
    if (dataset == "Peng_HumPDAC") {widths <- c(11.5, 4, 4)}
    plot <- get_summary_plot(analyses[[dataset]], widths)
    if (dataset == "Peng_HumPDAC") {ggplot2::ggsave(glue::glue("./results/plots/{dataset}.png"), plot, width=24, height=13.5, dpi=300)}
    else {ggplot2::ggsave(glue::glue("./results/plots/{dataset}.png"), plot, width=16, height=9, dpi=300)}}

# Figures 7, 8 and Supplementary Data : benchmark on real datasets
plots.real_benchmark <- get_plots.real_benchmark(benchmarks[benchmarks$is_real, ])
for (plot_name in names(plots.real_benchmark)) {
    ggplot2::ggsave(glue::glue("./results/plots/real_{plot_name}.png"), 
    plots.real_benchmark[[plot_name]], width=16, height=9, dpi=300)}

# Figure 9 and Supplementary Data : benchmark on synthetic datasets
plots.synthetic_benchmark <- get_plots.synthetic_benchmark(benchmarks[!benchmarks$is_real, ])
for (plot_name in names(plots.synthetic_benchmark)) {
    ggplot2::ggsave(glue::glue("./results/plots/synthetic_{plot_name}.png"),
    plots.synthetic_benchmark[[plot_name]], width=16, height=9, dpi=300)}