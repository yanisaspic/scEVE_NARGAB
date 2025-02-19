"Run this script to conduct a scEVE clustering analysis.

	2025/02/14 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
suppressPackageStartupMessages(library("sceve"))

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset)

params <- sceve::get_default_parameters()
params$figures_path <- glue::glue("./results/analyses/{dataset}")
params$sheets_path <- glue::glue("./results/analyses/{dataset}/{dataset}.xlsx")
sceve(data$expression.init, params, figures=TRUE, sheets=TRUE)