"Run this script to conduct a scEVE clustering analysis.

	2025/04/11 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
suppressPackageStartupMessages(library(feve))

dataset_label <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset_label)

params <- feve::get_parameters("scEVE")
params$figures_path <- glue::glue("./results/analyses/{dataset_label}")
params$sheets_path <- glue::glue("./results/analyses/{dataset_label}/{dataset_label}.xlsx")
feve::feve(data$dataset, params, figures=TRUE, sheets=TRUE)