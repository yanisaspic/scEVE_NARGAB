"Run this script to conduct a scEVE clustering analysis.

	2025/03/07 @yanisaspic"

source("./R/src/get_scrnaseq_data.R")
suppressPackageStartupMessages(library(feve))

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
data <- get_scrnaseq_data(dataset)

params <- feve::get_parameters("scEVE")
params$figures_path <- glue::glue("./results/analyses/{dataset}")
params$sheets_path <- glue::glue("./results/analyses/{dataset}/{dataset}.xlsx")
feve::feve(data$dataset, params, figures=TRUE, sheets=TRUE)