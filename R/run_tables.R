"Functions called to get tables of results.

	2025/03/07 @yanisaspic"

source("./R/src/plots/utils.R")
source("./R/src/plots/tables.R")
analyses <- get_analyses("./results/analyses")
benchmarks <- get_benchmarks("./results/benchmarks")

# Table 4 : characterization with CancerSea
cancer_data <- get_cancer_data(analyses$Darmanis_HumGBM)
cancer_data <- cancer_data[, c("C.L.1", "C.L.L.L", "mu")]
write.table(cancer_data, "./results/plots/CancerSEA.csv", row.names=TRUE, sep=" & ")

# Table 5 : sizes of the predictions
leftout_data <- get_leftout_data(benchmarks[benchmarks$is_real,])
write.table(leftout_data, "./results/plots/leftout.csv", row.names=TRUE, sep=" & ")

# Table 6 : contributions of the clustering methods
contributions_data <- get_contributions_data(analyses)
write.table(contributions_data, "./results/plots/contributions.csv", row.names=TRUE, sep=" & ")