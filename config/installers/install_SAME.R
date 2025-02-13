"Run this script to install the R dependencies required for the SAME clustering algorithm.

To prevent some issues with their usage, dependencies updated since the last update of SAME
are versionned according to the date of the last update (2021-04-04).

	2025/02/05 @yanisaspic"

source("./config/R_LIBS.R")

# if (!require("BiocManager")) install.packages("BiocManager", dependencies=TRUE)
if (!require("devtools")) install.packages("devtools", dependencies=TRUE)

print(.libPaths())

DOWNLOADS_PATH <- "./config/.tmp/"
REPOSITORY_URLS <- list("CRAN"="https://cran.r-project.org/src/contrib/Archive",
						"BioC"="https://mghp.osn.xsede.org/bir190004-bucket01/archive.bioconductor.org/packages/3.12/bioc/src/contrib")

f_inst <- function(url, pkg_archive) {
	path <- paste0(DOWNLOADS_PATH, pkg_archive)
	download.file(url, path)
	install.packages(path, type="source", repos=NULL, dependencies=TRUE)}

f_cran <- function(pkg_archive) {
	pkg <- unlist(strsplit(pkg_archive, split="_"))[1]
	url <- paste0(REPOSITORY_URLS[["CRAN"]], "/", pkg, "/", pkg_archive)
	f_inst(url, pkg_archive)}

f_bioc <- function(pkg_archive) {
	url <- paste0(REPOSITORY_URLS[["BioC"]], "/", pkg_archive)
	f_inst(url, pkg_archive)}

.libPaths(R_LIBS$SAME) # ___________ packageVersion("SAMEclustering") = 1.10
if (!require("Rcpp")) f_cran("Rcpp_1.0.6.tar.gz")
if (!require("Rtsne")) f_cran("Rtsne_0.15.tar.gz")
if (!require("rrcov")) f_cran("rrcov_1.5-5.tar.gz")
if (!require("SC3")) f_bioc("SC3_1.18.0.tar.gz")
if (!require("Seurat")) f_cran("Seurat_4.0.1.tar.gz")
if (!require("e1071")) f_cran("e1071_1.7-6.tar.gz")
if (!require("SingleCellExperiment")) f_bioc("SingleCellExperiment_1.12.0.tar.gz")
if (!require("SummarizedExperiment")) f_bioc("SummarizedExperiment_1.20.0.tar.gz")
if (!require("S4Vectors")) f_bioc("S4Vectors_0.28.1.tar.gz")
if (!require("doRNG")) f_cran("doRNG_1.8.2.tar.gz")
if (!require("cidr")) devtools::install_github("VCCRI/CIDR")
if (!require("SIMLR")) f_bioc("SIMLR_1.16.0.tar.gz")
if (!require("Matrix")) f_cran("Matrix_1.3-2.tar.gz")
if (!require("SAMEclustering")) devtools::install_github("yycunc/SAMEclustering")