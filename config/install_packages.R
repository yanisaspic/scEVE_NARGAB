"Run this script to install every R dependency.

	2025/02/05 @yanisaspic"

CRAN <- "http://cran.us.r-project.org"
if (!require("BiocManager")) install.packages("BiocManager", dependencies=TRUE)
if (!require("devtools")) install.packages("devtools", dependencies=TRUE)

#______________ packageVersion("sceve") = 0.0.0.9000
if (!require("sceve")) devtools::install_github("yanisaspic/sceve.package")

#_________________ packageVersion("scEFSC") = 0.1.0
if (!require("cidr")) devtools::install_github("VCCRI/CIDR")
if (!require("G1DBN")) {
	download.file("https://cran.r-project.org/src/contrib/Archive/G1DBN/G1DBN_3.1.1.tar.gz",
				  "./config/.tmp/G1DBN_3.1.1.tar.gz")
	install.packages("./config/.tmp/G1DBN_3.1.1.tar.gz", type="source", repos=NULL)}
if (!require("keras")) install.packages("keras", dependencies=TRUE, repos=CRAN)
if (!require("pcaReduce")) {
	download.file("https://github.com/JustinaZ/pcaReduce/raw/refs/heads/master/pcaReduce_1.0.tar.gz",
				  "./config/.tmp/pcaReduce_1.0.tar.gz")
	install.packages("./config/.tmp/pcaReduce_1.0.tar.gz", type="source", repos=NULL)}
if (!require("Rphenograph")) devtools::install_github("JinmiaoChenLab/Rphenograph")
if (!require("RaceID")) install.packages("RaceID", dependencies=TRUE, repos=CRAN)
if (!require("SC3")) BiocManager::install("SC3")
if (!require("scDHA")) install.packages("scDHA", dependencies=TRUE, repos=CRAN)
if (!require("SHARP")) devtools::install_github("shibiaowan/SHARP")
if (!require("SIMLR")) BiocManager::install("SIMLR")
if (!require("SINCERA")) devtools::install_github("xu-lab/SINCERA")
if (!require("scEFSC")) {
	download.file("https://github.com/Conan-Bian/scEFSC/raw/refs/heads/main/scEFSC_0.1.0.tar.gz",
				  "./config/.tmp/scEFSC_0.1.0.tar.gz")
	install.packages("./config/.tmp/scEFSC_0.1.0.tar.gz", type="source", repos=NULL)}

#_____________ packageVersion("SAFEclustering") = 2.0
if (!require("SAFEclustering")) devtools::install_github("yycunc/SAFEclustering")

# ___________ packageVersion("SAMEclustering") = 1.10
if (!require("SAMEclustering")) devtools::install_github("yycunc/SAMEclustering")

#_______ packageVersion("clusterExperiment") = 2.22.0
if (!require("howmany")) {
	download.file("https://cran.r-project.org/src/contrib/Archive/howmany/howmany_0.3-1.tar.gz",
				  "./config/.tmp/howmany_0.3-1.tar.gz")
	install.packages("./config/.tmp/howmany_0.3-1.tar.gz", type="source", repos=NULL)}
if (!require("clusterExperiment")) BiocManager::install("clusterExperiment")