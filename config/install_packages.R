"Run this script to install every R packages.

	2025/03/07 @yanisaspic"

CRAN <- "http://cran.us.r-project.org"
if (!require("BiocManager")) install.packages("BiocManager", dependencies=TRUE, repos=CRAN)
if (!require("devtools")) install.packages("devtools", dependencies=TRUE, repos=CRAN)

#______________ packageVersion("feve") = 1.0
if (!require("presto")) devtools::install_github('immunogenomics/presto')
if (!require("feve")) devtools::install_github("yanisaspic/feve@7b80758")

#_____________ packageVersion("SAFEclustering") = 2.0
if (!require("SAFEclustering")) {
	devtools::install_github("yycunc/SAFEclustering")
	download.file("https://github.com/yycunc/SAFEclustering/raw/refs/heads/master/gpmetis_and_shmetis_for_Linux/gpmetis",
				  "./config/dependencies/gpmetis")
	download.file("https://github.com/yycunc/SAFEclustering/raw/refs/heads/master/gpmetis_and_shmetis_for_Linux/shmetis",
				  "./config/dependencies/shmetis")}

# ___________ packageVersion("SAMEclustering") = 1.10
if (!require("SAMEclustering")) devtools::install_github("yycunc/SAMEclustering")

#_______ packageVersion("clusterExperiment") = 2.22.0
if (!require("howmany")) {
	download.file("https://cran.r-project.org/src/contrib/Archive/howmany/howmany_0.3-1.tar.gz",
				  "./config/dependencies/howmany_0.3-1.tar.gz")
	install.packages("./config/dependencies/howmany_0.3-1.tar.gz", type="source", repos=NULL)}
if (!require("clusterExperiment")) BiocManager::install("clusterExperiment")

# _______________________________________________Misc
if (!require("cancersea")) devtools::install_github("camlab-bioml/cancersea")
if (!require("pals")) install.packages("pals", dependencies=TRUE, repos=CRAN)
if (!require("SPARSim")) devtools::install_gitlab("sysbiobig/sparsim")
