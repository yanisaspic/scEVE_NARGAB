#!/bin/bash
#
#	Run this script to intiialize the datasets and the packages used in the scEVE paper.
#
#	2025/03/07 @yanisaspic

#_______________________________scRNA-seq datasets
chmod +x ./datasets/download_datasets.sh
./datasets/download_datasets.sh
python3 ./datasets/setup_datasets.py

#_______________________________ensemble algorithms
Rscript ./config/install_packages.R
chmod +x ./config/dependencies/gpmetis
chmod +x ./config/dependencies/shmetis

#_______________________________SAFE virtual environment
singularity build ./config/dependencies/SAFE.sif ./config/dependencies/SAFE.def