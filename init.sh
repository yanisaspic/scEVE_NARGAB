#!/bin/bash
#
#	Run this script to intiialize the datasets and the packages used in the scEVE paper.
#
#	2025/02/13 @yanisaspic

# chmod +x ./datasets/download_datasets.sh
# ./datasets/download_datasets.sh
# python3 ./datasets/setup_datasets.py
Rscript ./config/install_packages.R
chmod +x ./config/dependencies/gpmetis
chmod +x ./config/dependencies/shmetis