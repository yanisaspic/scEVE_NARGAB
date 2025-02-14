#!/bin/bash
#
#	Run this script to download some of the scRNA-seq datasets used in the scEVE paper.
#   cf. https://hemberg-lab.github.io/scRNA.seq.datasets/
#
#	2025/02/13 @yanisaspic

DOWNLOADS_PATH="./datasets/.tmp"

# Li (2017)__________________________
# accession: GSE81861
# cells: 561
# genes: 55,186
# clusters: 9
# sequencing: SMARTer
# doi: 10.1038/ng.3818
LI_PATH="$DOWNLOADS_PATH/Li"
mkdir $LI_PATH
wget -O $LI_PATH/data.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCell%5FLine%5FCOUNT%2Ecsv%2Egz'
gunzip $LI_PATH/data.csv.gz

# Baron (2016)_______________________
# accession: GSE84133
# cells: 8,569
# genes: 20,125
# clusters: 14
# sequencing: inDrop
# doi: 10.1016/j.cels.2016.08.011
BARON_PATH="$DOWNLOADS_PATH/Baron"
mkdir $BARON_PATH
wget -O $BARON_PATH/data.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84133&format=file'
tar -xvf $BARON_PATH/data.tar -C $BARON_PATH
gunzip $BARON_PATH/*.gz
rm $BARON_PATH/data.tar

# Tasic (2016)_______________________
# accession: GSE71585
# cells: 1,679
# genes: 24,057
# clusters: 18
# sequencing: SMARTer
# doi: 10.1038/nn.4216
TASIC_PATH="$DOWNLOADS_PATH/Tasic"
mkdir $TASIC_PATH
wget -O $TASIC_PATH/data.zip 'http://casestudies.brain-map.org/celltax/data/data_download.zip'
unzip $TASIC_PATH/data.zip -d $TASIC_PATH