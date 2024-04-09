#!/bin/bash
#
#	Run this script to download the files required to set-up the datasets used in the scEFSC paper.
#
#	2024/04/04 @yanisaspic


#######################################################################
# scEFSC datasets are downloaded according to the Hemberg Lab pipeline.
# see:    https://hemberg-lab.github.io/scRNA.seq.datasets/
#         https://github.com/hemberg-lab/scRNA.seq.datasets
#######################################################################

# Li (2017)__________________________
# accession: GSE81861
# cells: 561
# genes: 55,186
# clusters: 9
# sequencing: SMARTer
# doi: 10.1038/ng.3818
LI_DIR="./downloads/Li"
mkdir $LI_DIR
wget -O $LI_DIR/data.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCell%5FLine%5FCOUNT%2Ecsv%2Egz'
gunzip $LI_DIR/data.csv.gz

# Baron (2016)_______________________
# accession: GSE84133
# cells: 8,569
# genes: 20,125
# clusters: 14
# sequencing: inDrop
# doi: 10.1016/j.cels.2016.08.011
BARON_DIR="./downloads/Baron"
mkdir $BARON_DIR
wget -O $BARON_DIR/data.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84133&format=file'
tar -xvf $BARON_DIR/data.tar -C $BARON_DIR
gunzip $BARON_DIR/*.gz
find $BARON_DIR -type f | grep -v "human" | xargs rm

# Tasic (2016)_______________________
# accession: GSE71585
# cells: 1,679
# genes: 24,057
# clusters: 18
# sequencing: SMARTer
# doi: 10.1038/nn.4216
TASIC_DIR="./downloads/Tasic"
mkdir $TASIC_DIR
wget -O $TASIC_DIR/data.zip 'http://casestudies.brain-map.org/celltax/data/data_download.zip'
unzip $TASIC_DIR/data.zip -d $TASIC_DIR

# Camp (2017)________________________
# accession: GSE81252
# cells: 777
# genes: 19,020
# clusters: 7
# sequencing: SMARTer
# doi: 10.1038/nature22796
# /!\ FPKM matrix
CAMP_DIR="./downloads/Camp"
mkdir $CAMP_DIR
wget -O $CAMP_DIR/data_1.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81252&format=file&file=GSE81252%5Fdata%2Ecast%2Elog2%2Elineage%2Ecsv%2Egz'
wget -O $CAMP_DIR/data_2.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81252&format=file&file=GSE81252%5Fdata%2Ecast%2Elog2%2Eliverbud%2Ecsv%2Egz'
gunzip $CAMP_DIR/*.gz

# Lake (2016)________________________
# accession: phs000833.v3.p1
# cells: 777
# genes: 19,020
# clusters: 7
# sequencing: Fluidigm C1
# doi: 10.1126/science.aaf1204
# /!\ TPM matrix
LAKE_DIR="./downloads/Lake"
mkdir $LAKE_DIR
wget -O $LAKE_DIR/data.csv.gz 'http://genome-tech.ucsd.edu/public/Lake_Science_2016/Lake-2016_Gene_TPM.dat.gz'
wget -O $LAKE_DIR/annotations.txt 'http://genome-tech.ucsd.edu/public/Lake_Science_2016/Lake-2016_Gene_TPM_Sample-annotation.txt'
gunzip $LAKE_DIR/data.csv.gz