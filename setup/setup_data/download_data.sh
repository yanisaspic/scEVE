#!/bin/bash
#
#	Run this script to download the files required to set-up the datasets.
#
#	2024/04/04 @yanisaspic

# Darmanis (2017)____________________
# accession: GSE84465
# cells: 3,589
# genes: 23,460
# clusters: 7
# sequencing: Smart-Seq2
# doi: 10.1016/j.celrep.2017.10.030
#____________________________________
DARMANIS_DIR="./downloads/Darmanis"
mkdir $DARMANIS_DIR
wget -O $DARMANIS_DIR/metadata.tgz 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84465/miniml/GSE84465_family.xml.tgz'
wget -O $DARMANIS_DIR/data.csv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84465&format=file&file=GSE84465%5FGBM%5FAll%5Fdata%2Ecsv%2Egz'
tar -zxvf $DARMANIS_DIR/metadata.tgz -C $DARMANIS_DIR
gunzip $DARMANIS_DIR/data.csv.gz

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
#____________________________________
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
#____________________________________
BARON_DIR="./downloads/Baron"
mkdir $BARON_DIR
wget -O $BARON_DIR/data.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84133&format=file'
tar -xvf $BARON_DIR/data.tar -C $BARON_DIR
gunzip $BARON_DIR/*.gz
find $BARON_DIR -type f | grep -v "human" | xargs rm