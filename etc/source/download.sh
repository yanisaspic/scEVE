#!/bin/bash
#
#	Run this script to download the files required to set-up the datasets used in the scEFSC paper,
#	and to get the cancer signatures from CancerSea.
#
#	2024/05/23 @yanisaspic

DOWNLOADS_DIR="./etc/source/tmp"

# Cancer signatures (CancerSEA)_______
#
CANCER_SIGNATURES=("Angiogenesis" "Apoptosis" "Cell_Cycle" \
"Differentiation" "DNA_damage" "DNA_repair" \
"EMT" "Hypoxia" "Inflammation" "Invasion" \
"Metastasis" "Proliferation" "Quiescence" "Stemness")
CANCER_DIR="$DOWNLOADS_DIR/CancerSEA"
mkdir $CANCER_DIR
for signature in "${CANCER_SIGNATURES[@]}"
do
	wget -O $CANCER_DIR/$signature.txt "http://biocc.hrbmu.edu.cn/CancerSEA/download/signature/$signature.txt"
done

#######################################################################
# scEFSC datasets are downloaded according to the Hemberg Lab pipeline.
# see:    https://hemberg-lab.github.io/scRNA.seq.datasets/
#######################################################################

# Li (2017)__________________________
# accession: GSE81861
# cells: 561
# genes: 55,186
# clusters: 9
# sequencing: SMARTer
# doi: 10.1038/ng.3818
LI_DIR="$DOWNLOADS_DIR/Li"
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
BARON_DIR="$DOWNLOADS_DIR/Baron"
mkdir $BARON_DIR
wget -O $BARON_DIR/data.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84133&format=file'
tar -xvf $BARON_DIR/data.tar -C $BARON_DIR
gunzip $BARON_DIR/*.gz
rm $BARON_DIR/data.tar

# Tasic (2016)_______________________
# accession: GSE71585
# cells: 1,679
# genes: 24,057
# clusters: 18
# sequencing: SMARTer
# doi: 10.1038/nn.4216
TASIC_DIR="$DOWNLOADS_DIR/Tasic"
mkdir $TASIC_DIR
wget -O $TASIC_DIR/data.zip 'http://casestudies.brain-map.org/celltax/data/data_download.zip'
unzip $TASIC_DIR/data.zip -d $TASIC_DIR
