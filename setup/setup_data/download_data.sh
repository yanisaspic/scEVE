#!/bin/bash
#
#	Run this script to download the files required to set-up the datasets.
#
#	2024/04/03 @yanisaspic

### Darmanis: GSE84465
######################


### scEFSC datasets are downloaded according to the Hemberg Lab pipeline.
### see:    https://hemberg-lab.github.io/scRNA.seq.datasets/
###         https://github.com/hemberg-lab/scRNA.seq.datasets
#######################################################################

# Baron______________________________
# accession: GSE84133
# cells: 8,569
# genes: 20,125
# clusters: 14
# sequencing: inDrop
# doi: 10.1016/j.cels.2016.08.011"""
#____________________________________
BARON_DIR="./downloads/Baron"
mkdir $BARON_DIR
wget -O $BARON_DIR/data.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84133&format=file'
tar -xvf $BARON_DIR/data.tar -C $BARON_DIR
gunzip $BARON_DIR/*.gz
find $BARON_DIR -type f | grep -v "human" | xargs rm