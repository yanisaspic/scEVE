#!/bin/bash
#
#	Run this script to set-up the datasets used in the paper.
#
#	2024/05/23 @yanisaspic

module load r/4.3.1
module load python/3.9

CONFIG_DIR="./etc/config"
Rscript $CONFIG_DIR/install_requirements.R

SETUP_DIR="./etc/setup_data"
chmod +x $SETUP_DIR/real/download_source.sh
./$SETUP_DIR/real/download_source.sh
python3 $SETUP_DIR/real/setup_data.py
Rscript $SETUP_DIR/synthetic/generate_data.R
