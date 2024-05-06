#!/bin/bash
#
#	Run this script to set-up the datasets used in the paper.
#
#	2024/04/24 @yanisaspic

CONFIG_DIR="./etc/config"
Rscript $CONFIG_DIR/install_requirements.R

SETUP_DIR="./etc/setup_data"
chmod +x $SETUP_DIR/download_data.sh
./$SETUP_DIR/download_data.sh
python3 $SETUP_DIR/setup_data.py
