#!/bin/bash
#
#	Run this script to set-up the datasets used in the paper.
#
#	2024/09/20 @yanisaspic

SETUP_DIR="./etc/source"
chmod +x $SETUP_DIR/download.sh
./$SETUP_DIR/download.sh
python3 $SETUP_DIR/setup.py
