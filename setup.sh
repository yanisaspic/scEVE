#!/bin/bash
#
#	Run this script to set-up the datasets used in the paper.
#
#	2024/04/04 @yanisaspic

cd ./setup/setup_data
chmod +x ./download_data.sh
./download_data.sh
python3 setup_data.py