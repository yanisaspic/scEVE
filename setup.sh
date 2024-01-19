#!/bin/bash
#
#	Run this script to install all the required dependencies.
#
#	2023/05/04 @yanisaspic

pip3 install pipreqs
pipreqs --force --print . > config/requirements.txt
pip3 install -r config/requirements.txt

Rscript config/requirements.R
