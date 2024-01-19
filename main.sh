#!/bin/bash
#
#   Run the benchmark of the ensemblist method.
#
#   2023/12/13 @yanisaspic


### Optional: set-up datasets ###
#################################
for arg in "$@"; 
do
    case $arg in
        -gold)
        python3 ./src/get_data_from_samples.py
            # merge real Smart-Seq2 scRNA-seq measurements and their meta-data (<10k cells, 20k genes)
            ;;
    esac
done


DATASETS=("Darmanis_GBM")
METHODS=(Seurat monocle3 SHARP densityCut)


#################
#   ANALYSIS    #
#################
for data in "${DATASETS[@]%.csv}";  # %.csv: trim the extension
do
    rm ./results/figures/* ./results/tmp/*
    cp ./data/$data.csv ./results/tmp/ALL_DATA.csv
    Rscript ./src/main.R
    cp -R ./results $data
done
