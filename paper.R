"Run this script to generate the comparisons required for the scEVE paper.

	2024/03/18 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})
source("./scEVE.R")
source("./src/paper/methods.R")
source("./src/paper/metrics.R")

#_______________________________________________________________________showcase
expression.init <- read.csv("./data/datasets/Darmanis_HumGBM.csv", header=TRUE, row.names=1)
ground_truth <- get_ground_truth(expression.init)
output <- do_scEVE(expression.init)
figure_1.tree <- draw_tree(output$records, ground_truth)

#______________________________________________________________________benchmark
benchmark.scEFSC <- read.csv("./benchmark_scEFSC.csv", header=TRUE, row.names=1)
datasets <- c("./data/datasets/Baron_HumPan.csv",
              "./data/datasets/Li_HumCRC.csv",
              "./data/datasets/Tasic_MouBra.csv")
get_benchmark.scEVE()
benchmark.scEVE <- lapply(X=datasets, FUN=)