"Run this script to generate the benchmark of scEVE.

	2023/01/20 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(SummarizedBenchmark)
})

source("./scEVE.R")
source("./src/benchmark/individual_methods.R")
source("./src/benchmark/ensemble_methods.R")
source("./src/benchmark/metrics.R")

random_state <- 0
expression.init <- read.csv("./data/Darmanis_GBM.csv", header=TRUE, row.names=1)
records <- do_scEVE(expression.init, figures=TRUE)
# 