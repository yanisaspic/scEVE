"Run this script to benchmark scEVE to its component methods on real datasets.

	2024/05/23 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})

source("./scEVE.R")
source("./src/paper/data.R")
source("./src/paper/methods.R")
source("./src/paper/metrics.R")

dataset <- commandArgs(trailingOnly=TRUE)[[1]]
expression.init <- get_real_dataset(dataset)
ground_truth <- get_ground_truth(expression.init)

benchmark <- get_benchmark.dataset(expression.init, ground_truth, dataset,
                                   params=get_default_hyperparameters(), random_state=0)
write.csv(benchmark, glue("./results/{dataset}.csv"), row.names=FALSE)