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
cell_ids <- colnames(expression.init)
ground_truth <- get_ground_truth(cell_ids)

params <- get_default_hyperparameters()
params$figures_dir <- glue("./{dataset}")
params$records_file <- glue("./{dataset}.xlsx")
benchmark <- get_benchmark(expression.init, ground_truth, dataset, params,
                           random_state=0, save=TRUE, figures=TRUE)
write.csv(benchmark, glue("./benchmark/{dataset}.csv"), row.names=FALSE)