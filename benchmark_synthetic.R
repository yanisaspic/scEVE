"Run this script to benchmark scEVE to its component methods on synthetic datasets.

	2024/05/23 @yanisaspic"

source("./scEVE.R")
source("./src/paper/data.R")
source("./src/paper/methods.R")
source("./src/paper/metrics.R")
suppressPackageStartupMessages({library(glue)})

id <- as.numeric(commandArgs(trailingOnly=TRUE)[[1]])
HYPERPARAMETERS <- list(N_POPULATIONS=c(1,3,5,7,9),
                        DISTRIBUTION=c("uniform", "geometric"),
                        RANDOM_STATE=1:30)
HYPERPARAMETERS <- expand.grid(HYPERPARAMETERS)
dist <- HYPERPARAMETERS[id, "DISTRIBUTION"]
n_pop <- HYPERPARAMETERS[id, "N_POPULATIONS"]
seed <- HYPERPARAMETERS[id, "RANDOM_STATE"]
dataset <- glue("D{dist}_N{n_pop}_S{seed}")

data("PBMC_10X_param_preset")
parameters.init <- PBMC_10X_param_preset[[1]]
expression.init <- get_synthetic_dataset(n_pop, dist, seed, parameters.init)
ground_truth <- get_ground_truth(expression.init)

benchmark <- get_benchmark.dataset(expression.init, ground_truth, dataset,
                                   params=get_default_hyperparameters(), random_state=0)
write.csv(benchmark, glue("./results/{dataset}.csv"), row.names=FALSE)
