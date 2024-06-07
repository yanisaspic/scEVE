"Run this script to benchmark scEVE to its component methods on synthetic datasets.

	2024/06/07 @yanisaspic"

source("./scEVE.R")
source("./src/paper/data.R")
source("./src/paper/methods.R")
source("./src/paper/metrics.R")
suppressPackageStartupMessages({library(glue)})

HYPERPARAMETERS <- list(POPULATIONS=c(1,3,5,7,9),
                        BALANCED=c(TRUE, FALSE),
                        NESTED=c(TRUE, FALSE),
                        SEED=1:30)
HYPERPARAMETERS <- expand.grid(HYPERPARAMETERS)

id <- as.numeric(commandArgs(trailingOnly=TRUE)[[1]])
populations <- HYPERPARAMETERS[id, "POPULATIONS"]
balanced <- HYPERPARAMETERS[id, "BALANCED"]
nested <- HYPERPARAMETERS[id, "NESTED"]
seed <- HYPERPARAMETERS[id, "SEED"]

data("PBMC_10X_param_preset")
parameters.init <- PBMC_10X_param_preset[[1]]

dataset <- glue("P{populations}_B{ifelse(balanced, 'T', 'F')}_N{ifelse(nested, 'T', 'F')}_S{seed}")

expression.init <- get_synthetic_dataset(populations, balanced, nested, seed, parameters.init)
ground_truth <- get_ground_truth(expression.init)
benchmark <- get_benchmark.dataset(expression.init, ground_truth, dataset,
                                   params=get_default_hyperparameters(), random_state=0)
write.csv(benchmark, glue("./results/{dataset}.csv"), row.names=FALSE)