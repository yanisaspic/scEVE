"Run this script to benchmark scEVE to its component methods on synthetic datasets.

	2024/06/07 @yanisaspic"

source("./scEVE.R")
source("./src/paper/data.R")
source("./src/paper/methods.R")
source("./src/paper/metrics.R")
suppressPackageStartupMessages({library(glue)})

HYPERPARAMETERS <- list(POPULATIONS=c(1,3,5,7,9),
                        BALANCED=c(TRUE, FALSE),
                        RELATED=c(TRUE, FALSE),
                        SEED=1:30)
HYPERPARAMETERS <- expand.grid(HYPERPARAMETERS)

id <- as.numeric(commandArgs(trailingOnly=TRUE)[[1]])
populations <- HYPERPARAMETERS[id, "POPULATIONS"]
balanced <- HYPERPARAMETERS[id, "BALANCED"]
related <- HYPERPARAMETERS[id, "RELATED"]
seed <- HYPERPARAMETERS[id, "SEED"]
dataset <- glue("P{populations}_B{ifelse(balanced,'T','F')}_R{ifelse(related,'T','F')}_S{seed}")

data("PBMC_10X_param_preset")
parameters.init <- PBMC_10X_param_preset[[1]]
expression.init <- get_synthetic_dataset(populations, balanced, related, seed, parameters.init)
cell_ids <- colnames(expression.init)
ground_truth <- get_ground_truth(cell_ids)
gc()

params <- get_default_hyperparameters()
params$figures_dir <- glue("./{dataset}")
params$records_file <- glue("./{dataset}.xlsx")
benchmark <- get_benchmark(expression.init, ground_truth, dataset, params,
                           random_state=0, save=TRUE, figures=TRUE)
write.csv(benchmark, glue("./benchmark/{dataset}.csv"), row.names=FALSE)

file.rename(from=params$figures_dir, to=glue("./results/figures/{dataset}"))
file.rename(from=params$records_file, to=glue("./results/records/{dataset}.xlsx"))
