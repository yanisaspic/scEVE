"Run this script to compute the average similarity between the clustering methods
  integrated in scEVE on the related balanced synthetic datasets.

	2024/09/19 @yanisaspic"

source("./scEVE.R")
source("./src/paper/data.R")
source("./src/paper/methods.R")
source("./src/paper/metrics.R")
suppressPackageStartupMessages({library(glue)})

HYPERPARAMETERS <- list(POPULATIONS=c(1,3,5,7,9),
                        BALANCED=c(TRUE),
                        RELATED=c(TRUE),
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
results.individual <- get_results.individual(expression.init, params, random_state=0,
                                             save=FALSE, figures=FALSE)
similarity.individual <- get_similarity.individual(results.individual, dataset)
write.csv(similarity.individual, glue("./results/similarity/{dataset}.csv"), row.names=FALSE)