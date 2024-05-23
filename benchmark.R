"Run this script to generate the benchmark comparing scEVE to its component methods.

	2024/05/23 @yanisaspic"

source("./scEVE.R")
source("./src/paper/data.R")
source("./src/paper/methods.R")
source("./src/paper/metrics.R")
suppressPackageStartupMessages({library(glue)})

dataset <- commandArgs(trailingOnly = TRUE)[[1]]

compute_benchmark.dataset <- function(dataset, params, random_state) {
  #' Get the results of scEVE and the individual clustering methods it uses on a scRNA-seq dataset.
  #' These results include:
  #' - peakRAM (Mbytes): the maximum memory usage of the method.
  #' - time (s): the computation time in seconds.
  #' - ARI: a clustering performance metric.
  #' - NMI: a clustering performance metric.
  #'
  #' @param dataset: a valid scRNA-seq dataset label.
  #' @param params: a list of parameters.
  #' @param random_state: a numeric.
  #'
  #' @return a data.frame with 6 columns: 'peakRAM', 'time', 'ARI', 'NMI', 'method' and 'dataset'.
  #'
  expression.init <- get_expression(dataset)
  benchmark.list <- get_benchmark.dataset(expression.init, params, random_state)
  benchmark <- as.data.frame(do.call(rbind, benchmark.list))

  ground_truth <- get_ground_truth(expression.init)
  compute_metric <- function(metric) {sapply(X=benchmark$preds, FUN=get_metric, ground_truth, metric=metric)}
  for (metric in c("ARI", "NMI")) {benchmark[metric] <- compute_metric(metric)}

  benchmark["method"] <- rownames(benchmark)
  benchmark["dataset"] <- dataset
  benchmark$preds <- NULL
  for(col in c("peakRAM", "time", "ARI", "NMI")){benchmark[, col] <- as.numeric(benchmark[, col])}
  for(col in c("method", "dataset")){benchmark[, col] <- as.character(benchmark[, col])}
  return(benchmark)
}

benchmark <- compute_benchmark.dataset(dataset, params=get_default_hyperparameters(), random_state=0)
write.csv(benchmark, glue("./results/real/{dataset}.csv"), row.names=FALSE)
