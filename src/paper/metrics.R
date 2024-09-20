"Functions used to compute metrics for the benchmark.

	2024/09/19 @yanisaspic"

suppressPackageStartupMessages({
  library(aricode)  # packageVersion("aricode")==1.0.3
  library(entropy)
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
})
source("./src/scEVE/utils/misc.R")

get_ground_truth <- function(cell_ids) {
  #' Get a vector of labels corresponding to the ground truth from an expression matrix.
  #' The columns correspond to a unique cell id, and a cell id has the following structure:
  #' {ground_truth}_{i}, where i is a numeric.
  #' 
  #' @param expression.init: a vector of formatted cell ids.
  #' 
  #' @return a vector of ground truths.
  #' 
  get_ground_truth.cell_id <- function(cell_id) {
    tmp <- strsplit(cell_id, split="_")[[1]]
    ground_truth.cell_id <- paste(head(tmp, -1), collapse="_")
    return(ground_truth.cell_id)
  }
  
  ground_truth <- sapply(X=cell_ids, FUN=get_ground_truth.cell_id)

  ground_truth <- ground_truth[order(names(ground_truth))]
  # sort cells alphabetically to facilitate benchmarking.
  ground_truth <- factor(ground_truth)
  return(ground_truth)
}

get_metric <- function(preds, ground_truth, metric) {
  #' Compute a classification metric from a vector of preds and a ground truth.
  #'
  #' @param preds: a named vector: names are cell ids and values are cluster labels predicted.
  #' @param ground_truth: a named vector: names are cell ids and values are cluster labels of the authors.
  #' @param metric: a valid metric name. One of: 'ARI', 'AMI' and 'NMI'.
  #' 
  #' @return a numeric.
  #' 
  if (!identical(names(preds), names(ground_truth))) {
    stop("The cell ids of the predictions and the ground truth are not aligned.")
  }
  result <- get(metric)(preds, ground_truth)
  return(result)
}

get_distribution.population <- function(population, sheet.cells, ground_truth) {
  #' Get the distribution of ground truth labels for a given population.
  #' 
  #' @param population: a character.
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' @param ground_truth: a named vector: names are cell ids and values are cluster labels of the authors.
  #' 
  #' @return a vector of numerics.
  #' 
  cells_of_interest <- get_cells_of_interest(population, sheet.cells)
  data <- table(cells_of_interest, ground_truth[cells_of_interest])
  distribution <- apply(X=data, MARGIN=2, FUN=sum)
  return(distribution)
}

get_distributions.data <- function(sheet.cells) {
  #' Get the distributions of ground truth labels in all populations of a scEVE analysis.
  #' The output is ggplot-friendly.
  #' 
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' @param ground_truth: a named vector: names are cell ids and values are cluster labels of the authors.
  #' 
  #' @return a data.frame with three columns: 'n', 'ground_truth' and 'population'.
  #' 
  ground_truth <- get_ground_truth(rownames(sheet.cells))
  get_distributions.data.population <- function(population) {
    distribution.population <- get_distribution.population(population, sheet.cells, ground_truth)
    output <- melt(distribution.population, value.name="n")
    output$ground_truth <- rownames(output)
    output$population <- population
    return(output)
  }
  distributions.data <- lapply(X=colnames(sheet.cells), FUN=get_distributions.data.population)
  distributions.data <- do.call(rbind, distributions.data)
  return(distributions.data)
}

#' get_disruption.population <- function(distribution.population, init_entropy) {
#'   #' Get the normalized Shannon entropy of a distribution of ground truth labels.
#'   #' 
#'   #' @param distribution.population: a vector of numerics.
#'   #' @param max_entropy: the entropy of the initial population.
#'   #' 
#'   #' @return a numeric ranging from 0 to 1.
#'   #' 
#'   n <- sum(distribution.population > 0)
#'   entropy <- entropy::entropy(distribution.population, unit="log2")
#'   disruption <- entropy / init_entropy
#'   if (is.na(disruption)) {disruption <- 0}
#'     # this happens if there is a single ground truth label in the population.
#'   return(disruption)
#' }
#' 
#' get_disruption.dataset <- function(records) {
#'   #' Get the disruption and the consensus associated to each population of a dataset.
#'   #' 
#'   #' @param records: a named list of two data.frames: 'meta' and 'cells'.
#'   #' @param ground_truth: a named vector: names are cell ids and values are cluster labels of the authors.
#'   #' 
#'   #' @return a data.frame with two columns: 'consensus' and 'disruption'.
#'   #' 
#'   populations <- colnames(records$cells)
#'   cell_ids <- rownames(records$cells)
#'   ground_truth <- get_ground_truth(cell_ids)
#'   init_distribution <- table(ground_truth)
#'   init_entropy <- entropy::entropy(init_distribution)
#'   
#'   get_disruption.population.wrapper <- function(population) {
#'     distribution.population <- get_distribution.population(population, records$cells, ground_truth)
#'     disruption.population <- get_disruption.population(distribution.population, init_entropy)
#'     output <- list(disruption=disruption.population, consensus=records$meta[population, "consensus"])
#'   }
#'   
#'   disruption.dataset <- lapply(X=populations, FUN=get_disruption.population.wrapper)
#'   disruption.dataset <- do.call(rbind, disruption.dataset)
#'   return(disruption.dataset)
#' }
#' 
#' get_disruption <- function(datasets, path="./results/records") {
#'   #' Get a data.frame associating disruption and consensi values across multiple datasets.
#'   #' 
#'   #' @param datasets: a vector of dataset labels.
#'   #' @param path: a character.
#'   #' 
#'   #' @return a data.frame with three columns: 'consensus', 'disruption', and 'group'.
#'   #' 
#'   get_disruption.dataset.wrapper <- function(dataset) {
#'     records <- get_records(glue("{path}/{dataset}.xlsx"))
#'     disruption.dataset <- get_disruption.dataset(records)
#'   }
#'   disruptions <- lapply(X=datasets, FUN=get_disruption.dataset.wrapper)
#'   disruptions <- do.call(rbind, disruptions)
#'   
#'   disruptions <- as.data.frame(disruptions)
#'   for (col in c("disruption", "consensus")) {disruptions[, col] <- as.numeric(disruptions[, col])}
#'   return(disruptions)
#' }

get_similarity.individual <- function(results.individual, dataset) {
  #' Get the pairwise similarity (ARI and NMI) between clustering results 
  #' predicted by the methods integrated in scEVE.
  #' 
  #' @param results.individual: a nested list of three elements: 'peakRAM', 'time' and 'preds'.
  #' @param dataset: a character.
  #' 
  #' @return a data.frame with 5 columns: 'method_1', 'method_2', 'ARI', 'NMI' and 'dataset'.
  #' 
  combinations <- combn(names(results.individual), 2)
  combinations <- as.data.frame(t(combinations))
  colnames(combinations) <- c("method_1", "method_2")
  
  get_similarity.pairwise <- function(row, metric) {
    preds_1 <- results.individual[[row[["method_1"]]]]$preds
    preds_2 <- results.individual[[row[["method_2"]]]]$preds
    get_metric(preds_1, preds_2, metric)
  }
  
  for (metric in c("ARI", "NMI")) {
    combinations[, metric] <- apply(X=combinations, MARGIN=1,
                                    FUN=get_similarity.pairwise, metric=metric)
  }
  combinations[, "dataset"] <- dataset
  return(combinations)
}
