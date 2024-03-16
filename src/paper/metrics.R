"Functions used to compute metrics for the benchmark.

	2024/03/07 @yanisaspic"

suppressPackageStartupMessages({
  library(aricode)
})

get_ground_truth <- function(expression.init) {
  #' Get a vector of labels corresponding to the ground truth from an expression matrix.
  #' The columns correspond to a unique cell id, and a cell id has the following structure:
  #' {ground_truth}_{i}, where i is a numeric.
  #' 
  #' @param expression.init: a scRNA-seq count matrix:
  #' genes are rows | cells are cols.
  #' 
  #' @return a vector of ground truths.
  #' 
  get_ground_truth.cell_id <- function(cell_id) {
    tmp <- strsplit(cell_id, split="_")[[1]]
    ground_truth.cell_id <- paste(head(tmp, -1), collapse="_")
    return(ground_truth.cell_id)
  }
  
  cell_ids <- colnames(expression.init)
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