"Functions used to get the results (preds, time and peakRAM) for the benchmark.

	2024/02/28 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(scater)
  library(SeuratWrappers)
})
source("./scEVE.R")
source("./src/scEVE/trim.R")          # pre-processing
source("./src/scEVE/clusterings.R")   # individual methods
source("./src/benchmark/ensemble.R")  # ensemble methods

get_benchmark.method <- function(input, fun, random_state) {
  #' Get the results of a clustering method by applying its function to a specific input.
  #' These results include:
  #' - peakRAM (Mbytes): the maximum memory usage of the method.
  #' - time (s): the computation time in seconds.
  #' - preds: a named factor, where names are cells and values are cluster labels.
  #' 
  #' @param input: a scRNA-seq dataset. It can be a Seurat Object, a raw count or a log2 tpm matrix.
  #' The modality of the input is specific to each method.
  #' @param fun: the function corresponding to the measured clustering method.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #' 
  time.before <- Sys.time()
  memory_summary <- gc(reset=TRUE)
  peakRAM.before <- memory_summary[11] + memory_summary[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  preds <- fun(input, random_state)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  time.after <- Sys.time()
  memory_summary <- gc()
  peakRAM.after <- memory_summary[11] + memory_summary[12]
  
  results <- list(peakRAM = peakRAM.after - peakRAM.before,
                  time = as.numeric(time.after - time.before, units="secs"),
                  preds = preds)
  return(results)
}

get_input.method <- function(expression.init, method, n_HVGs) {
  #' Get the modality of the input specific to the method.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param method: a valid method name, i.e. one of:
  #' 'Seurat', 'monocle3', 'CIDR', 'SHARP', 'scLCA', 'densityCut', 'scCCESS.Kmeans', 'scCCESS.SIMLR'.
  #' @param n_HVGs: a numeric.
  #'
  #' @return a scRNA-seq dataset. It can be a Seurat Object, a raw-count or a log2-tpm matrix.
  #' The modality of the input is specific to each method.
  #'
  input <- get_expression.count(expression.init, n_HVGs)
  if(method %in% c("Seurat", "monocle3")){input <- get_SeurObj.count(input)}
  if(method %in% c("scCCESS.Kmeans", "scCCESS.SIMLR", "densityCut")){input <- log2(calculateTPM(input) + 1)}
  return(input)
}

get_benchmark.method.wrapper <- function(expression.init, method, n_HVGs, random_state) {
  #' Get the results of a clustering method on a scRNA-seq dataset.
  #'
  #' @param expression.init: a scRNA-seq raw count matrix, without selected genes:
  #' genes are rows | cells are cols.
  #' @param method: a valid method name, i.e. one of:
  #' 'Seurat', 'monocle3', 'CIDR', 'SHARP', 'scLCA', 'densityCut', 'scCCESS.Kmeans', 'scCCESS.SIMLR'.
  #' @param n_HVGs: a numeric.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #'
  input <- get_input.method(expression.init, method, n_HVGs)
  fun <- get(glue("do_{method}"))
  results <- get_benchmark.method(input, fun, random_state)
  return(results)
}

get_benchmark.scEVE <- function(expression.init, params, random_state) {
  #' Get the results of scEVE by applying it with a given set of parameters on a scRNA-seq raw count matrix.
  #' These results include:
  #' - peakRAM (Mbytes): the maximum memory usage of the method.
  #' - time (s): the computation time in seconds.
  #' - preds: a named factor, where names are cells and values are cluster labels.
  #' 
  #' @param expression.init: a scRNA-seq raw-count matrix:
  #' genes are rows | cells are cols.
  #' @param params: a list of parameters.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #' 
  time.before <- Sys.time()
  memory_summary <- gc(reset=TRUE)
  peakRAM.before <- memory_summary[11] + memory_summary[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  results <- do_scEVE(expression.init, params, figures=FALSE, random_state=random_state)
  preds <- results$preds
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  time.after <- Sys.time()
  memory_summary <- gc()
  peakRAM.after <- memory_summary[11] + memory_summary[12]
  
  results <- list(peakRAM = peakRAM.after - peakRAM.before,
                  time = as.numeric(time.after - time.before, units="secs"),
                  preds = preds)
  return(results)
}

get_benchmark.scEVE.default <- function(expression.init, random_state) {
  #' Get the results of scEVE with the default leftovers strategy.
  #' 
  #' @param expression.init: a scRNA-seq raw-count matrix:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #' 
  params <- get_default_hyperparameters()
  results <- get_benchmark.scEVE(expression.init, params, random_state)
  return(results)
}

get_benchmark.scEVE.naive <- function(expression.init, random_state) {
  #' Get the results of scEVE with the naive leftovers strategy.
  #' 
  #' @param expression.init: a scRNA-seq raw-count matrix:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #' 
  params <- get_default_hyperparameters()
  params$leftovers_strategy <- "naive"
  results <- get_benchmark.scEVE(expression.init, params, random_state)
  return(results)
}

get_benchmark.scEVE.weighted <- function(expression.init, random_state) {
  #' Get the results of scEVE with the weighted leftovers strategy.
  #' 
  #' @param expression.init: a scRNA-seq raw-count matrix:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #' 
  params <- get_default_hyperparameters()
  params$leftovers_strategy <- "weighted"
  results <- get_benchmark.scEVE(expression.init, params, random_state)
  return(results)
}