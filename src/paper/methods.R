"Functions used to get the results (preds, time and peakRAM) for the benchmark.

	2024/04/02 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(scater)
  library(SeuratWrappers)
})
source("./scEVE.R")
source("./src/scEVE/trim.R")          # pre-processing
source("./src/scEVE/clusterings.R")   # individual methods

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
  preds <- preds[order(names(preds))]
  # sort cells alphabetically to facilitate benchmarking.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  time.after <- Sys.time()
  memory_summary <- gc()
  peakRAM.after <- memory_summary[11] + memory_summary[12]
  
  results <- list(peakRAM = peakRAM.after - peakRAM.before,
                  time = as.numeric(time.after - time.before, units="secs"),
                  preds = preds)
  return(results)
}

get_input.individual_method <- function(expression.init, method, n_HVGs) {
  #' Get the input specific to an individual clustering method.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param method: a valid method name. In the JOBIM paper, it is one of: 'Seurat', 'monocle3', 'CIDR' and 'SHARP'.
  #' @param n_HVGs: a numeric.
  #'
  #' @return a scRNA-seq dataset. It can be a Seurat Object, a raw-count or a log2-tpm matrix.
  #' The modality of the input is specific to each method.
  #'
  input <- get_expression.count(expression.init, n_HVGs)
  if(method %in% c("Seurat", "monocle3")){input <- get_SeurObj.count(input)}
  if(method %in% c("densityCut")){input <- log2(calculateTPM(input) + 1)}
  return(input)
}

get_benchmark.individual_method.wrapper <- function(expression.init, method, n_HVGs, random_state) {
  #' Get the results of an individual clustering method on a scRNA-seq dataset.
  #'
  #' @param expression.init: a scRNA-seq raw count matrix, without selected genes:
  #' genes are rows | cells are cols.
  #' @param method: a valid method name. In the JOBIM paper, it is one of: 'Seurat', 'monocle3', 'CIDR' and 'SHARP'.
  #' @param n_HVGs: a numeric.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #'
  input <- get_input.individual_method(expression.init, method, n_HVGs)
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
  output <- do_scEVE(expression.init, params, figures=FALSE, random_state=random_state)
  preds <- output$preds
  preds <- preds[order(names(preds))]
  # sort cells alphabetically to facilitate benchmarking.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  time.after <- Sys.time()
  memory_summary <- gc()
  peakRAM.after <- memory_summary[11] + memory_summary[12]
  
  results <- list(peakRAM = peakRAM.after - peakRAM.before,
                  time = as.numeric(time.after - time.before, units="secs"),
                  preds = preds)
  return(results)
}