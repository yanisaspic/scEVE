"Functions used to get the results (preds, time and peakRAM) for the benchmark.

	2024/09/19 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(scater)
  library(SeuratWrappers)
})
source("./scEVE.R")
source("./src/scEVE/trim.R")          # pre-processing
source("./src/scEVE/clusterings.R")   # individual methods

get_results.method <- function(input, fun, random_state) {
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

get_results.individual_method.wrapper <- function(method, expression.init, n_HVGs, random_state) {
  #' Get the results of an individual clustering method on a scRNA-seq dataset.
  #'
  #' @param method: a valid method name. In the JOBIM paper, it is one of: 'Seurat', 'monocle3', 'CIDR' and 'SHARP'.
  #' @param expression.init: a scRNA-seq raw count matrix, without selected genes:
  #' genes are rows | cells are cols.
  #' @param n_HVGs: a numeric.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #'
  input <- get_input.individual_method(expression.init, method, n_HVGs)
  fun <- get(glue("do_{method}"))
  results <- get_results.method(input, fun, random_state)
  return(results)
}

get_results.scEVE <- function(expression.init, params, random_state, save, figures) {
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
  #' @param save: a boolean. If TRUE, the results are saved in a records file.
  #' @param figures: a boolean. If TRUE, intermediate figures are saved.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #' 
  time.before <- Sys.time()
  memory_summary <- gc(reset=TRUE)
  peakRAM.before <- memory_summary[11] + memory_summary[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output <- do_scEVE(expression.init, params, random_state=random_state,
                     save=save, figures=figures)
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

get_results.individual <- function(expression.init, params, random_state, save, figures) {
  #' Get the results of the individual clustering methods on a scRNA-seq dataset.
  #' These results include:
  #' - peakRAM (Mbytes): the maximum memory usage of the method.
  #' - time (s): the computation time in seconds.
  #' - preds: a named factor, where names are cells and values are cluster labels.
  #'
  #' @param expression.init: a scRNA-seq raw-count matrix:
  #' genes are rows | cells are cols.
  #' @param params: a list of parameters.
  #' @param random_state: a numeric.
  #' @param save: a boolean. If TRUE, the results are saved in a records file.
  #' @param figures: a boolean. If TRUE, intermediate figures are saved.
  #' 
  #' @return a nested list of three elements: 'peakRAM', 'time' and 'preds'.
  #' 
  individual_methods <- params$clustering_methods
  results <- lapply(X=individual_methods,
                    FUN=get_results.individual_method.wrapper,
                    expression.init=expression.init,
                    n_HVGs=params$n_HVGs,
                    random_state=random_state)
  names(results) <- c(individual_methods)
  return(results)
}

get_benchmark.individual <- function(expression.init, ground_truth, dataset, params, random_state) {
  #' Get the results of the individual clustering methods it uses on a scRNA-seq dataset.
  #' These results include:
  #' - peakRAM (Mbytes): the maximum memory usage of the method.
  #' - time (s): the computation time in seconds.
  #' - ARI: a clustering performance metric.
  #' - NMI: a clustering performance metric.
  #'
  #' @param expression.init: a scRNA-seq count matrix where rows are genes | cells are cols
  #' @param ground_truth: a vector of ground truths.
  #' @param dataset: a character.
  #' @param params: a list of parameters.
  #' @param random_state: a numeric.
  #'
  #' @return a data.frame with 6 columns: 'peakRAM', 'time', 'ARI', 'NMI', 'method' and 'dataset'.
  #'
  results.individual <- get_results.individual(expression.init, params, random_state,
                                               save, figures)
  benchmark <- as.data.frame(do.call(rbind, results.individual))
  
  compute_metric.preds <- function(preds, metric) {get_metric(preds, ground_truth[names(preds)], metric)}
    # some individual methods (e.g. SHARP, densityCut) won't cluster cells without HVGs expressed;
    # we only compare cells clustered by a method to prevent computational errors and ensure a fair comparison.
  compute_metric <- function(metric) {sapply(X=benchmark$preds, FUN=compute_metric.preds, metric=metric)}
  for (metric in c("ARI", "NMI")) {benchmark[metric] <- compute_metric(metric)}
  
  benchmark["method"] <- rownames(benchmark)
  benchmark["dataset"] <- dataset
  benchmark$preds <- NULL
  for(col in c("peakRAM", "time", "ARI", "NMI")){benchmark[, col] <- as.numeric(benchmark[, col])}
  for(col in c("method", "dataset")){benchmark[, col] <- as.character(benchmark[, col])}
  return(benchmark)
}

get_benchmark <- function(expression.init, ground_truth, dataset, params,
                          random_state, save=FALSE, figures=FALSE) {
  #' Get the results of scEVE and the individual clustering methods it uses on a scRNA-seq dataset.
  #' These results include:
  #' - peakRAM (Mbytes): the maximum memory usage of the method.
  #' - time (s): the computation time in seconds.
  #' - ARI: a clustering performance metric.
  #' - NMI: a clustering performance metric.
  #'
  #' @param expression.init: a scRNA-seq count matrix where rows are genes | cells are cols
  #' @param ground_truth: a vector of ground truths.
  #' @param dataset: a character.
  #' @param params: a list of parameters.
  #' @param random_state: a numeric.
  #' @param save: a boolean. If TRUE, the results of scEVE are saved in a records file.
  #' @param figures: a boolean. If TRUE, intermediate figures are saved.
  #'
  #' @return a data.frame with 6 columns: 'peakRAM', 'time', 'ARI', 'NMI', 'method' and 'dataset'.
  #'
  format_benchmark <- function(benchmark) {
    for (col in c("peakRAM", "time", "ARI", "NMI")) {benchmark[, col] <- as.numeric(benchmark[, col])}
    for (col in c("method", "dataset")) {benchmark[, col] <- as.character(benchmark[, col])}
    return(benchmark)
  }
  
  benchmark <- get_benchmark.individual(expression.init, ground_truth,
                                                   dataset, params, random_state)
  benchmark <- format_benchmark(benchmark)
  write.csv(benchmark, glue("./benchmark/{dataset}.csv"), row.names=FALSE)
    # save the individual methods performances early in case scEVE times out (e.g. Peng_HumPDAC)
  
  results.scEVE <- get_results.scEVE(expression.init, params, random_state, save, figures)
  for (metric in c("ARI", "NMI")) {
    results.scEVE[[metric]] <- get_metric(results.scEVE$preds, ground_truth, metric)}
  results.scEVE$method <- "scEVE"
  results.scEVE$dataset <- dataset
  results.scEVE$preds <- NULL
  
  row.scEVE <- rbind(results.scEVE)
  benchmark <- rbind(benchmark, row.scEVE)
  benchmark <- format_benchmark(benchmark)
  return(benchmark)
}