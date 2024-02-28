"Functions used to get individual clusterings, set up for the benchmark.

	2024/02/28 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(scater)
  library(SeuratWrappers)
  
  library(SAFEclustering)
  library(SAMEclustering)
  library(clusterExperiment)
})
source("./src/scEVE/trim.R")          # pre-processing
source("./src/scEVE/clusterings.R")   # individual methods
source("./src/benchmark/ensemble.R")  # ensemble methods

get_benchmark.method <- function(input, fun, random_state) {
  #' Get the results of a clustering method by applying its function to a specific input.
  #' These results include:
  #' - peakRAM (Mo): the maximum memory usage of the method.
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
  #' @param expression.init: a scRNA-seq raw count matrix, without selected genes.
  #' @param method: a valid method name, i.e. one of:
  #' 'Seurat', 'monocle3', 'CIDR', 'SHARP', 'scLCA', 'densityCut', 'scCCESS.Kmeans', 'scCCESS.SIMLR'.
  #' @param n_HVGs: a numeric.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #'
  input <- get_input.method(expression.init, method, n_HVGs)
  fun <- glue("do_{method}")
  results <- get_results.method(input, fun, random_state)
  return(results)
}

###########################
#     E N S E M B L E     #
#      M E T H O D S      #
###########################
setup_SAxE <- function(clusterings) {
  #' Set up a table of clusterings for SAFE and SAME methods.
  #' 
  #' @param clusterings: a data.frame where: rows are cells | cols are methods | cells are labels.
  #' 
  #' @return a matrix where: rows are methods | cols are cells | cells are a numeric.
  #' 
  setup_SAxE.col <- function(col) {as.numeric(as.factor(col))}
  result <- apply(X=clusterings, MARGIN=2, FUN=setup_SAxE.col)
  return(result)
}

benchmark_SAFE <- function(expression.count, SeurObj.count, random_state,
                           clustering_methods=c()) {
  #' Apply a SAFE clustering algorithm for the benchmark, with the default individual methods.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, with HVGs only:
  #' genes are rows | cells are cols.
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() have been applied already.
  #' @param random_state: a numeric.
  #' @param clustering_methods: a vector of valid method names. If empty, uses default SAFE methods.
  #' 
  #' @return a list of three elements: 'time', 'peakRAM' and 'labels'.
  #'
  start_time <- Sys.time()
  memory_data <- gc(reset=TRUE)
  max_memory_used.default <- memory_data[11] + memory_data[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (length(clustering_methods) == 0) {
    clusterings <- individual_clustering(expression.count)
    }
  else{
    SeurObj.count <- RunUMAP(SeurObj.count,
                             features=VariableFeatures(SeurObj.count),
                             seed.use=random_state)
    clusterings <- apply_clustering_algorithms(expression.count, SeurObj.count,
                                               clustering_methods, random_state)
    clusterings <- setup_SAxE(clusterings)
  }
  cluster.ensemble <- SAFE(clusterings)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  memory_data <- gc()
  results <- list(time=as.numeric(Sys.time()-start_time, units="secs"),
                  peakRAM=(memory_data[11] + memory_data[12]) - max_memory_used.default,
                  labels=labels)
  return(results)
}

benchmark_SAME <- function(expression.count, SeurObj.count, random_state,
                           clustering_methods=c(), criterion="BIC") {
  #' Apply a SAFE clustering algorithm for the benchmark, with the default individual methods.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, with HVGs only:
  #' genes are rows | cells are cols.
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() have been applied already.
  #' @param random_state: a numeric.
  #' @param clustering_methods: a vector of valid method names. If empty, uses default SAFE methods.
  #' @param criterion: a character indicating if the best clustering corresponds to the best 'AIC' or 'BIC'.
  #' 
  #' @return a list of three elements: 'time', 'peakRAM' and 'labels'.
  #'
  start_time <- Sys.time()
  memory_data <- gc(reset=TRUE)
  max_memory_used.default <- memory_data[11] + memory_data[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (length(clustering_methods) == 0) {
    clusterings <- individual_clustering(expression.count)
  }
  else{
    SeurObj.count <- RunUMAP(SeurObj.count,
                             features=VariableFeatures(SeurObj.count),
                             seed.use=random_state)
    clusterings <- apply_clustering_algorithms(expression.count, SeurObj.count,
                                               clustering_methods, random_state)
  }
  clusterings <- setup_SAxE(clusterings)
  cluster.ensemble <- SAMEclustering(Y=t(clusterings), rep=3, SEED=random_state)
  solution <- glue("{criterion}cluster")
  labels <- cluster.ensemble[[solution]]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  memory_data <- gc()
  results <- list(time=as.numeric(Sys.time()-start_time, units="secs"),
                  peakRAM=(memory_data[11] + memory_data[12]) - max_memory_used.default,
                  labels=labels)
  return(results)
}