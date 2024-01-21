"Functions used to get individual clusterings, set up for the benchmark.

	2024/01/21 @yanisaspic"

suppressPackageStartupMessages({
  library(scater)
  library(SeuratWrappers)
  
  library(cidr)
  library(scLCA)
  library(SHARP)
  library(SIMLR)
  library(Seurat)
  library(scCCESS)
  library(monocle3)
  library(densitycut)
})

get_SeurObj.count <- function(expression.init, n_HVGs) {
  #' Get a Seurat Object from a raw count expression matrix.
  #' This pre-processing step is required for Seurat and monocle3 methods.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param n_HVGs: a numeric.
  #' 
  #' @return a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() have been applied already.
  #' 
  SeurObj.count <- CreateSeuratObject(expression.init)
  SeurObj.count <- FindVariableFeatures(SeurObj.count, nfeatures=n_HVGs)
  SeurObj.count <- NormalizeData(SeurObj.count)
  SeurObj.count <- ScaleData(SeurObj.count, features=VariableFeatures(SeurObj.count))
  return(SeurObj.count)
}

get_expression.count <- function(expression.init, n_HVGs) {
  #' Get a raw count expression matrix with HVGs only.
  #' This pre-processing step is required for SHARP, CIDR, scLCA, scCCESS.Kmeans, scCCESS.SIMLR, densityCut.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param n_HVGs: a numeric.
  #' 
  #' @return a scRNA-seq dataset of raw count expression, with HVGs only:
  #' genes are rows | cells are cols.
  #' 
  SeurObj.count <- CreateSeuratObject(expression.init)
  SeurObj.count <- FindVariableFeatures(SeurObj.count, nfeatures=n_HVGs)
  expression.count <- expression.init[VariableFeatures(SeurObj.count),]
  return(expression.count)
}

benchmark_Seurat <- function(expression.init, n_HVGs, random_state) {
  #' Apply a Seurat clustering algorithm for the benchmark.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param n_HVGs: a numeric.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'time', 'memory' and 'labels'.
  #'
  start_time <- Sys.time()
  max_memory_used.default <- gc(reset=TRUE)[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # preprocessing
  ###############
  SeurObj.count <- get_SeurObj.count(expression.init, n_HVGs)
  SeurObj.count <- RunPCA(SeurObj.count,
                          features=VariableFeatures(SeurObj.count), 
                          seed.use=random_state)
  SeurObj.count <- FindNeighbors(SeurObj.count,
                                 features=VariableFeatures(SeurObj.count))
  # clustering
  ############
  SeurObj.count <- FindClusters(SeurObj.count,
                                random.seed=random_state)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  labels <- unname(Idents(SeurObj.count))
  results <- list(time=Sys.time()-start_time,
                  memory=gc()[12] - max_memory_used.default,
                  labels=labels)
  return(results)
}

benchmark_monocle3 <- function(expression.init, n_HVGs, random_state) {
  #' Apply a monocle3 clustering algorithm for the benchmark.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param n_HVGs: a numeric.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'time', 'memory' and 'labels'.
  #'
  start_time <- Sys.time()
  max_memory_used.default <- gc(reset=TRUE)[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # preprocessing
  ###############
  SeurObj.count <- get_SeurObj.count(expression.init, n_HVGs)
  SeurObj.count <- RunUMAP(SeurObj.count,
                           features=VariableFeatures(SeurObj.count),
                           seedd.use=random_state)
  CelDatSet <- as.cell_data_set(SeurObj.count)
  # clustering
  ############
  CelDatSet <- cluster_cells(CelDatSet, random_seed=random_state)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  labels <- unname(CelDatSet@clusters@listData$UMAP$clusters)
  results <- list(time=Sys.time()-start_time,
                  memory=gc()[12] - max_memory_used.default,
                  labels=labels)
  return(results)
}

benchmark_SHARP <- function(expression.init, n_HVGs, random_state) {
  #' Apply a SHARP clustering algorithm for the benchmark.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param n_HVGs: a numeric.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'time', 'memory' and 'labels'.
  #'
  start_time <- Sys.time()
  max_memory_used.default <- gc(reset=TRUE)[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # preprocessing
  ###############
  expression.count <- get_expression.count(expression.init, n_HVGs)
  # clustering
  ############
  labels <- SHARP(scExp=expression.count,
                  exp.type="count",
                  rN.seed=random_state)  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  labels <- results$pred_clusters
  results <- list(time=Sys.time()-start_time,
                  memory=gc()[12] - max_memory_used.default,
                  labels=labels)
  return(results)
}

benchmark_densityCut <- function(expression.init, n_HVGs, random_state) {
  #' Apply a SHARP clustering algorithm for the benchmark.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param n_HVGs: a numeric.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'time', 'memory' and 'labels'.
  #'
  start_time <- Sys.time()
  max_memory_used.default <- gc(reset=TRUE)[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # preprocessing
  ###############
  expression.count <- get_expression.count(expression.init, n_HVGs)
  log.expression.tpm <- log2(calculateTPM(expression.count) + 1)
  # clustering
  ############
  set.seed(random_state)
  data <- t(log.expression.tpm)
  results <- DensityCut(data, show.plot = FALSE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  labels <- results$cluster
  results <- list(time=Sys.time()-start_time,
                  memory=gc()[12] - max_memory_used.default,
                  labels=labels)
  return(results)
}
