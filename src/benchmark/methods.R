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
  
  library(SAFEclustering)
  library(SAMEclustering)
  library(clusterExperiment)
})
source("./src/scEVE/trim.R")        # get_n_HVGs
source("./src/scEVE/clusterings.R") # individual methods

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
  HVGs <- get_n_HVGs(expression.init, n=n_HVGs)
  expression.count <- expression.init[HVGs,]
  return(expression.count)
}

get_SeurObj.count <- function(expression.count) {
  #' Get a Seurat Object from a raw count expression matrix.
  #' This pre-processing step is required for Seurat and monocle3 methods.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, with HVGs only:
  #' genes are rows | cells are cols.
  #' @param n_HVGs: a numeric.
  #' 
  #' @return a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() have been applied already.
  #' 
  SeurObj.count <- CreateSeuratObject(expression.count)
  VariableFeatures(SeurObj.count) <- rownames(expression.count)
  SeurObj.count <- NormalizeData(SeurObj.count)
  SeurObj.count <- ScaleData(SeurObj.count,
                             features=VariableFeatures(SeurObj.count))
  return(SeurObj.count)
}
get_labels <- function(results){results$labels}
save_time <- function(results){rep(results$time, length(results$labels))}
save_peakRAM <- function(results){rep(results$peakRAM, length(results$labels))}



###########################
#   I N D I V I D U A L   #
#      M E T H O D S      #
###########################
benchmark_Seurat <- function(SeurObj.count, random_state) {
  #' Apply a Seurat clustering algorithm for the benchmark.
  #' 
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() have been applied already.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'time', 'peakRAM' and 'labels'.
  #'
  start_time <- Sys.time()
  memory_data <- gc(reset=TRUE)
  max_memory_used.default <- memory_data[11] + memory_data[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SeurObj.count <- RunPCA(SeurObj.count,
                          features=VariableFeatures(SeurObj.count), 
                          seed.use=random_state)
  SeurObj.count <- FindNeighbors(SeurObj.count,
                                 features=VariableFeatures(SeurObj.count))
  SeurObj.count <- FindClusters(SeurObj.count,
                                random.seed=random_state)
  labels <- unname(Idents(SeurObj.count))
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  memory_data <- gc()
  results <- list(time=as.numeric(Sys.time()-start_time, units="secs"),
                  peakRAM=(memory_data[11] + memory_data[12]) - max_memory_used.default,
                  labels=labels)
  return(results)
}

benchmark_monocle3 <- function(SeurObj.count, random_state) {
  #' Apply a monocle3 clustering algorithm for the benchmark.
  #' 
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() have been applied already.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'time', 'peakRAM' and 'labels'.
  #'
  start_time <- Sys.time()
  memory_data <- gc(reset=TRUE)
  max_memory_used.default <- memory_data[11] + memory_data[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SeurObj.count <- RunUMAP(SeurObj.count,
                           features=VariableFeatures(SeurObj.count),
                           seedd.use=random_state)
  CelDatSet <- as.cell_data_set(SeurObj.count)
  CelDatSet <- cluster_cells(CelDatSet, random_seed=random_state)
  labels <- unname(CelDatSet@clusters@listData$UMAP$clusters)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  memory_data <- gc()
  results <- list(time=as.numeric(Sys.time()-start_time, units="secs"),
                  peakRAM=(memory_data[11] + memory_data[12]) - max_memory_used.default,
                  labels=labels)
  return(results)
}

benchmark_SHARP <- function(expression.count, random_state) {
  #' Apply a SHARP clustering algorithm for the benchmark.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, with HVGs only:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'time', 'peakRAM' and 'labels'.
  #'
  start_time <- Sys.time()
  memory_data <- gc(reset=TRUE)
  max_memory_used.default <- memory_data[11] + memory_data[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  results <- SHARP(scExp=expression.count,
                  exp.type="count",
                  rN.seed=random_state)  
  labels <- results$pred_clusters
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  memory_data <- gc()
  results <- list(time=as.numeric(Sys.time()-start_time, units="secs"),
                  peakRAM=(memory_data[11] + memory_data[12]) - max_memory_used.default,
                  labels=labels)
  return(results)
}

benchmark_densityCut <- function(expression.count, random_state) {
  #' Apply a SHARP clustering algorithm for the benchmark.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, with HVGs only:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'time', 'peakRAM' and 'labels'.
  #'
  start_time <- Sys.time()
  memory_data <- gc(reset=TRUE)
  max_memory_used.default <- memory_data[11] + memory_data[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  set.seed(random_state)
  log.expression.tpm <- log2(calculateTPM(expression.count) + 1)
  data <- t(log.expression.tpm)
  results <- DensityCut(data, show.plot = FALSE)
  labels <- results$cluster
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  memory_data <- gc()
  results <- list(time=as.numeric(Sys.time()-start_time, units="secs"),
                  peakRAM=(memory_data[11] + memory_data[12]) - max_memory_used.default,
                  labels=labels)
  return(results)
}



###########################
#     E N S E M B L E     #
#      M E T H O D S      #
###########################
benchmark_scEVE <- function(expression.init, params, random_state) {
  #' Apply a SHARP clustering algorithm for the benchmark.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param params: a list of parameters.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'time', 'peakRAM' and 'labels'.
  #'
  start_time <- Sys.time()
  memory_data <- gc(reset=TRUE)
  max_memory_used.default <- memory_data[11] + memory_data[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  records <- do_scEVE(expression.init,
                      params=params,
                      figures=FALSE,
                      random_state=random_state)
  classification <- get_classification(records)
  labels <- classification$pred
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  memory_data <- gc()
  results <- list(time=as.numeric(Sys.time()-start_time, units="secs"),
                  peakRAM=(memory_data[11] + memory_data[12]) - max_memory_used.default,
                  labels=labels)
  return(results)
}

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