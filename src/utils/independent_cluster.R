"Functions called by get_independent_clusters.R to carry independent clusterings.

	2024/01/10 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(dplyr)
  library(scater)
  library(gridExtra)
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

do_Seurat <- function(SeurObj.count, random_state) {
  #' Apply a Seurat clustering algorithm.
  #' 
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() must have been applied already.
  #' @param random_state: a numeric.
  #' 
  #' @return a factor of cluster labels.
  #'
  SeurObj.count <- RunPCA(SeurObj.count,
                          features=VariableFeatures(SeurObj.count), 
                          seed.use=random_state)
  SeurObj.count <- FindNeighbors(SeurObj.count,
                                 features=VariableFeatures(SeurObj.count))
  SeurObj.count <- FindClusters(SeurObj.count,
                                random.seed=random_state)
  labels <- unname(Idents(SeurObj.count))
  return(labels)
}

do_monocle3 <- function(SeurObj.count, random_state) {
  #' Apply a monocle3 clustering algorithm.
  #' 
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() must have been applied already.
  #' @param random_state: a numeric.
  #' 
  #' @return a factor of cluster labels.
  #'
  CelDatSet <- as.cell_data_set(SeurObj.count)
  CelDatSet <- cluster_cells(CelDatSet, random_seed=random_state)
  labels <- unname(CelDatSet@clusters@listData$UMAP$clusters)
  return(labels)
}

do_SHARP <- function(expression.count, random_state) {
  #' Apply a SHARP clustering algorithm.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a vector of cluster labels.
  #'
  results <- SHARP(scExp=expression.count,
                   exp.type="count",
                   rN.seed=random_state)
  labels <- results$pred_clusters
  return(labels)
}

do_scCCESS.Kmeans <- function(log.expression.tpm, random_state) {
  #' Apply a scCCESS with K-means clustering algorithm. 
  #' Hyperparameters from Yu et al. (2022). 
  #' 
  #' @param log.expression.tpm: a scRNA-seq dataset of log2 TPM expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a vector of cluster labels.
  #'
  data <- t(log.expression.tpm)
  
  estimation <- estimate_k(
    data,
    krange=2:25,
    ensemble_sizes=20,
    cluster_func=function(x, centers) {
      kmeans(x, centers)
    },
    seed=random_state)
  
  labels <- ensemble_cluster(
    data,
    ensemble_sizes=20,
    cluster_func = function(x) {
      kmeans(x, estimation$ngroups)
    },
    seed=random_state)
  
  return(labels)
}

do_scCCESS.SIMLR <- function(log.expression.tpm, random_state) {
  #' Apply a scCCESS with SIMLR clustering algorithm. 
  #' Hyperparameters from Yu et al. (2022). 
  #' 
  #' @param log.expression.tpm: a scRNA-seq dataset of log2 TPM expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a vector of cluster labels.
  #'
  data <- t(log.expression.tpm)
  
  estimation <- estimate_k(
    data,
    krange=2:25,
    ensemble_sizes=20,
    cluster_func = function(x, centers) {
      SIMLR_Large_Scale(t(x), c=centers, kk=15)
    },
    seed=random_state)
  
  labels <- ensemble_cluster(
    data,
    ensemble_sizes=20,
    cluster_func = function(x) {
      SIMLR_Large_Scale(t(x), c=estimation$ngroups, kk=15)
    },
    seed=random_state)
  
  return(labels)
}

do_CIDR <- function(expression.count, random_state) {
  #' Apply a CIDR clustering algorithm.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a vector of cluster labels.
  #'
  set.seed(random_state)
  SinCelData <- scDataConstructor(as.matrix(expression.count))
  SinCelData <- determineDropoutCandidates(SinCelData)
  SinCelData <- wThreshold(SinCelData)
  SinCelData <- scDissim(SinCelData)
  SinCelData <- scPCA(SinCelData, plotPC=FALSE)
  SinCelData <- nPC(SinCelData)
  
  labels <- scCluster(SinCelData)@clusters
  return(labels)
}

do_densityCut <- function(log.expression.tpm, random_state) {
  #' Apply a densityCut clustering algorithm.
  #' 
  #' @param log.expression.tpm: a scRNA-seq dataset of log2 TPM expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a vector of cluster labels.
  #'
  set.seed(random_state)
  data <- t(log.expression.tpm)
  labels <- DensityCut(data, show.plot = FALSE)$cluster
  return(labels)
}

do_scLCA <- function(expression, random_state) {
  #' Apply a scLCA clustering algorithm.
  #' 
  #' @param expression: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric random seed.
  #' 
  #' @return a vector of cluster labels.
  #'
  set.seed(random_state)
  labels <- myscLCA(expression.count)[[1]]
  return(labels)
}

get_clusters <- function(expression.count, SeurObj.count, clustering_methods, random_state) {
  #' Apply multiple independent methods of clustering on a scRNA-seq dataset. 
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() must have been applied already.
  #' @param clustering_methods: a vector of valid clustering method names, i.e. in: 
  #' 'Seurat', 'monocle3', 'CIDR', 'SHARP', 'scLCA', 'densityCut', 'scCCESS.Kmeans', 'scCCESS.SIMLR'.
  #' @param random_state: a numeric.
  #' 
  #' @return a data.frame: cells are rows | methods are cols | a cell is a label.
  #' 
  log.expression.tpm <- log2(calculateTPM(expression.count) + 1)
  inputs <- list(
    Seurat=SeurObj.count,
    monocle3=SeurObj.count,
    SHARP=expression.count, # nb: SHARP is the only method explicitly requiring no trimming ?
    CIDR=expression.count,
    scLCA=expression.count,
    scCCESS.Kmeans=log.expression.tpm,
    scCCESS.SIMLR=log.expression.tpm,
    densityCut=log.expression.tpm
  )
  
  get_clusters.method <- function(method){
    fun=get(glue("do_{method}"))
    labels <- fun(inputs[[method]], random_state)
    labels <- as.character(labels)
    labels <- paste0(glue("{method}_"), labels)
    return(labels)
  }
  results <- lapply(X=clustering_methods, FUN=get_clusters.method)
  gc()
  
  clusters <- do.call(cbind, results)
  rownames(clusters) <- colnames(expression.count)
  colnames(clusters) <- clustering_methods
  return(clusters)
}

get_clusters_plots <- function(SeurObj, clusters) {
  #' Get a DimPlot with the cells colored according to the clusters for each method used.
  #'
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() must have been applied already.
  #' @param clusters: a data.frame: cells are rows | methods are cols | a cell is a label.
  #'
  #' @return a DimPlot
  #'
  get_cluster_number <- function(cluster_label){
    cluster_number <- strsplit(cluster_label, split="_")[[1]][2]
  }
  clusters <- apply(clusters, c(1,2), get_cluster_number)
  
  get_plot <- function(method){
    cells <- rownames(clusters)
    labels <- factor(clusters[, method])
    SeurObj@active.ident <- setNames(labels, cells)
    plot <- do_DimPlot(SeurObj) + ggtitle(method)
    plot <- plot + 
      theme(
        plot.title=element_text(hjust=0.5, margin=margin(1, 0, 0, 0))
      )
  }
  
  clustering_methods <- colnames(clusters)
  plots <- lapply(X=clustering_methods, FUN=get_plot)
  names(plots) <- clustering_methods
  return(plots)
}

get_independent_clusters <- function(data.loop, population, params, figures, random_state) {
  #' Conduct the independent clusters step.
  #' 
  #' @param data.loop: a list of three data.frames: 'expression.loop' and 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param population: a character.
  #' @param params: a list of parameters, with 'clustering_methods'.
  #' @param figures: a boolean. If TRUE, draw figures summarizing the independent clustering step. 
  #' @param random_state: a numeric.
  #'
  #' @return a data.frame where: rows are cells | cols are clusterings | cells are labels.
  #'
  clusters <- get_clusters(data.loop$expression.loop, 
                           data.loop$SeurObj.loop, 
                           params$clustering_methods, 
                           random_state)
  
  if (figures) {
    pdf(file = glue("./figures/{population}_independent_clusters.pdf"))
    clusters_plots <- get_clusters_plots(data.loop$SeurObj.loop,
                                         clusters)
    composite_clusters_plot <- do.call(grid.arrange, clusters_plots)
    print(composite_clusters_plot)
    dev.off()
  }
  
  return(clusters)
}
