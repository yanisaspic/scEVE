"Functions called by get_independent_clusters.R to carry independent clusterings.

	2024/01/10 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(dplyr)
  
  library(cidr)
  library(scLCA)
  library(SHARP)
  library(SIMLR)
  library(Seurat)
  library(scater)
  library(scCCESS)
  library(monocle3)
  library(densitycut)
})

do_Seurat <- function(SeurObj.count, seed) {
  #' Apply a Seurat clustering algorithm.
  #' 
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() must have been applied already.
  #' @param seed: a numeric random seed.
  #' 
  #' @return a factor of cluster labels.
  #'
  SeurObj.count <- RunPCA(SeurObj.count,
                    features=VariableFeatures(SeurObj.count), 
                    seed.use=seed)
  SeurObj.count <- FindNeighbors(SeurObj.count,
                           features=VariableFeatures(SeurObj.count))
  SeurObj.count <- FindClusters(SeurObj.count,
                          random.seed=seed)
  labels <- unname(Idents(SeurObj.count))
  return(labels)
}

do_monocle3 <- function(SeurObj.count, seed) {
  #' Apply a monocle3 clustering algorithm.
  #' 
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() must have been applied already.
  #' @param seed: a numeric random seed.
  #' 
  #' @return a factor of cluster labels.
  #'
  CelDatSet <- as.cell_data_set(SeurObj.count)
  CelDatSet <- cluster_cells(CelDatSet, random_seed=seed)
  labels <- unname(CelDatSet@clusters@listData$UMAP$clusters)
  return(labels)
}

do_SHARP <- function(expression.count, seed) {
  #' Apply a SHARP clustering algorithm.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param seed: a numeric random seed.
  #' 
  #' @return a vector of cluster labels.
  #'
  results <- SHARP(scExp=expression.count,
                   exp.type="count",
                   rN.seed=seed)
  labels <- results$pred_clusters
  return(labels)
}

do_scCCESS.Kmeans <- function(log.expression.tpm, seed) {
  #' Apply a scCCESS with K-means clustering algorithm. 
  #' Hyperparameters from Yu et al. (2022). 
  #' 
  #' @param log.expression.tpm: a scRNA-seq dataset of log2 TPM expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param seed: a numeric random seed.
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
    seed=SEED)
  
  labels <- ensemble_cluster(
    data,
    ensemble_sizes=20,
    cluster_func = function(x) {
      kmeans(x, estimation$ngroups)
      },
    seed=SEED)
  
  return(labels)
}

do_scCCESS.SIMLR <- function(log.expression.tpm, seed) {
  #' Apply a scCCESS with SIMLR clustering algorithm. 
  #' Hyperparameters from Yu et al. (2022). 
  #' 
  #' @param log.expression.tpm: a scRNA-seq dataset of log2 TPM expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param seed: a numeric random seed.
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
    seed=SEED)
  
  labels <- ensemble_cluster(
    data,
    ensemble_sizes=20,
    cluster_func = function(x) {
      SIMLR_Large_Scale(t(x), c=estimation$ngroups, kk=15)
      },
    seed=SEED)
  
  return(labels)
}

do_CIDR <- function(expression.count, seed=NA) {
  #' Apply a CIDR clustering algorithm.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param seed: a numeric random seed.
  #' 
  #' @return a vector of cluster labels.
  #'
  SinCelData <- scDataConstructor(as.matrix(expression.count))
  SinCelData <- determineDropoutCandidates(SinCelData)
  SinCelData <- wThreshold(SinCelData)
  SinCelData <- scDissim(SinCelData)
  SinCelData <- scPCA(SinCelData, plotPC=FALSE)
  SinCelData <- nPC(SinCelData)
  
  labels <- scCluster(SinCelData)@clusters
  return(labels)
}

do_densityCut <- function(log.expression.tpm, seed=NA) {
  #' Apply a densityCut clustering algorithm.
  #' 
  #' @param log.expression.tpm: a scRNA-seq dataset of log2 TPM expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param seed: a numeric random seed.
  #' 
  #' @return a vector of cluster labels.
  #'
  data <- t(log.expression.tpm)
  labels <- DensityCut(data, show.plot = FALSE)$cluster
  return(labels)
}

do_scLCA <- function(expression, seed=NA) {
  #' Apply a scLCA clustering algorithm.
  #' 
  #' @param expression: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param seed: a numeric random seed.
  #' 
  #' @return a vector of cluster labels.
  #'
  labels <- myscLCA(expression.count)[[1]]
  return(labels)
}

get_clusters <- function(expression.count, SeurObj.count, methods, seed) {
  #' Apply independently multiple methods of clustering algorithms on a scRNA-seq dataset. 
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() must have been applied already.
  #' @param methods: a vector of valid clustering method names, i.e. in: 
  #' 'Seurat', 'monocle3', 'CIDR', 'SHARP', 'scLCA', 'densityCut', 'scCCESS.Kmeans', 'scCCESS.SIMLR'.
  #' @param seed: a numeric random seed.
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
    labels <- fun(inputs[[method]], seed)
    labels <- as.character(labels)
    labels <- paste0(glue("{method}_"), labels)
    return(labels)
  }
  results <- lapply(X=methods, FUN=get_clusters.method)
  gc()
  
  clusters <- do.call(cbind, results)
  rownames(clusters) <- colnames(expression.count)
  colnames(clusters) <- methods
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
  
  methods <- colnames(clusters)
  plots <- lapply(X=methods, FUN=get_plot)
  names(plots) <- methods
  return(plots)
}

get_independent_clusters <- function(expression.count, SeurObj.count, methods, seed=0, draw=FALSE) {
  #' Conduct the independent clusters step.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() must have been applied already.
  #' @param methods: a vector of valid clustering method names, i.e. in: 
  #' 'Seurat', 'monocle3', 'CIDR', 'SHARP', 'scLCA', 'densityCut', 'scCCESS.Kmeans', 'scCCESS.SIMLR'.
  #' @param seed: a numeric random seed.
  #' @param draw: a boolean. If TRUE, a summary figure is drawn and saved.
  #' 
  clusters <- get_clusters(expression.count, SeurObj.count, methods, seed)
  output <- list(clusters=clusters, step="seeds")
  
  if (draw) {
    pdf(file = glue("./results/figures/{POPULATION}_independent_clusters.pdf"))
    clusters_plots <- get_clusters_plots(SeurObj.count,
                                         clusters)
    composite_clusters_plot <- do.call(grid.arrange, clusters_plots)
    print(composite_clusters_plot)
    dev.off()
  }
  
  return(output)
}