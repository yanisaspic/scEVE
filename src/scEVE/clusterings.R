"Functions used to get individual clusterings.
  In the scEVE JOBIM paper, 4 methods are used: Seurat, monocle3, SHARP and densityCut.

	2024/04/02 @yanisaspic"

suppressPackageStartupMessages({
  library(densitycut)
  library(dplyr)
  library(glue)
  library(gridExtra) # packageVersion("gridExtra")==2.3
  library(monocle3)
  library(scater)
  library(Seurat)
  library(SeuratWrappers)
  library(SHARP)
})

format_preds <- function(cells, labels) {
  #' Format the predictions of a method to obtain a standard output.
  #' Seurat and monocle3 are already formatted, so the function is not required.
  #' 
  #' @param data: a vector of cells names.
  #' @param labels: a vector of cluster labels.
  #' 
  #' @return a named factor, where names are cells and values are cluster labels.
  #' 
  preds <- setNames(labels, cells)
  preds <- factor(preds)
  return(preds)
}

do_Seurat <- function(SeurObj.count, random_state) {
  #' Apply a Seurat clustering algorithm.
  #' 
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' VariableFeatures() and ScaleData() must have been applied already.
  #' @param random_state: a numeric.
  #' 
  #' @return a named factor, where names are cells and values are cluster labels.
  #'
  SeurObj.count <- RunPCA(SeurObj.count,
                          features=VariableFeatures(SeurObj.count), 
                          seed.use=random_state)
  SeurObj.count <- FindNeighbors(SeurObj.count,
                                 features=VariableFeatures(SeurObj.count))
  SeurObj.count <- FindClusters(SeurObj.count,
                                random.seed=random_state)
  preds <- Idents(SeurObj.count)
  return(preds)
}

do_monocle3 <- function(SeurObj.count, random_state) {
  #' Apply a monocle3 clustering algorithm.
  #' 
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' VariableFeatures() and ScaleData() must have been applied already.
  #' @param random_state: a numeric.
  #' 
  #' @return a named factor, where names are cells and values are cluster labels.
  #'
  if (!"umap" %in% names(SeurObj.count@reductions)) {
    SeurObj.count <- RunUMAP(SeurObj.count, features=VariableFeatures(SeurObj.count),
                             seed.use=random_state) 
  }
  CelDatSet <- as.cell_data_set(SeurObj.count)
  CelDatSet <- cluster_cells(CelDatSet, random_seed=random_state)
  preds <- CelDatSet@clusters@listData$UMAP$clusters
  return(preds)
}

do_SHARP <- function(expression.count, random_state) {
  #' Apply a SHARP clustering algorithm.
  #' 
  #' @param expression.count: a scRNA-seq dataset of raw count expression, with selected genes:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a named factor, where names are cells and values are cluster labels.
  #'
  results <- SHARP(scExp=expression.count,
                   exp.type="count",
                   n.cores = 1,
                   rN.seed=random_state)
  preds <- format_preds(cells=colnames(expression.count), labels=results$pred_clusters)
  return(preds)
}

do_densityCut <- function(log.expression.tpm, random_state) {
  #' Apply a densityCut clustering algorithm.
  #' 
  #' @param log.expression.tpm: a scRNA-seq dataset of log2 TPM expression, with selected genes:
  #' genes are rows | cells are cols.
  #' @param random_state: a numeric.
  #' 
  #' @return a named factor, where names are cells and values are cluster labels.
  #'
  set.seed(random_state)
  data <- t(log.expression.tpm)
  results <- DensityCut(data, show.plot = FALSE)
  preds <- format_preds(cells=colnames(log.expression.tpm), labels=results$cluster)
  return(preds)
}

apply_clustering_algorithms <- function(data.loop, clustering_methods, random_state) {
  #' Apply multiple independent methods of clustering on a scRNA-seq dataset. 
  #' 
  #' @param data.loop: a list of three data.frames: 'expression.loop', 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param clustering_methods: a vector of valid clustering method names.
  #' In the scEVE JOBIM paper, 4 methods are used: 'Seurat', 'monocle3', 'SHARP' and 'densityCut'.
  #' @param random_state: a numeric.
  #' 
  #' @return a data.frame: cells are rows | methods are cols | a cell is a label.
  #' 
  log.expression.tpm <- log2(calculateTPM(data.loop$expression.loop) + 1)
  inputs <- list(
    Seurat=data.loop$SeurObj.loop, monocle3=data.loop$SeurObj.loop,
    SHARP=data.loop$expression.loop,
    densityCut=log.expression.tpm
  )
  
  get_clusterings.method <- function(method){
    fun=get(glue("do_{method}"))
    preds <- fun(inputs[[method]], random_state)
    return(preds)
  }
  clusterings <- sapply(X=clustering_methods, FUN=get_clusterings.method)
  gc()
  
  rename.col <- function(method) {
    rename.pred <- function(pred){glue("{method}_{pred}")}
    col <- sapply(X=clusterings[, method], FUN=rename.pred)
    return(col)
  }
  clusterings <- sapply(X=clustering_methods, FUN=rename.col)
  return(clusterings)
}

get_clusterings_plots <- function(SeurObj, clusterings) {
  #' Get a DimPlot with the cells colored according to the clusterings for each method used.
  #'
  #' @param SeurObj.count: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() must have been applied already.
  #' @param clusterings: a data.frame: cells are rows | methods are cols | a cell is a label.
  #'
  #' @return a DimPlot
  #'
  get_cluster_number <- function(cluster_label){
    cluster_number <- strsplit(cluster_label, split="_")[[1]][2]
  }
  clusterings <- apply(clusterings, c(1,2), get_cluster_number)
  
  get_plot <- function(method) {
    cells <- rownames(clusterings)
    labels <- factor(clusterings[, method])
    SeurObj@active.ident <- setNames(labels, cells)
    plot <- do_DimPlot(SeurObj) + ggtitle(method)
    plot <- plot + theme(plot.title=element_text(hjust=0.5, margin=margin(1, 0, 0, 0)))
    return(plot)
  }
  
  clustering_methods <- colnames(clusterings)
  plots <- lapply(X=clustering_methods, FUN=get_plot)
  names(plots) <- clustering_methods
  return(plots)
}

get_clusterings <- function(data.loop, population, params, figures, random_state) {
  #' Get individual clusterings.
  #' 
  #' @param data.loop: a list of three data.frames: 'expression.loop' and 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param population: a character.
  #' @param params: a list of parameters, with 'clustering_methods'.
  #' @param figures: a boolean. If TRUE, draw figures summarizing the independent clustering step. 
  #' @param random_state: a numeric.
  #'
  #' @return a data.frame where: rows are cells | cols are clusterings | cells are labels.
  #'
  clusterings <- apply_clustering_algorithms(data.loop, params$clustering_methods, random_state)
  if (figures) {
    pdf(file = glue("./figures/{population}_clusterings.pdf"))
    clusterings_plots <- get_clusterings_plots(data.loop$SeurObj.loop, clusterings)
    composite_clusterings_plot <- do.call(grid.arrange, clusterings_plots)
    print(composite_clusterings_plot)
    dev.off()
  }
  return(clusterings)
}
