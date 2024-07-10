"Functions called to trim the data specific to an iteration.

	2024/04/02 @yanisaspic"

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(qpdf)
  library(SCpubr)
  library(Seurat)
})
source("./src/scEVE/utils/misc.R")

init_records <- function(expression.init) {
  #' Get a named list of three data.frames, with: 
  #' - 'cells': rows are cells | cols are populations. 
  #' The value of the cell i,j corresponds to the probability that a cell i belongs to a population j.
  #' - 'markers': rows are genes | cols are populations. 
  #' The value of the cell i,j is 1 if a cell is a marker of a population. Else, it is 0.
  #' - 'meta': rows are populations | cols are 'consensus', 'parent', 'n' and 'to_do'.
  #' 
  #' @param expression.init: the base scRNA-seq dataset: rows are genes | cols are cells.
  #' 
  #' @return a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' 
  cells <- data.frame(
    C=as.numeric(rep(1, ncol(expression.init))),
    row.names=colnames(expression.init)
    )
  markers <- data.frame(
    C=as.numeric(rep(0, nrow(expression.init))),
    row.names=rownames(expression.init)
    )
  meta <- data.frame(
    consensus=0,
    parent="?",
    n=as.numeric(ncol(expression.init)),
    to_dig=as.numeric(1),
    row.names="C"
    )
  records <- list(cells=cells, markers=markers, meta=meta)
  return(records)
}

get_undug_population <- function(records) {
  #' Get a population label that has not been analyzed yet.
  #' This function is required to iterate the algorithm (scEDEA-L).
  #'
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' 
  #' @return a character.
  #' 
  meta <- records$meta
  undug_populations <- rownames(meta[meta$to_dig==TRUE, ])
  population <- undug_populations[1]
  return(population)
}

get_n_HVGs <- function(expression, n) {
  #' Get the n most variable genes in the data, by using Seurat.
  #'
  #' @param expression: a scRNA-seq dataset: rows are genes | cols are cells.
  #' @param n: the number of highly variable genes to sample.
  #'
  #' @return a vector of n genes
  #'
  SeurObj <- CreateSeuratObject(expression)
  SeurObj <- FindVariableFeatures(SeurObj, nfeatures = n)
  return(VariableFeatures(SeurObj))
}

get_ranked_genes <- function(expression) {
  #' Get a data.frame, where: rows are cells | ranks are cols | genes are in cells.
  #' The nth first elements of a row i correspond to the names of the nth most expressed genes of a cell i.
  #' 
  #' @param expression: a scRNA-seq dataset: rows are genes | cols are cells.
  #' 
  #' @return a data.frame, where: rows are cells | ranks are cols | genes are in cells.
  #' 
  get_ranks.cell <- function(expr.cell){rank(-expr.cell, ties.method="random")} # most expressed gene is 1
  get_genes.cell <- function(named_ranks.cell){names(sort(named_ranks.cell))}
  get_ranked_genes.cell <- function(expr.cell){get_genes.cell(get_ranks.cell(expr.cell))}
  ranked_genes <- apply(X=expression, MARGIN=2, FUN=get_ranked_genes.cell)
  return(ranked_genes)
}

filter_expressed_genes <- function(ranked_genes, expression) {
  #' Set NA values to genes not expressed in the ranked_genes data.frame.
  #' 
  #' @param ranked_genes: a data.frame, where: rows are cells | ranks are cols | genes are in cells.
  #' @param expression: a scRNA-seq dataset: rows are genes | cols are cells.
  #' 
  #' @return a data.frame, where: rows are cells | ranks are cols | genes or NA are in cells.
  #' 
  n_genes <- nrow(ranked_genes)
  count_expressed_genes.cell <- function(expr.cell){sum(expr.cell>0)}
  n_expressed_genes <- apply(X=expression, MARGIN=2, FUN=count_expressed_genes.cell)
  for (cell in colnames(ranked_genes)) {
    missing_genes <- (n_expressed_genes[cell]+1):n_genes
    ranked_genes[missing_genes, cell] <- NA
  }
  return(ranked_genes)
}

get_expression.count <- function(expression.init, n_HVGs) {
  #' Get a raw count expression matrix with HVGs only.
  #' This pre-processing step is required for SHARP, CIDR and scLCA.
  #' It is also required to get a log2 tpm matrix for scCCESS.Kmeans, scCCESS.SIMLR, densityCut.
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
  expression.count <- expression.count[, colSums(expression.count)>0]
    # remove cells with 0 gene expression on HVGs
  return(expression.count)
}

get_SeurObj.count <- function(expression.count) {
  #' Get a Seurat Object from a raw count expression matrix with HVGs only.
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

trim_data <- function(expression.init, population, records, params,
                      figures, random_state, SeurObj.init) {
  #' Trim a scRNA-seq dataset to keep only the cells of interest and a subset of their HVGs.
  #' 
  #' @param expression.init: the base scRNA-seq dataset: rows are genes | cols are cells.
  #' @param population: a character.
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param params: a list of parameters, with 'n_HVGs'.
  #' @param figures: a boolean. If TRUE, draw figures summarizing the trimming step. 
  #' If FALSE, SeurObj.init is not required.
  #' @param random_state: a numeric.
  #' @param SeurObj.init: a SeuratObject generated from the base scRNA-seq dataset. 
  #' A U-MAP has already been applied on the object.
  #' 
  #' @return a list of three data.frames: 'expression.loop' and 'SeurObj.loop', and 'ranked_genes.loop'.
  #' 
  data.loop <- list()
  cells_of_interest <- get_cells_of_interest(population, records$cells)
  
  if (length(cells_of_interest)>=100) { 
    # sncells=100 by default (SHARP).

    # trim the matrix to account for the cells of interest only and their HVGs
    #####################################################
    expression.loop <- expression.init[, cells_of_interest]
    expression.loop <- get_expression.count(expression.loop, params$n_HVGs)
    
    # set-up a SeuratObject for downstream methods and visualizations
    #####################################################
    SeurObj.loop <- get_SeurObj.count(expression.loop)
    SeurObj.loop <- RunUMAP(SeurObj.loop, features=VariableFeatures(SeurObj.loop),
                            seed.use=random_state)

    # rank and filter the expressed genes
    #####################################
    ranked_genes.loop <- get_ranked_genes(expression.loop)
    ranked_genes.loop <- filter_expressed_genes(ranked_genes.loop, expression.loop)
    
    # store the data for the iteration
    ##################################
    data.loop <- list(expression.loop=expression.loop,
                      SeurObj.loop=SeurObj.loop,
                      ranked_genes.loop=ranked_genes.loop)
  }
  
  if (figures) {
    trim_plot <- do_DimPlot(SeurObj.init, cells.highlight=cells_of_interest) + 
      ggtitle(population)
    trim_plot <- trim_plot +
      theme_bw() +
      theme(panel.background=element_rect(fill="lightgrey"),
            axis.title=element_blank(), legend.position="bottom")
    
    pdf(file = glue("{params$figures_dir}/{population}_trim.pdf"))
    print(trim_plot)
    dev.off()
  }
  return(data.loop)
}
