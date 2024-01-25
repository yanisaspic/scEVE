"Misc. functions called by multiple R scripts.

	2024/01/23 @yanisaspic"

suppressPackageStartupMessages({
  library(qpdf)
  library(dplyr)
  library(readxl)
  library(openxlsx)
  
  library(glue)
  library(scales)
  library(igraph)
  library(Seurat)
  library(SCpubr)
  library(ggplot2)
  library(cowplot)
})

init_records <- function(expression.init, params) {
  #' Get a named list of three data.frames, with: 
  #' - 'cells': rows are cells | cols are populations. 
  #' The value of the cell i,j corresponds to the probability that a cell i belongs to a population j.
  #' - 'markers': rows are genes | cols are populations. 
  #' The value of the cell i,j is 1 if a cell is a marker of a population. Else, it is 0.
  #' - 'meta': rows are populations | cols are 'consensus', 'parent', 'n' and 'to_do'.
  #' 
  #' @param expression.init: the base scRNA-seq dataset: rows are genes | cols are cells.
  #' @param params: a list of parameters, with 'root_consensus'.
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
    parent=as.character(NA),
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

get_existing_pdfs <- function(population) {
  #' Get the names of the existing pdf files (intermediate figures) w.r.t a population.
  #'
  #' @param population: a character.
  #'
  #' @return a vector of filenames.
  #'
  files <- c()
  for (cat in c("trim", "clusterings", "seeds", "genes")) {
    filename <- glue("./figures/{population}_{cat}.pdf")
    if (file.exists(filename)) {
      files <- c(files, filename)
    }
  }
  return(files)
}

merge_pdfs <- function(population) {
  #' Merge a group of pdf files together.
  #' 
  #' @param population: a character.
  #' 
  files <- get_existing_pdfs(population)
  pdf_combine(input = files, output = glue("./figures/{population}.pdf"))
  unlink(files)
}

update_records <- function(records, seeds, population) {
  #' Add the results of a loop to the records.
  #'
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has five keys: 'clusters', 'consensus', 'cells', 'genes' and 'markers'.
  #' @param population: a character.
  #'
  #' @return a named list of three data.frames: 'cells', 'meta' and 'markers'.
  #'
  for (i in 1:length(seeds)) {
    s <- seeds[[i]]
    cluster_label <- glue("{population}{i}")
    
    # meta
    ######
    row.meta <- c(round(s$consensus,2), population, length(s$cells), TRUE)
    records$meta[cluster_label,] <- row.meta
    records$meta$consensus <- as.numeric(records$meta$consensus)
    
    # cells & markers
    #################
    col.cells <- rownames(records$cells) %in% s$cells
    records$cells[, cluster_label] <- as.numeric(col.cells)
    col.markers <- rownames(records$markers) %in% s$markers
    records$markers[, cluster_label] <- as.numeric(col.markers)
  }
  return(records)
}

trim_data <- function(expression.init, 
                      population, 
                      records, 
                      params,
                      figures,
                      random_state,
                      SeurObj.init) {
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
  cells_of_interest <- rownames(records$cells)[records$cells[, population]==1]
  
  if (length(cells_of_interest)>=100) {
    
    # get discriminant transcriptome of the subpopulation
    #####################################################
    expression.loop <- expression.init[, cells_of_interest]
    HVGs.loop <- get_n_HVGs(expression.loop, n=params$n_HVGs)
    expression.loop <- expression.loop[HVGs.loop,]
    
    # set-up data for SeurObj-dependent clustering algorithms
    #########################################################
    SeurObj.loop <- CreateSeuratObject(expression.loop)
    VariableFeatures(SeurObj.loop) <- HVGs.loop
    SeurObj.loop <- NormalizeData(SeurObj.loop)
    SeurObj.loop <- ScaleData(SeurObj.loop,
                              features=VariableFeatures(SeurObj.loop))
    SeurObj.loop <- RunUMAP(SeurObj.loop,
                            features=VariableFeatures(SeurObj.loop),
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
    trim_plot <- do_DimPlot(SeurObj.init, cells.highlight=cells_of_interest)
    pdf(file = glue("./figures/{population}_trim.pdf"))
    print(trim_plot)
    dev.off()
  }
    
  return(data.loop)
}

get_classification <- function(records, ground_truth=FALSE) {
  #' Get a data.frame associating every unique cell to its most informative cluster label.
  #' 
  #' @param records: a list of three data.frames: 'meta', 'cells' and 'markers'.
  #' @param ground_truth: a boolean. If TRUE, cell names must follow the {ground_truth}_{i} format.
  #' Add a column corresponding to the ground truth labels.
  #' 
  #' @return a data.frame with two columns: 'cell' and 'pred'. 
  #' If ground_truth is TRUE, there is a third column 'ground'.
  #' 
  all_cells <- rownames(records$cells)
  classification <- data.frame(cell=all_cells, pred="C")
  for (cluster in colnames(records$cells)) {
    is_in_cluster <- records$cells[, cluster]==1
    cells_in_cluster <- all_cells[is_in_cluster]
    classification[classification$cell %in% cells_in_cluster, "pred"] <- cluster
  }
  
  if (ground_truth) {
    ground_truth <- sapply(strsplit(classification$cell, split="_"), "[", 1)
    classification$ground <- ground_truth
  }
  return(classification)
}

get_piecharts.population <- function(population, classification, records) {
  #' Get the number of individuals from each ground truth group within a population.
  #'
  #' @param population: a character.
  #' @param classification: a data.frame with three columns: 'cell', 'pred' and 'ground'. 
  #' @param records: a list of three data.frames: 'meta', 'cells' and 'markers'.
  #'
  #' @return a vector of numeric.
  #' 
  cells_of_interest <- rownames(records$cells[records$cells[population]==1,])
  data <- classification[classification$cell %in% cells_of_interest,]
  
  data <- data %>% count(ground)
  piechart <- setNames(object = data$n, data$ground)
  ground_groups <- unique(classification$ground)
  missing_groups <- setdiff(ground_groups, data$ground)
  
  piechart[missing_groups] <- 0
  piechart <- piechart[ground_groups]
  return(piechart)
}

get_piecharts <- function(records, classification) {
  #' For each population, get the number of individuals composing it, split by ground truth groups.
  #' 
  #' @param records: a list of three data.frames: 'meta', 'cells' and 'markers'.
  #' @param classification: a data.frame with three columns: 'cell', 'pred' and 'ground'. 
  #' 
  #' @return a list of vectors.
  #' The vector elements are ordered and correspond to the number of individuals in each ground truth group.
  #'
  pred_groups <- colnames(records$cells)
  piecharts <- lapply(X=pred_groups,
                      FUN=get_piecharts.population,
                      records=records,
                      classification=classification)
  names(piecharts) <- (pred_groups)
  return(piecharts)
}

draw_scEVE <- function(records) {
  #' Draw a graph with pie charts as vertices.
  #' The vertices are the cell population and the edges their hierarchy.
  #' 
  #' @param records: a list of three data.frames: 'meta', 'cells' and 'markers'.
  #' 
  #' @return a graph.
  #'
  data <- data.frame(parent=records$meta$parent,
                     child=row.names(records$meta), 
                     consensus=records$meta$consensus,
                     size=records$meta$n)
  hierarchy_graph <- graph.data.frame(data[-1,], directed=TRUE)
  E(hierarchy_graph)$label <- data[-1,]$consensus
  nodes_order <- vertex_attr(hierarchy_graph)$name
  
  # set nodes pie charts w.r.t. ground truth
  ##########################################
  classification <- get_classification(records, ground_truth = TRUE)
  ground_groups <- unique(classification$ground)
  colors <- hue_pal()(length(ground_groups))
  piecharts <- get_piecharts(records, classification)
  piecharts <- piecharts[nodes_order]

  g <- plot(hierarchy_graph,
       layout=layout_as_tree,
       root=-1,
       vertex.shape="pie",
       vertex.pie=piecharts,
       vertex.pie.color=list(colors))
  
  return(g)
}
