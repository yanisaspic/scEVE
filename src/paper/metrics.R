"Functions used to compute metrics for the benchmark.

	2024/04/02 @yanisaspic"

suppressPackageStartupMessages({
  library(aricode)  # packageVersion("aricode")==1.0.3
  library(ggplot2)
  library(RColorBrewer)
  library(reshape2)
})
source("./src/scEVE/utils/misc.R")

get_ground_truth <- function(expression.init) {
  #' Get a vector of labels corresponding to the ground truth from an expression matrix.
  #' The columns correspond to a unique cell id, and a cell id has the following structure:
  #' {ground_truth}_{i}, where i is a numeric.
  #' 
  #' @param expression.init: a scRNA-seq count matrix:
  #' genes are rows | cells are cols.
  #' 
  #' @return a vector of ground truths.
  #' 
  get_ground_truth.cell_id <- function(cell_id) {
    tmp <- strsplit(cell_id, split="_")[[1]]
    ground_truth.cell_id <- paste(head(tmp, -1), collapse="_")
    return(ground_truth.cell_id)
  }
  
  cell_ids <- colnames(expression.init)
  ground_truth <- sapply(X=cell_ids, FUN=get_ground_truth.cell_id)

  ground_truth <- ground_truth[order(names(ground_truth))]
  # sort cells alphabetically to facilitate benchmarking.
  ground_truth <- factor(ground_truth)
  return(ground_truth)
}

get_metric <- function(preds, ground_truth, metric) {
  #' Compute a classification metric from a vector of preds and a ground truth.
  #'
  #' @param preds: a named vector: names are cell ids and values are cluster labels predicted.
  #' @param ground_truth: a named vector: names are cell ids and values are cluster labels of the authors.
  #' @param metric: a valid metric name. One of: 'ARI', 'AMI' and 'NMI'.
  #' 
  #' @return a numeric.
  #' 
  if (!identical(names(preds), names(ground_truth))) {
    stop("The cell ids of the predictions and the ground truth are not aligned.")
  }
  result <- get(metric)(preds, ground_truth)
  return(result)
}

get_distribution.population <- function(population, sheet.cells, ground_truth) {
  #' Get the distribution of ground truth labels for a given population.
  #' 
  #' @param population: a character.
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' @param ground_truth: a named vector: names are cell ids and values are cluster labels of the authors.
  #' 
  #' @return a vector of numerics.
  #' 
  cells_of_interest <- get_cells_of_interest(population, sheet.cells)
  data <- table(cells_of_interest, ground_truth[cells_of_interest])
  distribution <- apply(X=data, MARGIN=2, FUN=sum)
  return(distribution)
}

draw_tree <- function(records, ground_truth) {
  #' Draw a hierarchical graph with pie charts as vertices.
  #' The vertices are the cell population and the edges their hierarchy.
  #' The consensus of a population is written over the edge to its parent.
  #' The ground truth labels are represented by the colors in the pie charts.
  #' 
  #' @param records: a named list of three data.frames: 'meta', 'cells' and 'markers'.
  #' @param ground_truth: a named vector: names are cell ids and values are cluster labels of the authors.
  #' 
  #' @return a graph.
  #'

  #_________________________________________________________________________data
  data <- data.frame(parent=records$meta$parent,
                     child=row.names(records$meta), 
                     consensus=round(records$meta$consensus,2))
  tree <- graph.data.frame(data[-1,], directed=TRUE)
  ordered_nodes <- vertex_attr(tree)$name
  
  #_________________________________________________________________________draw
  E(tree)$label <- data[-1,]$consensus
  E(tree)$label.color <- "black"
  E(tree)$label.cex <- 1.5
  
  V(tree)$label.color <- "black"
  V(tree)$label.cex <- 1.5
  V(tree)$label.degree <- 2*pi/3
  V(tree)$label.dist <- 0.6
  
  colors <- brewer.pal(n=length(levels(ground_truth)), name="Dark2")
    # Dark2 is color-blind friendly
  distributions <- lapply(X=ordered_nodes, FUN=get_distribution.population,
                          sheet.cells=records$cells, ground_truth=ground_truth)
  plot(tree, layout=layout_as_tree, root=-1, vertex.shape="pie", vertex.size=40,
       vertex.pie=distributions, vertex.pie.color=list(colors))
  # this figure is manually screen captured and annotated for the manuscript.
}

draw_heatmap <- function(preds, ground_truth) {
  #' Draw a heatmap where rows are cluster predictions and columns are clusters assigned by authors.
  #' 
  #' @param preds: a named vector: names are cell ids and values are cluster labels predicted.
  #' @param ground_truth: a named vector: names are cell ids and values are cluster labels of the authors.
  #' 

  #_________________________________________________________________________data
  data <- table(preds, ground_truth)
  maximum_n <- apply(X=data, MARGIN=2, FUN=sum)
  data.ggplot2 <- melt(data)
  colnames(data.ggplot2) <- c("output_clusters", "expected_clusters", "n")
  get_proportion.row <- function(row) {as.numeric(row["n"]) / maximum_n[row["expected_clusters"]]}
  data.ggplot2$proportion <- apply(X=data.ggplot2, MARGIN=1, FUN=get_proportion.row)
  
  #_________________________________________________________________________draw
  p <- ggplot(data.ggplot2, aes(x=output_clusters, y=expected_clusters)) +
    geom_tile(aes(fill=proportion), color="white", lwd=1.5, linetype=1) +
    geom_text(aes(label=n), color="black", size=6) +
    coord_fixed() +
    scale_fill_gradient(low="#FFFFFF", high="#666666") +
    guides(fill="none") +
    theme_bw() +
    theme(axis.text = element_text(size=12)) +
    theme(axis.title = element_text(size=15))
  return(p)
}
