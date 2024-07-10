"Functions used to identify seeds, i.e. groups of cells unambiguously assigned together.
In the papers, we refer to seeds as consensus clusters. Rules refer to overlaps in the paper.

	2024/05/24 @yanisaspic"

suppressPackageStartupMessages({
  library(arules)
  library(dplyr)
  library(glue)
  library(igraph)
  library(rlang)
  library(SCpubr)
})

get_transactions <- function(clusterings) {
  #' Get transactions from a table of clusterings.
  #' 
  #' @param clusterings: a data.frame where: rows are cells | cols are clusterings | cells are labels.
  #' 
  #' @return an object of the class 'transactions'.
  #'
  key <- hash(clusterings)
  # key is generated to prevent overlap with parallel runs.
  clusterings_path <- glue("./null/{key}.tmp")
  write.table(clusterings, file=clusterings_path, col.names=FALSE, row.names=FALSE)
  transactions <- read.transactions(clusterings_path)
  file.remove(clusterings_path)
  return(transactions)
}

get_rules <- function(transactions, min_support, min_confidence) {
  #' Get the rules associating two clusters above a confidence and a support threshold.
  #'
  #' @param transactions: a formal class transactions object.
  #' @param min_support: a numeric between 0 and 1
  #' @param min_confidence: a numeric between 0 and 1.
  #' 
  #' @return a data.frame with three columns: 'A', 'B' and 'confidence'.
  #' 
  data <- apriori(transactions, 
                  parameter = list(support=min_support,
                                   confidence=min_confidence,
                                   minlen=2,
                                   maxlen=2,
                                   target="rules")
  )
  out <- DATAFRAME(data)
  
  get_item <- function(x){substr(x, 2, nchar(x)-1)}
  rules <- list(
    A=sapply(X=as.character(out$LHS), FUN=get_item),
    B=sapply(X=as.character(out$RHS), FUN=get_item),
    confidence=as.numeric(out$confidence)
  )
  rules <- do.call(cbind, rules)
  rules <- data.frame(rules, row.names = NULL)
  return(rules)
}

get_symmetric_rules <- function(rules) {
  #' Filter and merge the symmetric rules into a single one. 
  #' The confidence of the new rule is the minimal confidence of its two symmetric rules.
  #' 
  #' @param rules: a data.frame with three columns: 'A', 'B' and 'confidence'.
  #' 
  #' @return a data.frame with three columns: 'A', 'B' and 'confidence'.
  #' 
  rules$confidence <- as.numeric(rules$confidence)
  rules <- rules[order(rules$confidence, decreasing = TRUE),]
  # decreasing sort: duplicate keys are symmetric rules with minimum confidence
  
  get_key.row <- function(row) {
    elements <- sort(row[c("A", "B")])
    key <- paste(elements, collapse=".")}
  keys <- apply(X=rules, MARGIN=1, FUN=get_key.row)
  
  rules <- rules[duplicated(keys),]
  return(rules)
}

get_max_consensus <- function(methods) {
  #' Get the value of the maximum consensus expected w.r.t. the number of methods n.
  #' It corresponds to the number of edges in a complete non-directed graph with n nodes.
  #' 
  #' @param methods: a vector of clustering methods used.
  #' 
  #' @return a numeric
  #'
  n_methods <- length(methods)
  n_edges <- n_methods * (n_methods - 1) / 2
  return(n_edges)
}

get_hint.component <- function(connected_component,
                               max_consensus) {
  #' Get a cluster hint from a connected component.
  #' It is a group of clusters that overlap across multiple methods.
  #'
  #' @param connected_component: a connected_component with 'name' and 'confidence' attributes.
  #' @param max_consensus: a numeric.
  #' 
  #' @return a list with two keys: 'consensus' and 'clusters'.
  #' 
  clusters <- names(V(connected_component))
  edges_confidence <- as.numeric(E(connected_component)$confidence)
  consensus <- sum(edges_confidence)
  hint <- list(clusters=clusters, consensus=consensus/max_consensus)
  return(hint)
}

get_hints <- function(rules, max_consensus) {
  #' Use the connected components of a graph to identify raw seeds.
  #'
  #' @param rules: a data.frame with three columns: 'A', 'B' and 'confidence'.
  #' @param max_consensus: a numeric.
  #' 
  #' @return a nested list, where each sub-list has two keys: 'consensus' and 'clusters'.
  #' 
  overlap_graph <- graph.data.frame(rules)
  connected_components <- decompose.graph(overlap_graph)
  hints <- lapply(X=connected_components, 
                  FUN=get_hint.component,
                  max_consensus=max_consensus)
  return(hints)
}

get_method.cluster <- function(cluster){
  #' Get the name of the method used to define an independent cluster.
  #' 
  #' @param cluster: a character in the format {method}_{i}.
  #' 
  #' @return a character.
  #' 
  method <- strsplit(cluster, split="_")[[1]][1]
  return(method)
}

has_independent_clusters_only <- function(hint) {
  #' Check if a raw seed has independent clusters only, i.e. no more than one cluster per method.
  #' 
  #' @param hint: a list with two keys: 'consensus' and 'members'.
  #' 
  #' @return a logical
  #'
  clusters <- hint$clusters
  methods <- as.character(sapply(X=clusters, FUN=get_method.cluster))
  result <- length(unique(methods)) == length(clusters) 
  return(result)
}

filter_conflictual_hints <- function(hints, expected_consensus) {
  #' Remove conflictual hints, i.e:
  #' - hints with two dependent members, or 
  #' - hints below or equal to a consensus threshold.
  #' 
  #' @param hints: a nested list, where each sub-list has two keys: 'consensus' and 'clusters'.
  #' @param expected_consensus: a numeric.
  #' 
  #' @return a nested list, where each sub-list has two keys: 'consensus' and 'clusters'
  #' 
  is_more_consensual <- function(hint){hint$consensus > expected_consensus}
  is_consensual_hint <- function(hint){is_more_consensual(hint) & has_independent_clusters_only(hint)}
  consensual_hints <- Filter(f=is_consensual_hint, x=hints)
  return(consensual_hints)
}

get_cells.hint <- function(hint, clusterings) {
  #' Get the names of the cells at the intersection of all the clusters of a hint.
  #' 
  #' @param hint: a list with two keys: 'consensus' and 'clusters'.
  #' @param clusterings: a data.frame: cells are rows | methods are cols | a cell is a label.
  #' 
  #' @return a vector of cell names.
  #' 
  data <- clusterings
  methods <- as.character(sapply(X=hint$clusters, FUN=get_method.cluster))
  data <- data[, methods]
  
  is_in_hint.cluster <- function(cluster){as.numeric(cluster %in% hint$clusters)}
  data <- apply(X=data, MARGIN = c(1,2), FUN=is_in_hint.cluster)
  data <- data[rowSums(data)==ncol(data),]
  
  cells <- rownames(data)
  return(cells)
}

get_cells <- function(hints, clusterings) {
  #' Get a nested list of cells with two keys: 'hint' and 'all'.
  #' The sub-list 'hint' is the cells w.r.t. hint.
  #' The vector 'all' is the aggregation of the cells, with duplicates.
  #' 
  #' @param hints: a nested list, where each sub-list has two keys: 'consensus' and 'clusters'.
  #' @param clusterings: a data.frame: cells are rows | methods are cols | a cell is a label.
  #' 
  #' @return a nested list with two keys: 'hint' and 'all'.
  #' 
  cells_by_hint <- lapply(X=hints,
                          FUN=get_cells.hint,
                          clusterings=clusterings)
  cells <- list(hints=cells_by_hint, all=unlist(cells_by_hint))
  return(cells)
}

add_cells_to_hints <- function(hints, cells) {
  #' Add the cells specific to a hint. A hint with a list of specific cells is called a 'seed'.
  #'
  #' @param hints: a nested list, where each sub-list has two keys: 'consensus' and 'clusters'.
  #' @param cells: a nested list with two keys: 'meta' and 'all'.
  #'
  #' @return a nested list, where each sub-list has three keys: 'consensus', 'cells' and 'clusters'.
  #'
  all_cells <- cells[["all"]]
  unspecific_cells <- unique(all_cells[duplicated(all_cells)])
  for (i in 1:length(hints)) {
    cells_of_hint <- cells$hints[[i]]
    specific_cells_of_hint <- setdiff(cells_of_hint, unspecific_cells)
    hints[[i]]$cells <- specific_cells_of_hint
  }
  return(hints)
}

add_seeds <- function(SeurObj, seeds, population) {
  #' Add a seed factor to a Seurat Object.
  #' 
  #' @param SeurObj: a Seurat Object of raw count expression, with selected genes.
  #' FindVariableFeatures() and ScaleData() must have been applied already.
  #' @param seeds: a nested list, where each sub-list has three keys: 'consensus', 'cells' and 'clusters'.
  #' @param population: a character.
  #' 
  #' @return: a Seurat Object where seeds are accessible with SeurObj$seeds.
  #'
  get_cells.seed <- function(i) {
    cells_of_seed <- seeds[[i]]$cells
    label <- glue("{population}.{i}")
    if (i == length(seeds)) {label <- paste(label, "(leftovers)")}
    out <- setNames(rep(label, length(cells_of_seed)), cells_of_seed)
    return(out)
  }
  cells <- lapply(X=1:length(seeds),
                  FUN=get_cells.seed)
  cells <- do.call(c, cells)
  cells <- as.factor(cells) 
  SeurObj$seeds <- cells
  return(SeurObj)
}

get_minimal_support <- function(expression.init, clusterings, params) {
  #' Get the minimal support expected in a hint w.r.t. the global minimal support.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param clusterings: a data.frame where: rows are cells | cols are clusterings | cells are labels.
  #' @param params: a list of parameters.
  #' 
  #' @return a numeric.
  #' 
  n_cells.init <- ncol(expression.init)
  expected_min_n_cells <- params$min_prop_cells * n_cells.init
  n_cells.loop <- nrow(clusterings)
  minimal_support <- expected_min_n_cells / n_cells.loop
  return(minimal_support)
}

draw_seeds <- function(data.loop, seeds, population, params) {
  #' Draw the seeds identified with a U-MAP.
  #'
  #' @param data.loop: a list of three data.frames: 'expression.loop' and 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param seeds: a nested list, where each sub-list has three keys: 'consensus', 'cells' and 'clusters'.
  #' @param population: a character.
  #' @param params: a named list with 'figures_dir'.
  #' 
  data.loop$SeurObj.loop <- add_seeds(data.loop$SeurObj.loop, seeds, population)
  pdf(file = glue("{params$figures_dir}/{population}_seeds.pdf"))
  seeds_plot <- do_DimPlot(data.loop$SeurObj.loop,
                           split.by="seeds",
                           legend.position="none")
  
  for (i in 1:length(seeds_plot)) {
    seeds_plot[[i]][[1]] <- seeds_plot[[i]][[1]] +
      theme_bw() +
      theme(plot.title=element_text(hjust=0.5, margin=margin(1, 0, 0, 0)),
            panel.background=element_rect(fill="lightgrey"),
            legend.position="none",
            axis.title=element_blank())
  }
  
  print(seeds_plot)
  dev.off()
}

get_seeds <- function(expression.init, data.loop, clusterings, params, records, population, figures) {
  #' Get consensual seeds from a table of cluster labels.
  #'
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param data.loop: a list of three data.frames: 'expression.loop' and 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param clusterings: a data.frame where: rows are cells | cols are clusterings | cells are labels.
  #' @param params: a list of parameters.
  #' @param records: a list of three data.frames: 'meta', 'cells', 'markers'.
  #' @param population: a character.
  #' @param figures: a boolean. If TRUE, draw figures summarizing the seed digging.
  #'
  #' @return a nested list, where each sub-list has three keys: 'consensus', 'cells' and 'clusters'.
  #' 
  
  # hyperparameters
  #################
  minimal_support <- get_minimal_support(expression.init, clusterings, params)
  threshold_relative_consensus <- max(records$meta[population, "consensus"], params$root_consensus)
  max_absolute_consensus <- get_max_consensus(methods=colnames(clusterings))
  
  # association mining
  #####################
  transactions <- get_transactions(clusterings)
  rules <- get_rules(transactions, minimal_support, min_confidence = .5)
  rules <- get_symmetric_rules(rules)
  hints <- get_hints(rules, max_absolute_consensus)
  
  # seeds
  #######
  seeds <- list()
  while (TRUE) {
    if(length(hints)==0){break()}
    
    # filter conflictual hints
    ##########################
    hints <- filter_conflictual_hints(hints, threshold_relative_consensus)
    if(length(hints)==0){break()}
    
    # identify unambiguous cells for each hint
    ##########################################
    cells.hints <- get_cells(hints, clusterings)
    seeds <- add_cells_to_hints(hints, cells.hints)
    has_cells <- function(seed){length(seed$cells)>1}
    seeds <- Filter(f=has_cells, x=seeds)
    
    if(length(seeds)==0){break()}
    consensi <- sapply(seeds, "[[", "consensus")
    seeds <- seeds[order(-consensi)]
    
    # generate a leftover seed
    ##########################
    cells.seeds <- unlist(sapply(seeds, "[[", "cells"))
    missing_cells <- setdiff(rownames(clusterings), cells.seeds)
    leftover_seed <- list(
      consensus=0,
      clusters=c(),
      cells=missing_cells
    )
    seeds[[length(seeds) + 1]] <- leftover_seed
    
    # draw the seeds
    ################
    if (figures) {draw_seeds(data.loop, seeds, population, params)}
    
    break() 
  }
  return(seeds)
}
