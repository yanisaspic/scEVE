"Functions called to generate results where the leftover cells are soft-clustered, w.r.t. their gene expression.

	2024/02/20 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})
source("./src/scEVE/trim.R")
source("./src/scEVE/genes.R")
source("./src/scEVE/leftovers_strategy/default.R")

update_leftover_seed <- function(seeds, population, sheet.cells, data.loop,
                                 occurrences.population, params) {
  #' Update the cells and the markers of the leftover seed to account for cells displaced after the soft-clustering.
  #' 
  #' @param seeds: a nested list, where each sub-list has five keys: 'clusters', 'consensus', 'cells', 'genes' and 'markers'.
  #' @param population: a character.
  #' @param data.loop: a list of three data.frames: 'expression.loop' and 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' @param occurrences.population: a data.frame where: genes are rows | sampling effort is cols | cells are occurrences.
  #' @param params: a list of parameters, with 'min_likelihood'.
  #' 
  #' @return a nested list, where each sub-list has five keys: 'clusters', 'consensus', 'cells', 'genes' and 'markers'.
  #'
  leftover_subpopulation <- glue("{population}{length(seeds)}")
  leftover_seed <- seeds[[length(seeds)]]
  leftover_seed$cells <- get_cells_of_interest(leftover_subpopulation, sheet.cells, params)
  leftover_seed$occurrences <- get_occurrences(data.loop$ranked_genes.loop[, leftover_seed$cells])
  
  # the markers kept must be over-represented prior, and after updating the seed.
  markers_after_update <- get_markers.seed(leftover_seed, occurrences.population)
  leftover_seed$markers <- intersect(leftover_seed$markers, markers_after_update) 
  leftover_seed$specific_markers <- intersect(leftover_seed$markers, leftover_seed$specific_markers)
  
  seeds[[length(seeds)]] <- leftover_seed
  return(seeds)
}

get_likelihoods.cell.naive <- function(ranked_genes.cell, seeds) {
  #' Get the membership likelihoods of a leftover cell according to the genes it expresses.
  #' It corresponds to the percentage of specific markers expressed by the cell.
  #' The likelihoods are normalized so that their sum is equal to 1.
  #'
  #' @param ranked_genes.cell: a vector of genes expressed.
  #' @param seeds: a nested list, where each sub-list has five keys: 'consensus', 'cells', 'clusters', 'markers' and 'specific_markers'.
  #' 
  #' @return a vector of numeric.
  #' 
  get_likelihood <- function(seed) {
    specific_markers_expressed <- intersect(ranked_genes.cell, seed$specific_markers)
    likelihood <- length(specific_markers_expressed) / length(seed$specific_markers)
    return(likelihood)
  }
  likelihoods.cell <- sapply(X=seeds, FUN=get_likelihood)
  likelihoods.cell <- likelihoods.cell / sum(likelihoods.cell)
  return(likelihoods.cell)
}

get_leftover_likelihoods <- function(ranked_genes, fun_likelihoods.cell, seeds, population) {
  #' Get a sheet of cell membership likelihood w.r.t. populations specific to cells.
  #' The value i,j in the results indicates the likelihood of a cell i belonging to the population j.
  #' The leftover cells are soft-clustered to the existing seeds according to the specific markers they express.
  #' The likelihood is a value between 0 and 1 related to the percentage of specific markers expressed by a leftover cell.
  #'
  #' @param ranked_genes: a data.frame where: ranks are rows | cells are cols | cells are genes.
  #' @param fun_likelihoods.cell: a function namespace. It is specific to the strategy used to compute the likelihood ("naive" or "weighted).
  #' @param seeds: a nested list, where each sub-list has five keys: 'consensus', 'cells', 'clusters', 'markers' and 'specific_markers'.
  #' @param population: a character.
  #'
  #' @return a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #'
  leftover_likelihoods <- apply(X=ranked_genes, MARGIN=2, FUN=fun_likelihoods.cell, seeds=seeds)
  name_subpopulation <- function(i){glue("{population}{i}")}
  rownames(leftover_likelihoods) <- sapply(X=1:nrow(leftover_likelihoods),
                                           FUN=name_subpopulation)
  leftover_likelihoods <- t(leftover_likelihoods)
  return(leftover_likelihoods)
}

get_sheet.cells.soft <- function(records, seeds, population, data.loop, params) {
  #' Get a sheet of cell membership likelihood w.r.t. populations.
  #' The value i,j in the results indicates the likelihood of a cell i belonging to the population j.
  #' The leftover cells are soft-clustered to the existing seeds according to the specific markers they express.
  #' The likelihood is a value between 0 and 1 related to the percentage of specific markers expressed by a leftover cell.
  #' The markers of the leftover cells are updated as cells are sequentially displaced.
  #'
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has five keys: 'consensus', 'cells', 'clusters', 'markers' and 'specific_markers'.
  #' @param population: a character.
  #' @param data.loop: a list of three data.frames: 'expression.loop' and 'SeurObj.loop', and 'ranked_genes.loop'.
  #' 
  #' @return a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  sheet.cells <- get_sheet.cells.default(records, seeds, population)
  occurrences.population <- get_occurrences(data.loop$ranked_genes.loop)
  leftover_cells.init <- seeds[[length(seeds)]]$cells
  leftover_cells.previous <- c()
  leftover_cells.current <- leftover_cells.init
  
  if (params$leftovers_strategy=="naive") {fun_likelihoods.cell <- get_likelihoods.cell.naive}
  if (params$leftovers_strategy=="weighted") {fun_likelihoods.cell <- get_likelihoods.cell.weighted}
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # prune the leftover seed to generate a homogeneous population w.r.t. markers
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  while (!identical(leftover_cells.previous, leftover_cells.current)) {
    ranked_genes.current <- data.loop$ranked_genes.loop[, leftover_cells.current]
    leftover_likelihoods.current <- get_leftover_likelihoods(ranked_genes.current,
                                                             fun_likelihoods.cell, 
                                                             seeds, population)
    
    seeds <- update_leftover_seed(seeds, population, leftover_likelihoods.current,
                                  data.loop, occurrences.population, params)
    # the leftover seed is updated (cells, markers and specific_markers)
    # to account for cells displaced after likelihood calculation.
    
    leftover_cells.previous <- leftover_cells.current
    leftover_cells.current <- seeds[[length(seeds)]]$cells
  } #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # the likelihoods are calculated once with the stable updated marker genes
  ranked_genes.init <- data.loop$ranked_genes.loop[, leftover_cells.init]
  leftover_likelihoods.init <- get_leftover_likelihoods(ranked_genes.init,
                                                        fun_likelihoods.cell,
                                                        seeds, population)
  sheet.cells[leftover_cells.init,] <- leftover_likelihoods.init
  
  return(sheet.cells)
}

plot_clusters <- function(sheet.cells, SeurObj.init){}