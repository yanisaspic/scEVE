"Functions called to generate results with the default strategy: leftover cells are in the leftover seed.

	2024/02/20 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})
source("./src/scEVE/trim.R")
source("./src/scEVE/genes.R")
source("./src/scEVE/leftovers_strategy/default.R")

get_likelihoods.cell <- function(ranked_genes.cell, seeds) {
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

get_sheet.cells.soft <- function(records, seeds, population, data.loop) {
  #' Get a sheet of cell membership likelihood w.r.t. populations.
  #' The value i,j in the results indicates the likelihood of a cell i belonging to the population j.
  #' The leftover cells are soft-clustered to the existing seeds according to the specific markers they express.
  #' The likelihood is a value between 0 and 1 related to the percentage of specific markers expressed by a leftover cell.
  #'
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has five keys: 'consensus', 'cells', 'clusters', 'markers' and 'specific_markers'.
  #' @param population: a character.
  #' @param data.loop: a list of three data.frames: 'expression.loop' and 'SeurObj.loop', and 'ranked_genes.loop'.
  #' 
  #' @return a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  sheet.cells <- get_sheet.cells.default(records, seeds, population)
  
  leftover_cells <- seeds[[length(seeds)]]$cells
  leftover_likelihoods <- apply(X=data.loop$ranked_genes.loop[, leftover_cells],
                                MARGIN=2,
                                FUN=get_likelihoods.cell,
                                seeds=seeds)
  
  name_subpopulation <- function(i){glue("{population}{i}")}
  rownames(leftover_likelihoods) <- sapply(X=1:nrow(leftover_likelihoods), FUN=name_subpopulation)
  leftover_likelihoods <- t(leftover_likelihoods)
  
  sheet.cells[leftover_cells,] <- leftover_likelihoods
  print(sheet.cells)
  return(sheet.cells)
}

get_sheet.meta.soft <- function(){}
update_markers <- function(data.loop, sheet.cells){}
plot_clusters <- function(sheet.cells, SeurObj.init){}