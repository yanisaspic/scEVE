"Functions called to generate results where the leftover cells are soft-clustered,
  w.r.t. their gene expression.

	2024/02/20 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(pracma)
})
source("./src/scEVE/genes.R")
source("./src/scEVE/utils/misc.R")
source("./src/scEVE/leftovers_strategy/default.R")

update_leftover_seed <- function(seeds, population, sheet.cells, data.loop) {
  #' Update the cells and the markers of the leftover seed to account for cells displaced after the soft-clustering.
  #' 
  #' @param seeds: a nested list, where each sub-list has five keys: 'clusters', 'consensus', 'cells', 'genes' and 'markers'.
  #' @param population: a character.
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' @param data.loop: a list of four data.frames: 'expression.loop', 'occurrences.loop', 'SeurObj.loop', and 'ranked_genes.loop'.
  #' 
  #' @return a nested list, where each sub-list has five keys: 'clusters', 'consensus', 'cells', 'genes' and 'markers'.
  #'
  leftover_subpopulation <- glue("{population}{length(seeds)}")
  leftover_seed <- seeds[[length(seeds)]]
  leftover_seed$cells <- get_cells_of_interest(leftover_subpopulation, sheet.cells)
  leftover_seed$occurrences <- get_occurrences(data.loop$ranked_genes.loop[, leftover_seed$cells])
  
  markers_after_update <- get_markers.seed(leftover_seed, data.loop$occurrences.loop)
  leftover_seed$markers <- intersect(leftover_seed$markers, markers_after_update) 
  leftover_seed$specific_markers <- intersect(leftover_seed$markers, leftover_seed$specific_markers)
  # the markers kept after the update are also over-represented before the update.
  
  seeds[[length(seeds)]] <- leftover_seed
  return(seeds)
}

get_likelihoods.cell.naive <- function(ranked_genes.cell, seeds) {
  #' Get the membership likelihoods of a leftover cell according to the genes it expresses.
  #' It corresponds to the percentage of specific markers expressed.
  #' The likelihoods are normalized so that their sum is equal to 1.
  #'
  #' @param ranked_genes.cell: a vector of genes expressed.
  #' @param seeds: a nested list, where each sub-list has five keys: 'consensus', 'cells', 'clusters', 'markers' and 'specific_markers'.
  #' 
  #' @return a vector of numeric.
  #' 
  get_likelihood.seed <- function(seed) {
    specific_markers_expressed <- intersect(ranked_genes.cell, seed$specific_markers)
    likelihood <- length(specific_markers_expressed) / length(seed$specific_markers)
    return(likelihood)
  }
  likelihoods.cell <- sapply(X=seeds, FUN=get_likelihood.seed)
  return(likelihoods.cell)
}

get_likelihoods.cell.weighted <- function(ranked_genes.cell, seeds) {
  #' Get the membership likelihoods of a leftover cell according to the genes it expresses.
  #' It corresponds to the AUC of the cumulative percentage of specific markers expressed w.r.t. the effort.
  #' The likelihoods are normalized so that their sum is equal to 1.
  #' 
  #' @param ranked_genes.cell: a vector of genes expressed.
  #' @param seeds: a nested list, where each sub-list has five keys: 'consensus', 'cells', 'clusters', 'markers' and 'specific_markers'.
  #' 
  #' @return a vector of numeric.
  #' 
  maximal_effort <- sum(!is.na(ranked_genes.cell))
  get_likelihoods.effort <- function(effort) {get_likelihoods.cell.naive(ranked_genes.cell[1:effort], seeds)}
  likelihoods.efforts <- sapply(X=1:maximal_effort, get_likelihoods.effort)
  # get the likelihoods with a variable effort (i.e. a variable number of top genes sampled)
  
  likelihoods.efforts <- cbind(0, likelihoods.efforts) # -> area calculation requires at least two points
  get_likelihoods.seed <- function(i) {trapz(x=0:maximal_effort, y=likelihoods.efforts[i,])}
  likelihoods.cell <- sapply(X=1:nrow(likelihoods.efforts), FUN=get_likelihoods.seed)
  # the average likelihoods regardless of effort is the area under curve of a plot
  # where the x-axis is the effort and the y-axis is the likelihood.
  
  return(likelihoods.cell)
}

get_leftover_likelihoods <- function(ranked_genes, fun_likelihoods.cell, seeds, population) {
  #' Get a sheet of cell membership likelihood w.r.t. populations specific to cells.
  #' The value i,j in the results indicates the likelihood of a cell i belonging to the population j.
  #' The leftover cells are soft-clustered to the existing seeds according to the specific markers they express.
  #' The likelihoods are normalized for each cell so that their sum is equal to 1.
  #' Thus, the likelihood is a value between 0 and 1 related to the percentage of specific markers expressed by a leftover cell.
  #'
  #' @param ranked_genes: a data.frame where: ranks are rows | cells are cols | cells are genes.
  #' @param fun_likelihoods.cell: a function namespace. It is specific to the strategy used to compute the likelihood ("naive" or "weighted).
  #' @param seeds: a nested list, where each sub-list has five keys: 'consensus', 'cells', 'clusters', 'markers' and 'specific_markers'.
  #' @param population: a character.
  #'
  #' @return a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #'
  leftover_likelihoods <- apply(X=ranked_genes, MARGIN=2, FUN=fun_likelihoods.cell, seeds=seeds)
  total_likelihoods <- apply(X=leftover_likelihoods, MARGIN=2, FUN=sum)
  get_normalized_likelihoods <- function(cell) {leftover_likelihoods[, cell] / total_likelihoods[cell]}
  normalized_likelihoods <- sapply(X=colnames(leftover_likelihoods), FUN=get_normalized_likelihoods)
  # normalize the likelihoods so that their sum is equal to 1 for each cell.
  
  replace_infinites <- function(x) {ifelse(is.finite(x), x, 0)}
  normalized_likelihoods <- apply(X=normalized_likelihoods, MARGIN=c(1,2), FUN=replace_infinites)
  # if no marker gene is expressed by a cell, the normalized leftover likelihoods is not finite.
  # these values are replaced with 0: no marker gene is expressed at all, so all the likelihoods are 0.
  
  name_subpopulation <- function(i){glue("{population}{i}")}
  rownames(normalized_likelihoods) <- sapply(X=1:nrow(normalized_likelihoods),
                                             FUN=name_subpopulation)
  normalized_likelihoods <- t(normalized_likelihoods)
  return(normalized_likelihoods)
}

update_all_seeds <- function(seeds, population, data.loop, sheet.cells) {
  #' Update the cells and the occurrences of all the seeds to account for the cells dispatched after soft-clustering.
  #' This function must be called before generating the meta and markers sheets.
  #' 
  #' @param seeds: a nested list, where each sub-list has five keys: 'clusters', 'consensus', 'cells', 'genes' and 'markers'.
  #' @param population: a character.
  #' @param data.loop: a list of four data.frames: 'expression.loop', 'occurrences.loop', 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  #' @return a nested list, where each sub-list has five keys: 'clusters', 'consensus', 'cells', 'genes' and 'markers'.
  #'
  name_subpopulation <- function(i){glue("{population}{i}")}
  update_seed <- function(i) {
    seed <- seeds[[i]]
    subpopulation <- name_subpopulation(i)
    seed$cells <- get_cells_of_interest(subpopulation, sheet.cells)
    seed$occurrences <- get_occurrences(data.loop$ranked_genes.loop[, seed$cells])
    return(seed)
  }
  
  seeds <- lapply(X=1:length(seeds), FUN=update_seed)
  return(seeds)
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
  #' @param data.loop: a list of four data.frames: 'expression.loop', 'occurrences.loop', 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param params: a list of parameters, with 'leftovers_strategy'.
  #' 
  #' @return a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  sheet.cells <- get_sheet.cells.default(records, seeds, population)
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
    
    seeds <- update_leftover_seed(seeds, population, leftover_likelihoods.current, data.loop)
    # the leftover seed is updated (cells, markers and specific_markers)
    # to account for cells displaced w.r.t. their maximum likelihood
    
    leftover_cells.previous <- leftover_cells.current
    leftover_cells.current <- seeds[[length(seeds)]]$cells
  } #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ranked_genes.init <- data.loop$ranked_genes.loop[, leftover_cells.init]
  leftover_likelihoods.init <- get_leftover_likelihoods(ranked_genes.init,
                                                        fun_likelihoods.cell,
                                                        seeds, population)
  # the likelihoods are calculated again after generating
  # a pruned set of marker genes for the leftover seed
  
  sheet.cells[leftover_cells.init,] <- leftover_likelihoods.init
  return(sheet.cells)
}