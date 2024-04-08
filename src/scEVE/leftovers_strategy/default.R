"Functions called to generate cell clusters with the default strategy: 
  leftover cells are in the leftover seed.

	2024/02/20 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})

get_sheet.cells.default <- function(records, seeds, population) {
  #' Get a sheet of cell membership likelihoods w.r.t. populations.
  #' The value i,j in the results indicates the likelihood of a cell i belonging to the population j.
  #' The likelihood is binary: leftover cells are in the leftover seed.
  #' 
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has four keys: 'consensus', 'cells', 'clusters' and 'markers'.
  #' @param population: a character.
  #' 
  #' @return a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  cells <- rownames(records$cells)
  get_col <- function(i) {
    seed <- seeds[[i]]
    cluster_label <- glue("{population}{i}")
    col <- as.numeric(cells %in% seed$cells)
    return(col)
  }
  cols <- lapply(X=1:length(seeds), FUN = get_col)
  
  name_subpopulation <- function(i){glue("{population}{i}")}
  sheet.cells <- do.call(cbind, cols)
  rownames(sheet.cells) <- cells
  colnames(sheet.cells) <- name_subpopulation(1:length(seeds))
  return(sheet.cells)
}