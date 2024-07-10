"Misc. functions called by multiple steps of the pipeline.

	2024/04/10 @yanisaspic"

get_resolution <- function(population) {
  #' Get the resolution level of a predicted population.
  #' 
  #' @param population: a character.
  #' 
  #' @return an integer.
  #' 
  resolution <- length(strsplit(population, split=".", fixed=TRUE)[[1]])
  return(resolution)
}

get_populations_at_resolution <- function(sheet.cells, target_resolution) {
  #' Get all the populations at a specific resolution.
  #' The root population is resolution 1, and its children populations are resolution 2, etc.
  #' 
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' @param target_resolution: an integer.
  #' 
  #' @return a vector of population labels.
  #' 
  populations <- colnames(sheet.cells)
  resolutions <- sapply(X=populations, FUN=get_resolution)
  populations_at_resolution <- populations[resolutions==target_resolution]
  return(populations_at_resolution)
}

get_max_resolution <- function(sheet.cells) {
  #' Get the maximum resolution of a scEVE clustering analysis.
  #' 
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  #' @return an integer.
  #' 
  populations <- colnames(sheet.cells)
  resolutions <- sapply(X=populations, FUN=get_resolution)
  max_resolution <- max(resolutions)
  return(max_resolution)
}

get_cells_of_interest <- function(population, sheet.cells) {
  #' Get all the cells of interest for the current iteration.
  #' It corresponds to the cells most likely to belong to a target population.
  #' 
  #' @param population: a character.
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  #' @return a vector of cell labels.
  #' 

  if (population=="C") {return(rownames(sheet.cells))}  
  # all cells belong to the root population
  
  cell_is_interesting <- function(cell_likelihoods) {
    max_likelihood <- max(cell_likelihoods)
    most_likely_population <- names(which.max(cell_likelihoods))
    is_interesting <- (most_likely_population == population) &
      (max_likelihood > 0) 
  }
  # max_likelihood threshold ensures that cells outside of 'population' are not included.
  
  target_resolution <- get_resolution(population)
  sheet.cells <- sheet.cells[, get_populations_at_resolution(sheet.cells, target_resolution)]
  cells_are_interesting <- apply(X=sheet.cells, MARGIN=1, FUN=cell_is_interesting)
  cells_of_interest <- rownames(sheet.cells[cells_are_interesting, ])
  return(cells_of_interest)
}

get_occurrences <- function(ranked_genes) {
  #' Get the genes sampled w.r.t. effort, i.e. the number of top ranking genes sampled.
  #' 
  #' @param ranked_genes: a data.frame where: ranks are rows | cells are cols | cells are genes.
  #' 
  #' @return a data.frame where: genes are rows | sampling effort is cols | cells are number of occurrences. 
  #'
  max_effort <- max(apply(X=!is.na(ranked_genes), MARGIN=2, FUN=sum))
  all_genes <- sort(unique(as.vector(ranked_genes)))
  
  get_occurrences_at_effort <- function(effort) {
    data <- ranked_genes[1:effort,]
    occurrences_at_effort <- table(as.vector(data))
    unsampled_genes <- setdiff(all_genes, names(occurrences_at_effort))
    occurrences_at_effort[unsampled_genes] <- 0
    return(occurrences_at_effort[all_genes])
  }
  
  results <- lapply(X=1:max_effort, FUN=get_occurrences_at_effort)
  occurrences <- as.data.frame(do.call(cbind, results))
  occurrences <- as.data.frame(occurrences)
  colnames(occurrences) <- 1:max_effort
  occurrences[] <- lapply(occurrences, as.numeric)
  
  return(occurrences)
}

get_records <- function(path) {
  #' Load the records of a scEVE analysis in a list.
  #' 
  #' @param path: a character.
  #' 
  #' @return a named list with: 'meta', 'cells' and 'markers'.
  #' 
  sheet_names <- getSheetNames(path)
  get_sheet <- function(sheet_name) {read.xlsx(path, sheet=sheet_name, rowNames=TRUE)}
  sheets <- sapply(X=sheet_names, FUN=get_sheet)
  return(sheets)
}
