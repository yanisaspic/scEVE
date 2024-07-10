"Functions called to report the results of a scEVE analysis.

	2024/02/20 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})
source("./src/scEVE/utils/misc.R")
source("./src/scEVE/leftovers_strategy/default.R")
source("./src/scEVE/markers_strategy/default.R")

get_existing_pdfs <- function(population, params) {
  #' Get the names of the existing pdf files (intermediate figures) w.r.t a population.
  #'
  #' @param population: a character.
  #' @param params: a named list with 'figures_dir'.
  #'
  #' @return a vector of filenames.
  #'
  files <- c()
  for (cat in c("trim", "clusterings", "seeds", "genes")) {
    filename <- glue("{params$figures_dir}/{population}_{cat}.pdf")
    if (file.exists(filename)) {
      files <- c(files, filename)
    }
  }
  return(files)
}

merge_pdfs <- function(population, params) {
  #' Merge a group of pdf files together.
  #' 
  #' @param population: a character.
  #' @param params: a named list with 'figures_dir'.
  #' 
  files <- get_existing_pdfs(population, params)
  pdf_combine(input = files, output = glue("{params$figures_dir}/{population}.pdf"))
  unlink(files)
}

get_sheet.cells <- function(records, seeds, population, data.loop, params) {
  #' Get a sheet of cell membership likelihood w.r.t. populations.
  #' The value i,j in the results indicates the likelihood of a cell i belonging to the population j.
  #'
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has four keys: 'consensus', 'cells', 'clusters' and 'markers'.
  #' @param population: a character.
  #' @param data.loop: a list of four data.frames: 'expression.loop', 'occurrences.loop', 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param params: a list of parameters, with 'leftovers_strategy'.
  #' In the scEVE JOBIM paper, 1 strategy exists:
  #' + default: leftover cells stay in the leftover seed.
  #' 
  #' @return a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  sheet.cells <- get_sheet.cells.default(records, seeds, population)
  records$cells <- cbind(records$cells, sheet.cells)
  records$cells <- apply(X=records$cells, MARGIN=c(1,2), FUN=as.numeric)
  return(records$cells)
}

get_sheet.meta <- function(records, seeds, population) {
  #' Get a sheet of metadata w.r.t. populations.
  #' 
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has five keys: 'consensus', 'cells', 'clusters', 'markers' and 'specific_markers'.
  #' @param population: a character.
  #'
  #' @return a data.frame with four columns: 'consensus', 'parent', 'n', and 'to_dig'.
  #'
  name_subpopulation <- function(i){glue("{population}.{i}")}
  get_row <- function(i) {
    subpopulation <- name_subpopulation(i)
    row <- c(consensus=seeds[[i]]$consensus, parent=population, 
                  n=length(seeds[[i]]$cells), to_dig=TRUE)
    return(row)
  }
  
  rows <- lapply(X=1:length(seeds), FUN=get_row)
  sheet.meta <- do.call(rbind, rows)
  rownames(sheet.meta) <- name_subpopulation(1:length(seeds))
  records$meta <- rbind(records$meta, sheet.meta)
  records$meta$consensus <- as.numeric(records$meta$consensus)
  return(records$meta)
}

get_sheet.markers <- function(records, seeds, population, params, occurrences.population) {
  #' Get a sheet of cell membership likelihood w.r.t. populations.
  #' The value i,j in the results indicates the likelihood of a cell i belonging to the population j.
  #'
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has four keys: 'consensus', 'cells', 'clusters' and 'markers'.
  #' @param population: a character.
  #' @param params: a list of parameters, with 'markers_strategy'. 
  #' In the scEVE JOBIM paper, 1 strategy exists:
  #' + default: markers are reported in a binary matrix. If marker i is over-represented in population j, it is 1.
  #' @param occurrences.population: a data.frame where: genes are rows | sampling effort is cols | cells are occurrences.
  #' 
  #' @return a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  sheet.markers <- get_sheet.markers.default(records, seeds, population)
  records$markers <- cbind(records$markers, sheet.markers)
  records$markers <- apply(X=records$markers, MARGIN=c(1,2), FUN=as.numeric)
  return(records$markers)
}

update_records <- function(records, seeds, population, data.loop, params) {
  #' Add the results of a loop to the records.
  #'
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has five keys: 'clusters', 'consensus', 'cells', 'genes' and 'markers'.
  #' @param population: a character.
  #' @param data.loop: a list of four data.frames: 'expression.loop', 'occurrences.loop', 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param params: a list of parameters, with 'leftovers_strategy'.
  #' In the scEVE JOBIM paper, 1 strategy exists:
  #' + default: leftover cells stay in the leftover seed.
  #'
  #' @return a named list of three data.frames: 'cells', 'meta' and 'markers'.
  #'
  sheet.cells <- get_sheet.cells(records, seeds, population, data.loop, params)
  sheet.markers <- get_sheet.markers(records, seeds, population, params, data.loop$occurrences.loop)
  sheet.meta <- get_sheet.meta(records, seeds, population)
  records <- list(cells=sheet.cells, meta=sheet.meta, markers=sheet.markers)
  return(records)
}

get_leaves <- function(sheet.cells) {
  #' Get a named vector associating each cell to its most informative cluster label.
  #' 
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  #' @return a named vector where cells are names and cluster labels are values.
  #' 

  if (is.vector(sheet.cells)) {
    leaves <- setNames(object = rep("C", length(sheet.cells)), nm = names(sheet.cells))
    leaves <- leaves[order(names(leaves))]
    return(leaves)
    # sheet.cells is a vector at the 1st resolution only:
    # all cells belong to the root population C.
  }
    
  cells <- rownames(sheet.cells)
  max_resolution <- get_max_resolution(sheet.cells)
  
  get_cells_of_interest.wrapper <- function(population) {
    cells_of_interest <- get_cells_of_interest(population, sheet.cells)
    labels <- rep(population, length(cells_of_interest))
    output <- setNames(labels, cells_of_interest)
    return(output)
  }
  get_cells_at_resolution <- function(resolution) {
    populations_at_resolution <- get_populations_at_resolution(sheet.cells, resolution)
    cells_at_resolution <- lapply(X=populations_at_resolution, FUN=get_cells_of_interest.wrapper)
    output <- unlist(cells_at_resolution)
    return(output)
  }
  
  # get the labels of every cell with a bottom-up approach.
  labels <- lapply(X=max_resolution:1, FUN=get_cells_at_resolution)
  labels <- unlist(labels)
  leaves <- labels[!duplicated(names(labels))]
  # if a cell has multiple labels, duplicates are lower resolution ones.
  
  leaves <- leaves[order(names(leaves))]
  # sort cells alphabetically to facilitate benchmarking.
  return(leaves)
}

get_leaves_at_resolution <- function(sheet.cells, resolution) {
  #' Get a named vector associating each cell to its cluster label at a given resolution.
  #' 
  #' @param sheet.cells: a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' @param resolution: a numeric. Resolution 1 corresponds to the root cluster C.
  #'
  #' @return a named vector where cells are names and cluster labels are values.
  #' 
  sheet.cells_at_resolution <- sheet.cells[, nchar(colnames(sheet.cells))<=resolution]
  leaves_at_resolution <- get_leaves(sheet.cells_at_resolution)
  return(leaves_at_resolution)
}