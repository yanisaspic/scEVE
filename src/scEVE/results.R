"Functions called to report the results of a scEVE analysis.

	2024/02/20 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})
source("./src/scEVE/utils/misc.R")
source("./src/scEVE/leftovers_strategy/default.R")
source("./src/scEVE/leftovers_strategy/soft.R")
source("./src/scEVE/markers_strategy/default.R")
source("./src/scEVE/markers_strategy/weighted.R")

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

get_sheet.cells <- function(records, seeds, population, data.loop, params) {
  #' Get a sheet of cell membership likelihood w.r.t. populations.
  #' The value i,j in the results indicates the likelihood of a cell i belonging to the population j.
  #'
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has four keys: 'consensus', 'cells', 'clusters' and 'markers'.
  #' @param population: a character.
  #' @param data.loop: a list of four data.frames: 'expression.loop', 'occurrences.loop', 'SeurObj.loop', and 'ranked_genes.loop'.
  #' @param params: a list of parameters, with 'leftovers_strategy'.
  #' Currently, 2 strategies exist:
  #' + default: leftover cells stay in the leftover seed.
  #' + naive: leftover cells are soft-clustered w.r.t. markers they express, regardless of their expression level.
  #' 
  #' @return a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  if (params$leftovers_strategy=="default") {sheet.cells <- 
    get_sheet.cells.default(records, seeds, population)}
  else {sheet.cells <- 
    get_sheet.cells.soft(records, seeds, population, data.loop, params)}
  
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
  name_subpopulation <- function(i){glue("{population}{i}")}
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
  #' @param params: a list of parameters, with 'markers_strategy'. Currently, 2 strategies exist:
  #' + default: markers are reported in a binary matrix. If marker i is over-represented in population j, it is 1.
  #' + weighted: markers are reported with a value between 0 and 1. It corresponds to the F1 score.
  #' @param occurrences.population: a data.frame where: genes are rows | sampling effort is cols | cells are occurrences.
  #' 
  #' @return a data.frame where rows are cells | cols are populations | values are membership likelihood.
  #' 
  if (params$markers_strategy=="default") {sheet.markers <- 
    get_sheet.markers.default(records, seeds, population)}
  else {sheet.markers <- 
    get_sheet.markers.weighted(records, seeds, population, occurrences.population)}
  
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
  #' @param params: a list of parameters with 'leftovers_strategy'.
  #' Currently, 2 strategies exist:
  #' + default: leftover cells stay in the leftover seed.
  #' + naive: leftover cells are soft-clustered w.r.t. markers they express, regardless of their expression level.
  #'
  #' @return a named list of three data.frames: 'cells', 'meta' and 'markers'.
  #'
  sheet.cells <- get_sheet.cells(records, seeds, population, data.loop, params)
  
  if (params$leftovers_strategy != "default") {
    seeds <- update_all_seeds(seeds, population, data.loop, sheet.cells)
    draw_seeds(data.loop, seeds, population)
    draw_genes(data.loop, seeds, population)
  }
  # leftover cells have been soft-clustered and some cells have been displaced;
  # the cells and the occurrences of each seed must be updated.
  
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
  cells <- rownames(sheet.cells)
  labels <- colnames(sheet.cells)
  max_resolution <- max(nchar(labels))
  
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
  return(leaves)
}
