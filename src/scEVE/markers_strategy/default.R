"Functions called to generate markers with the default strategy:
  if a marker is statistically over-represented, it is encoded in a binary, regardless of the effect.

	2024/02/20 @yanisaspic"

get_sheet.markers.default <- function(records, seeds, population) {
  #' Get a sheet of marker genes w.r.t. populations.
  #' The value i,j is 1 if the gene i is marker for the population j. It is 0 otherwise.
  #' 
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has four keys: 'consensus', 'cells', 'clusters' and 'markers'.
  #' @param population: a character.
  #' 
  #' @return a data.frame where rows are genes | cols are populations | values are binary.
  #' 
  genes <- rownames(records$markers)
  frame_markers.seed <- function(seed) {as.numeric(genes %in% seed$markers)}
  sheet.markers <- lapply(X=seeds, FUN=frame_markers.seed)
  sheet.markers <- do.call(cbind, sheet.markers)
  
  name_subpopulation <- function(i){glue("{population}.{i}")}
  colnames(sheet.markers) <- sapply(X=1:length(seeds), FUN=name_subpopulation)
  rownames(sheet.markers) <- genes
  return(sheet.markers)
}
