"Functions called to generate results where the marker genes are weighted,
  w.r.t. their specificity to a population.

	2024/02/26 @yanisaspic"

f1_score <- function(TP, FN, FP) {
  #' Get the F1-score. Its formula is: 2*TP / (2*TP + FN + FP), with:
  #'
  #' @TP: a numeric. The number of cells expressing marker i within population j.
  #' @FN: a numeric. The number of cells not expressing marker i within population j.
  #' @FP: a numeric. The number of cells expressing marker i outside of population j.
  #' 
  #' @return a numeric between 0 and 1.
  #' 
  F1 <- 2*TP / (2*TP + FN + FP)
  F1 <- round(F1, 2)
  return(F1)
}

add_genes_to_f1_scores.seed <- function(f1_scores.seed, genes) {
  #' Add a set of genes to a named vector of F1 scores.
  #' Their value is 0.
  #'
  #' @param f1_scores.seed: a named vector of f1-scores. The names are genes.
  #' @param genes: a vector of characters.
  #' 
  #' @return a named vector of f1-scores. The names are genes.
  #' 
  f1_scores.seed[genes] <- 0
  # the vector is sorted alphabetically for the downstream cbind:
  f1_scores.seed <- f1_scores.seed[order(names(f1_scores.seed))]
  return(f1_scores.seed)
}

get_f1_scores.seed <- function(seed, occurrences.population) {
  #' Get the F1-scores of every marker genes of a seed. Its value is between 0 and 1, 
  #' and the higher it is, the better a marker gene is to characterize a subpopulation.
  #' 
  #' @param seed: a list with six keys: 'consensus', 'cells', 'clusters', 'markers', 'specific_markers' and 'occurrences'.
  #' @param occurrences.population: a data.frame where: genes are rows | sampling effort is cols | cells are occurrences.
  #'
  #' @return a named vector of f1-scores. The names are marker genes.
  #' 
  occ.seed <- seed$occurrences[seed$markers,]
  occ.pop <- occurrences.population[seed$markers,]
  
  TP <- occ.seed[, ncol(occ.seed)]
  data <- data.frame(TP=TP,
                     FN=length(seed$cells) - TP,
                     FP=occ.pop[, ncol(occ.pop)] - TP)
  rownames(data) <- seed$markers
  
  get_f1_score.marker <- function(row) {f1_score(row['TP'], row['FN'], row['FP'])}
  f1_scores.seed <- apply(X=data, MARGIN=1, FUN=get_f1_score.marker)
  
  non_markers <- setdiff(rownames(occurrences.population), seed$markers)
  f1_scores.seed <- add_genes_to_f1_scores.seed(f1_scores.seed, non_markers)
  return(f1_scores.seed)
}

get_sheet.markers.weighted <- function(records, seeds, population, occurrences.population) {
  #' Get a sheet of marker genes w.r.t. populations.
  #' The value i,j is the F1-score of marker i w.r.t. population j (see get_f1_score.gene).
  #' 
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param seeds: a nested list, where each sub-list has six keys: 'consensus', 'cells', 'clusters', 'markers', 'specific_markers' and 'occurrences'.
  #' @param population: a character.
  #' @param occurrences.population: a data.frame where: genes are rows | sampling effort is cols | cells are occurrences.
  #' 
  #' @return a data.frame where rows are genes | cols are populations | values are binary.
  #' 
  get_f1_scores.seed.wrapper <- function(i) {get_f1_scores.seed(seeds[[i]], occurrences.population)}
  f1_scores <- lapply(X=1:length(seeds), FUN=get_f1_scores.seed.wrapper)
  
  non_HVGs <- setdiff(rownames(records$markers), rownames(occurrences.population))
  f1_scores <- lapply(X=f1_scores, FUN=add_genes_to_f1_scores.seed, genes=non_HVGs)
  sheet.markers <- do.call(cbind, f1_scores)
  
  name_subpopulation <- function(i){glue("{population}{i}")}
  colnames(sheet.markers) <- sapply(X=1:length(seeds), FUN=name_subpopulation)
  return(sheet.markers)
}
