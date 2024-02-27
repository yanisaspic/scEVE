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

get_classification <- function(records, ground_truth=FALSE) {
  #' Get a data.frame associating every unique cell to its most informative cluster label.
  #' 
  #' @param records: a list of three data.frames: 'meta', 'cells' and 'markers'.
  #' @param ground_truth: a boolean. If TRUE, cell names must follow the {ground_truth}_{i} format.
  #' Add a column corresponding to the ground truth labels.
  #' 
  #' @return a data.frame with two columns: 'cell' and 'pred'. 
  #' If ground_truth is TRUE, there is a third column 'ground'.
  #' 
  all_cells <- rownames(records$cells)
  classification <- data.frame(cell=all_cells, pred="C")
  for (cluster in colnames(records$cells)) {
    is_in_cluster <- records$cells[, cluster]==1
    cells_in_cluster <- all_cells[is_in_cluster]
    classification[classification$cell %in% cells_in_cluster, "pred"] <- cluster
  }
  
  if (ground_truth) {
    ground_truth <- sapply(strsplit(classification$cell, split="_"), "[", 1)
    classification$ground <- ground_truth
  }
  return(classification)
}

get_piecharts.population <- function(population, classification, records) {
  #' Get the number of individuals from each ground truth group within a population.
  #'
  #' @param population: a character.
  #' @param classification: a data.frame with three columns: 'cell', 'pred' and 'ground'. 
  #' @param records: a list of three data.frames: 'meta', 'cells' and 'markers'.
  #'
  #' @return a vector of numeric.
  #' 
  cells_of_interest <- rownames(records$cells[records$cells[population]==1,])
  data <- classification[classification$cell %in% cells_of_interest,]
  
  data <- data %>% count(ground)
  piechart <- setNames(object = data$n, data$ground)
  ground_groups <- unique(classification$ground)
  missing_groups <- setdiff(ground_groups, data$ground)
  
  piechart[missing_groups] <- 0
  piechart <- piechart[ground_groups]
  return(piechart)
}

get_piecharts <- function(records, classification) {
  #' For each population, get the number of individuals composing it, split by ground truth groups.
  #' 
  #' @param records: a list of three data.frames: 'meta', 'cells' and 'markers'.
  #' @param classification: a data.frame with three columns: 'cell', 'pred' and 'ground'. 
  #' 
  #' @return a list of vectors.
  #' The vector elements are ordered and correspond to the number of individuals in each ground truth group.
  #'
  pred_groups <- colnames(records$cells)
  piecharts <- lapply(X=pred_groups,
                      FUN=get_piecharts.population,
                      records=records,
                      classification=classification)
  names(piecharts) <- (pred_groups)
  return(piecharts)
}

draw_scEVE <- function(records) {
  #' Draw a graph with pie charts as vertices.
  #' The vertices are the cell population and the edges their hierarchy.
  #' 
  #' @param records: a list of three data.frames: 'meta', 'cells' and 'markers'.
  #' 
  #' @return a graph.
  #'
  data <- data.frame(parent=records$meta$parent,
                     child=row.names(records$meta), 
                     consensus=records$meta$consensus,
                     size=records$meta$n)
  hierarchy_graph <- graph.data.frame(data[-1,], directed=TRUE)
  E(hierarchy_graph)$label <- data[-1,]$consensus
  nodes_order <- vertex_attr(hierarchy_graph)$name
  
  # set nodes pie charts w.r.t. ground truth
  ##########################################
  classification <- get_classification(records, ground_truth = TRUE)
  ground_groups <- unique(classification$ground)
  colors <- hue_pal()(length(ground_groups))
  piecharts <- get_piecharts(records, classification)
  piecharts <- piecharts[nodes_order]
  
  g <- plot(hierarchy_graph,
            layout=layout_as_tree,
            root=-1,
            vertex.shape="pie",
            vertex.pie=piecharts,
            vertex.pie.color=list(colors))
  
  return(g)
}