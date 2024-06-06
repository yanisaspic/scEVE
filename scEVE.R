"Main functions to run a scEVE analysis.

	2024/04/03 @yanisaspic"

suppressPackageStartupMessages({
  library(openxlsx) # packageVersion("openxlsx")==4.2.5.2
})

source("./src/scEVE/trim.R")
source("./src/scEVE/genes.R")
source("./src/scEVE/seeds.R")
source("./src/scEVE/results.R")
source("./src/scEVE/utils/misc.R")
source("./src/scEVE/clusterings.R")

get_default_hyperparameters <- function() {
  #' Get the default hyperparameters for the scEVE algorithm. They include:
  #' - n_HVGs: number of HVGs used at each iteration.
  #' - min_prop_cells: minimum n_cells threshold to consider an overlap, w.r.t. the full dataset.
  #' - root_consensus: minimum consensus threshold to consider a consensus cluster.
  #' - clustering_methods: a vector of valid method names. In the scEVE JOBIM paper, 4 methods are used:
  #' Seurat, monocle3, SHARP and densityCut.
  #'
  #' -> leftovers_strategy: a valid strategy to handle leftover cells. In the scEVE JOBIM paper, 1 strategy exists:
  #' + default: leftover cells are in the leftover seed
  #'
  #' -> markers_strategy: a valid strategy to report marker genes. In the scEVE JOBIM paper, 1 strategy exists:
  #' + default: markers are reported in a binary matrix. If marker i is over-represented in population j, it is 1.
  #'
  #' @return a list of hyperparameters.
  #'
  params <- list(
    n_HVGs=500, # see Theis et al.
    min_prop_cells=0.001, # rare cells subpopulation: 1/1000
    root_consensus=0.34, # 0.17: >2 methods overlapping with 4 total methods. (5: 0.10)
    clustering_methods=c("Seurat", "monocle3", "SHARP", "densityCut"), # see Yu et al: 4 fastest methods. (5: CIDR)
    leftovers_strategy="default",
    markers_strategy="default"
  )
  return(params)
}

scEVE.iteration <- function(expression.init, population, records, params,
                            figures, random_state, SeurObj.init, save) {
  #' Conduct an iteration of the scEVE algorithm.
  #'
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param population: a character.
  #' @param records: a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #' @param params: a list of parameters.
  #' @param figures: a boolean. If TRUE, draw figures summarizing the iterative clustering of populations.
  #' @param random_state: a numeric.
  #' @param SeurObj.init: a SeuratObject generated from the base scRNA-seq dataset.
  #' A U-MAP has already been applied on the object.
  #' @param save: if TRUE, a records file is saved.
  #'
  #' @return a named list of three data.frames: 'cells', 'markers' and 'meta'.
  #'

  while (TRUE) {
    data.loop <- trim_data(expression.init, population, records, params,
                           figures, random_state, SeurObj.init)
    if (length(data.loop)==0){break()}
    # _______________________________________if too little cells, do not try to cluster

    clusterings <- get_clusterings(data.loop, population, params, figures, random_state)
    seeds <- get_seeds(expression.init, data.loop, clusterings, params,
                       records, population, figures)
    if (length(seeds)==0){break()}
    # ______________________________if too little consensus for every seed, do not characterize

    data.loop$occurrences.loop <- get_occurrences(data.loop$ranked_genes.loop)
    seeds <- get_genes(data.loop, seeds, params, population, figures)
    if (length(seeds)==0){break()}
    # _____________________________if too little characterization for any seed, do not report

    records <- update_records(records, seeds, population, data.loop, params)
    if (save) {write.xlsx(records, "./records.xlsx", rowNames=TRUE)}
    break()
  }

  if (figures){merge_pdfs(population)}
  return(records)
}

do_scEVE <- function(expression.init, params=get_default_hyperparameters(),
                     random_state=0, figures=TRUE, save=TRUE) {
  #' Conduct a scRNA-seq clustering analysis with the scEVE algorithm.
  #'
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param params: a list of parameters.
  #' @param random_state: a numeric.
  #' @param figures: a boolean. If TRUE, draw figures summarizing the iterative clustering of populations.
  #' @param save: if TRUE, a records file is saved.
  #'
  #' @return a list with two elements:
  #' - records: a list of three data.frames: 'meta', 'cells', 'markers'. It corresponds to the generated xlsx.
  #' - labels: a named factor, where names are cells and values are cluster labels.
  #'

  #_________________________________________________________________________init
  records <- init_records(expression.init)
  if (figures) {
    dir.create("figures")
    SeurObj.init <- CreateSeuratObject(expression.init)
    SeurObj.init <- FindVariableFeatures(SeurObj.init, nfeatures=params$n_HVGs)
    SeurObj.init <- NormalizeData(SeurObj.init)
    SeurObj.init <- ScaleData(SeurObj.init, features=VariableFeatures(SeurObj.init))
    SeurObj.init <- RunUMAP(SeurObj.init,
                            features=VariableFeatures(SeurObj.init),
                            seed.use=random_state)
  } else {SeurObj.init <- NA}

  #____________________________________________________________________main loop
  population <- "C"
  while (!is.na(population)) {
    records <- scEVE.iteration(expression.init, population, records, params,
                               figures, random_state, SeurObj.init, save)
    records$meta[population, "to_dig"] <- FALSE
    population <- get_undug_population(records)
    gc()
  }

  #_______________________________________________________________________finale
  is_marker <- function(row) {sum(row) > 0}
  records$markers <- records$markers[apply(X=records$markers, MARGIN=1, FUN=is_marker),]
  if (save) {write.xlsx(records, "./records.xlsx", rowNames=TRUE)}
  results <- list(records=records, preds=factor(get_leaves(records$cells)))
  return(results)
}