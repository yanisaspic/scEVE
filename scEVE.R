"Main functions to run a scEVE analysis.

	2024/01/19 @yanisaspic"

source("./src/utils/misc.R")
source("./src/utils/genes.R")
source("./src/utils/seeds.R")
source("./src/utils/independent_cluster.R")

get_default_hyperparameters <- function() {
  #' Get the default hyperparameters for the scEVE algorithm. They include:
  #' - n_HVGs: number of HVGs used at each iteration.
  #' - min_prop_cells: minimum cells threshold to keep a seed, w.r.t. the full dataset.
  #' - root_consensus: minimum consensus threshold to consider a seed.
  #' - clustering_methods: a vector of valid method names. Currently, 8 methods are implemented:
  #' Seurat, monocle3, SHARP, densityCut, CIDR, scLCA, scCCESS.Kmeans, scCCESS.SIMLR.
  #' @return a list of hyperparameters.
  #' 
  params <- list(
    n_HVGs=500,
    min_prop_cells=0.001, # rare cells subpopulation: 1/1000
    root_consensus=0.17, # 0.17: >2 methods
    clustering_methods=c("Seurat", "monocle3", "SHARP", "densityCut")
  )
  return(params)
}

do_scEVE <- function(expression.init, 
                     params=get_default_hyperparameters(),
                     figures=TRUE,
                     random_state=0) {
  #' Conduct a scRNA-seq cluster analysis with the scEVE algorithm.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param params: a list of parameters.
  #' @param figures: a boolean. If TRUE, draw figures summarizing the iterative clustering of populations.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three data.frames: 'meta', 'cells', 'markers'.
  #' 
  
  ## I N I T
  records <- init_records(expression.init, params)
  SeurObj.init <- NA

  if (figures) { #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dir.create("figures")
    SeurObj.init <- CreateSeuratObject(expression.init)
    SeurObj.init <- FindVariableFeatures(SeurObj.init, nfeatures=params$n_HVGs)
    SeurObj.init <- NormalizeData(SeurObj.init)
    SeurObj.init <- ScaleData(SeurObj.init, features=VariableFeatures(SeurObj.init))
    SeurObj.init <- RunUMAP(SeurObj.init,
                            features=VariableFeatures(SeurObj.init),
                            seed.use=random_state)
  } #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## M A I N
  population <- get_undug_population(records)
  records$meta[population, "to_dig"] <- FALSE
  
  while (!is.na(population)) { #++++++++++++++++++++++++++++++++++++++++++++++++
    
    while (TRUE) { #============================================================
      data.loop <- trim_data(expression.init, population, records, params, 
                             figures, random_state, SeurObj.init)
      
      if (length(data.loop)==0){break()}
      clusters <- get_independent_clusters(data.loop, population, params, 
                                           figures, random_state)
      seeds <- get_seeds(expression.init, data.loop, clusters, params, records,
                         population, figures)
      
      if (length(seeds)==0){break()}
      seeds <- get_genes(data.loop, seeds, params, population, figures) 
      
      break()
    } #=========================================================================

    merge_pdfs(population)
    if (length(seeds) > 0) {
      records <- update_records(records, seeds, population)
      write.xlsx(records, "./records.xlsx", rowNames=TRUE)
    }
    
    seeds <- list()
    population <- get_undug_population(records)
    records$meta[population, "to_dig"] <- FALSE
  } #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  is_marker <- as.logical(apply(X=records$markers, MARGIN=1, FUN=sum))
  records$markers <- records$markers[is_marker,]
  write.xlsx(records, "./records.xlsx", rowNames=TRUE)
  return(records)
}

expression.init <- read.csv("./data/Darmanis_GBM.csv", row.names=1, header=TRUE)
records <- do_scEVE(expression.init, figures=TRUE)