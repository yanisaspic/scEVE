"Main functions to run a scEVE analysis.

	2024/02/20 @yanisaspic"

source("./src/scEVE/trim.R")
source("./src/scEVE/genes.R")
source("./src/scEVE/seeds.R")
source("./src/scEVE/results.R")
source("./src/scEVE/clusterings.R")

get_default_hyperparameters <- function() {
  #' Get the default hyperparameters for the scEVE algorithm. They include:
  #' - n_HVGs: number of HVGs used at each iteration.
  #' - min_prop_cells: minimum cells threshold to keep a seed, w.r.t. the full dataset.
  #' - root_consensus: minimum consensus threshold to consider a seed.
  #' - clustering_methods: a vector of valid method names. Currently, 8 methods are implemented:
  #' Seurat, monocle3, SHARP, densityCut, CIDR, scLCA, scCCESS.Kmeans, scCCESS.SIMLR.
  #' - leftovers_strategy: a valid strategy to handle leftover cells. Currently, 2 strategies exist:
  #' + default: leftover cells are in the leftover seed
  #' + soft: a leftover cell is soft-clustered w.r.t. marker genes it expresses
  #' - min_likelihood: minimum likelihood expected to include a cell in a subsequent iteration
  #' 
  #' @return a list of hyperparameters.
  #' 
  params <- list(
    n_HVGs=500,
    min_prop_cells=0.001, # rare cells subpopulation: 1/1000
    root_consensus=0.17, # 0.17: >2 methods
    clustering_methods=c("Seurat", "monocle3", "SHARP", "densityCut"),
    leftovers_strategy="default",
    min_likelihood=0
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
  records <- init_records(expression.init)
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
  
  while (!is.na(population)) { #++++++++++++++++++++++++++++++++++++++++++++++++
    records$meta[population, "to_dig"] <- FALSE
    
    while (TRUE) { #============================================================
      data.loop <- trim_data(expression.init, population, records, params, 
                             figures, random_state, SeurObj.init)
      
      if (length(data.loop)==0){break()}
      clusterings <- get_clusterings(data.loop, population, params, figures,
                                     random_state)
      seeds <- get_seeds(expression.init, data.loop, clusterings, params, 
                         records, population, figures)
      
      if (length(seeds)==0){break()}
      seeds <- get_genes(data.loop, seeds, params, population, figures)
      
      break()
    } #=========================================================================
    gc()

    if (length(seeds)>0) {
      records <- update_records(records, seeds, population, data.loop, params)
      write.xlsx(records, "./records.xlsx", rowNames=TRUE)
    }
    if (figures){merge_pdfs(population)}
    
    seeds <- list()
    population <- get_undug_population(records)
  } #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  is_marker <- as.logical(apply(X=records$markers, MARGIN=1, FUN=sum))
  records$markers <- records$markers[is_marker,]
  write.xlsx(records, "./records.xlsx", rowNames=TRUE)
  
  pdf(file="./results.pdf")
  fig <- draw_scEVE(records)
  print(fig)
  dev.off()
  
  return(records)
}
