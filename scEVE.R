"Main functions to run a scEVE analysis.

	2024/01/19 @yanisaspic"

source("./src/utils/misc.R")
source("./src/utils/seeds.R")
source("./src/utils/independent_cluster.R")

get_default_hyperparameters <- function() {
  #' Get the default hyperparameters for the scEVE algorithm. They include:
  #' - n_HVGs: number of HVGs used at each iteration.
  #' - min_prop_cells: minimum cells threshold to keep a seed, w.r.t. the full dataset.
  #' - min_prop_HVGs: minimum significant HVGs threshold to keep a population clustering, w.r.t. the HVGs sampled.
  #' - root_consensus: minimum consensus threshold to consider a seed.
  #' 
  #' @return a list of hyperparameters.
  #' 
  params <- list(
    n_HVGs=500,
    min_prop_cells=0.001,
    root_consensus=0.17, # 0.17: >2 methods
    min_prop_HVGs=0.5
  )
  return(params)
}

do_scEVE <- function(expression.init, 
                     clustering_methods=c("Seurat", "monocle3", "SHARP", "densityCut"), 
                     params=get_default_hyperparameters(),
                     figures=TRUE,
                     random_state=0) {
  #' Conduct a scRNA-seq cluster analysis with the scEVE algorithm.
  #' 
  #' @param expression.init: a scRNA-seq dataset of raw count expression, without selected genes:
  #' genes are rows | cells are cols.
  #' @param clustering_methods: a vector of valid method names. Currently, 8 methods are implemented:
  #' - Seurat, monocle3, SHARP, densityCut, CIDR, scLCA, scCCESS.Kmeans, scCCESS.SIMLR.
  #' @param params: a list of parameters.
  #' @param figures: a boolean. If TRUE, draw figures summarizing the iterative clustering of populations.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three data.frames: 'meta', 'cells', 'markers'.
  #' 
  
  # init
  ########
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
  } #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  # main
  ########
  population <- get_undug_population(records)
  
  while (!is.na(population)) { #++++++++++++++++++++++++++++++++++++++++++++++++
    records$meta[population, "to_dig"] <- FALSE
    
    data.loop <- trim_data(expression.init, 
                           population, 
                           records, 
                           params, 
                           figures, 
                           random_state,
                           SeurObj.init)
    
    clusters <- get_independent_clusters(data.loop,
                                         population,
                                         clustering_methods,
                                         figures,
                                         random_state)
    
    seeds <- get_seeds(expression.init,
                       data.loop,
                       clusters,
                       params,
                       records,
                       population,
                       figures)
    
    population <- get_undug_population(records)
  } #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  return(seeds)
}

expression.init <- read.csv("./results/tmp/ALL_DATA.csv", row.names=1, header=TRUE)
test <- do_scEVE(expression.init, figures=FALSE)

while (!is.na(POPULATION)) {
  
  if(STEP=="genes"){source("./src/get_genes.R")}
  merge_pdfs(POPULATION)
  
  if (length(SEEDS)>0) {
    RECORDS <- update_records(RECORDS, SEEDS, POPULATION)
    write.xlsx(RECORDS, "./results/records.xlsx", rowNames=TRUE)
    SEEDS <- list()
  }
  
  POPULATION <- get_undug_population(RECORDS)
}

f <- apply(X=RECORDS$markers, MARGIN=1, FUN=sum)
RECORDS$markers <- RECORDS$markers[as.logical(f),]
write.xlsx(RECORDS, "./results/records.xlsx", rowNames=TRUE)
