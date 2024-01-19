"Identify seeds, i.e. groups of cells unambiguously assigned to a unique meta-clusters.

  This script expects the following GLOBAL VARIABLES:
  - RECORDS: a list of three data.frames: 'meta', 'cells' and 'markers'.
  - CLUSTERS: a data.frame: cells are rows | methods are cols | a cell is a label.
  - POPULATION: a character.
  - PARAMS: a named list of hyper-parameters, with 'min_prop_cells' and 'min_confidence'.
  - DATASET.LOOP: a subset scRNA-seq dataset: genes as rows | cells as cols.
  - SEUROBJ.LOOP: a SeuratObject generated from DATASET.LOOP.

	2024/01/10 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(arules)
})
source("./src/utils/seeds.R")

### Get a transactions matrix ###
#################################
clusters_path <- glue("./results/tmp/clusters/{POPULATION}.csv")
write.table(CLUSTERS,
            file=clusters_path,
            col.names=FALSE,
            row.names=FALSE
)
transactions <- read.transactions(clusters_path)

### Get meta-populations from independent clusters ###
######################################################
total_n_cells <- ncol(DATASET.INIT)
min_n_cells <- PARAMS[["min_prop_cells"]] * total_n_cells
loop_n_cells <- ncol(DATASET.LOOP)
min_support <- min_n_cells / loop_n_cells
consensus.max <- get_consensus.max(METHODS)

rules <- get_rules(transactions, min_support, PARAMS[["min_confidence"]])
rules <- get_symmetric_rules(rules)
hints <- get_hints(rules, consensus.max)

while(TRUE) {
  ### filter out conflictual hints ###
  ####################################
  if(length(hints)==0){break()}
  expected_consensus <- max(PARAMS$root_consensus,
                            RECORDS[["meta"]][POPULATION, "consensus"])
  hints <- filter_conflictual_hints(hints, expected_consensus)

  ### get seeds: i.e. identify unambiguous cells for each hint ###
  ################################################################
  if(length(hints)==0){break()}
  cells.hints <- get_cells(hints, CLUSTERS)
  SEEDS <- add_cells_to_hints(hints, cells.hints)
  has_cells <- function(pop){length(pop$cells)>1}
  SEEDS <- Filter(f=has_cells, x=SEEDS)
  consensi <- sapply(SEEDS, "[[", "consensus")
  SEEDS <- SEEDS[order(-consensi)]  # descending sort
  
  ### generate a leftover seed: naive approach ###
  ################################################
  if(length(SEEDS)==0){break()}
  cells.seeds <- unlist(sapply(SEEDS, "[[", "cells"))
  missing_cells <- setdiff(colnames(DATASET.LOOP), cells.seeds)
  leftover_seed <- list(
    consensus=0,
    clusters=c(),
    cells=missing_cells
  )
  SEEDS[[length(SEEDS) + 1]] <- leftover_seed
  
  ### draw the seeds ###
  ######################
  SEUROBJ.LOOP <- add_seeds(SEUROBJ.LOOP, SEEDS)
  pdf(file = glue("./results/figures/{POPULATION}_seeds.pdf"))
  seeds_plot <- do_DimPlot(SEUROBJ.LOOP,
                           split.by="seeds",
                           legend.position="none")
  print(seeds_plot)
  dev.off()
  
  STEP <- "genes"
  break()
}

### clear seed-specific memory space ###
########################################
rm(
  rules,
  hints,
  cells.hints,
  cells.seeds,
  missing_cells,
  transactions, 
  total_n_cells,
  min_n_cells,
  min_support,
  expected_consensus,
  loop_n_cells,
  consensus.max,
  clusters_path,
  seeds_plot,
  consensi,
  leftover_seed
)
gc()
