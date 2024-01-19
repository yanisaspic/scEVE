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

while(TRUE) {

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
