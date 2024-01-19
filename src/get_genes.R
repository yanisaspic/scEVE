"Identify characteristic genes of a cell population.

  This script expects the following GLOBAL VARIABLES:
  - SEEDS: a nested list, where each sub-list has three keys: 'consensus', 'cells' and 'members'.
  - RANKED_GENES.LOOP: a subset scRNA-seq dataset: genes as rows | cells as cols.
  - PARAMS: a named list of hyper-parameters, with 'p-value'.

	2024/01/10 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(grid)
})

source("./src/utils/genes.R")

### Effort-Dependent enrichment Analysis: EDeA ###
##################################################
SEEDS <- add_genes_to_seeds(RANKED_GENES.LOOP, SEEDS)
occurrences.loop <- get_occurrences(RANKED_GENES.LOOP)
MARKERS <- get_markers(SEEDS, occurrences.loop, PARAMS)
SEEDS <- add_specific_markers(SEEDS, MARKERS)

### remove uninformative seeds ###
##################################
HVGs_used <- unique(MARKERS$all)
if (length(HVGs_used) < PARAMS[["min_prop_HVGs"]] * PARAMS[["n_HVGs"]]){SEEDS <- list()}

if (length(SEEDS)>0) {
  
  # boxplot of sampling effort, and upsetplot of markers, per seed #
  ##################################################################
  efforts.frame <- get_efforts.frame(RANKED_GENES.LOOP, SEEDS)
  efforts.plot <- get_efforts.plot(efforts.frame)
  markers.plot <- get_markers.plot(MARKERS)

  pdf(file = glue("./results/figures/{POPULATION}_genes.pdf"))
  composite_genes_plot <- do.call(grid.arrange, list(efforts=efforts.plot, markers=markers.plot))
  print(composite_genes_plot)
  dev.off() 
  
  STEP <- "consensus_cluster"
}

### clear genes-specific memory space ###
#########################################
rm(
  efforts.frame,
  efforts.plot,
  occurrences.loop,
  markers.plot,
  composite_genes_plot
)
gc()
