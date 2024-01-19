"Trim a scRNA-seq 'DATASET.INIT' to generate a subset scRNA-seq 'DATASET.LOOP'.
  - the n most highly variable genes are kept.
  - the cells in the population of interest are kept.
  
  This script expects the following GLOBAL VARIABLES:
  - DATASET.INIT: a complete scRNA-seq dataset: genes as rows | cells as cols.
  - PARAMS: a named list of hyper-parameters, with 'n_HVGs'.
  - POPULATION: a character label of a group of cells.
  - SEED: a numeric random seed.
  - SEUROBJ.INIT: a SeuratObject with a U-MAP applied.
  
	2024/01/10 @yanisaspic"

set.seed(SEED)

suppressPackageStartupMessages({
  library(glue)
  library(Seurat)
  library(SCpubr)
  library(ggplot2)
  library(cowplot)
})
source("./src/utils/misc.R")

cells <- RECORDS[["cells"]]

### plot the cells of interest ###
##################################
cells_of_interest <- rownames(cells)[cells[, POPULATION] > 0.95]  # 0.95: soft-clustering
trim_plot <- do_DimPlot(SEUROBJ.INIT, cells.highlight=cells_of_interest)
pdf(file = glue("./results/figures/{POPULATION}_trim.pdf"))
print(trim_plot)
dev.off()

### trim genes if downstream analysis is relevant (i.e. enough cells) ###
#########################################################################
while(TRUE) {
  # use the default hyperparameters of the packages
  # to identify a minimum number of cells
  if (length(cells_of_interest) < 100) {
    # SHARP(sncells=100)
    print(trim_plot)
    dev.off()
    break()
  }
  
  # get the trimmed dataset #
  ###########################
  DATASET.LOOP <- DATASET.INIT[, cells_of_interest]
  HVGs.loop <- get_n_HVGs(DATASET.LOOP, n=PARAMS[["n_HVGs"]])
  DATASET.LOOP <- DATASET.LOOP[HVGs.loop, ]
  
  # set-up for SeuratObject-dependent clustering algorithms #
  ###########################################################
  SEUROBJ.LOOP <- CreateSeuratObject(DATASET.LOOP)
  VariableFeatures(SEUROBJ.LOOP) <- HVGs.loop
  SEUROBJ.LOOP <- NormalizeData(SEUROBJ.LOOP)
  SEUROBJ.LOOP <- ScaleData(SEUROBJ.LOOP, 
                            features=VariableFeatures(SEUROBJ.LOOP))
  SEUROBJ.LOOP <- RunUMAP(SEUROBJ.LOOP, 
                          features=VariableFeatures(SEUROBJ.LOOP), 
                          seed.use=SEED)
  
  # rank the genes expression cell-wise #
  #######################################
  RANKED_GENES.LOOP <- get_ranked_genes(DATASET.LOOP)
  RANKED_GENES.LOOP <- filter_expressed_genes(RANKED_GENES.LOOP, DATASET.LOOP)
  
  STEP <- "independent_cluster"
  break()
}

# clear trim-specific memory space #
####################################
rm(cells,
   cells_of_interest,
   HVGs.loop,
   trim_plot
)
gc()
