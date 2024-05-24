"Functions used to get the scRNA-seq datasets used in the scEVE paper.

	2024/05/23 @yanisaspic"

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(SPARSim)
  library(TMExplorer) # packageVersion("TMExplorer")==1.10.0
})


#___________________________________________________________________________real

get_cell_ids <- function(cell_labels) {
  #' Get a vector of unique cell ids from a vector of cell labels.
  #' 
  #' @param cell_labels: a vector of characters.
  #' 
  #' @return a vector of characters.
  #' 
  cell_ids <- c()
  label_counter <- list()
  for (label in unique(cell_labels)) {label_counter[[label]] <- 1}
  for (n in 1:length(cell_labels)) {
    cell_label <- cell_labels[n]
    cell_ids[n] <- glue("{cell_label}_{label_counter[[cell_label]]}")
    label_counter[[cell_label]] <- label_counter[[cell_label]] + 1
  }
  return(cell_ids)
}

get_scRNAseq_matrix.dataset <- function(dataset) {
  #' Get a scRNA-seq matrix from a TMExplorer dataset.
  #' 
  #' @param dataset: a dataset queried with TMExplorer, with a ground truth and raw counts.
  #' 
  #' @return a scRNA-seq raw count matrix where columns are cells, and 
  #' all cells are associated to a unique id, with the following structure: {label}_{n}.
  #' The label corresponds to the ground truth of the authors.
  #' 
  removed_labels <- c("?", "NA", "", "Doublets")
  cell_labels <- colData(dataset)
  cell_labels <- cell_labels[!cell_labels$label %in% removed_labels, , drop=FALSE]
  # remove cells with missing ground truth labels or doublets
  
  cell_order <- rownames(cell_labels)
  counts <- counts(dataset)
  counts <- counts[, cell_order]
  # align the cells of the counts and the labels
  
  if (!identical(rownames(cell_labels), colnames(counts))) {
    stop("The cell ids of the counts and the labels are not aligned.")
  }
  
  cell_ids <- get_cell_ids(cell_labels$label)
  colnames(counts) <- cell_ids
  counts <- data.frame(counts)
  return(counts)
}

get_scRNAseq_matrix.accession <- function(accession) {
  #' Get a standard scRNA-seq matrix from an accession number recognized by TMExplorer.
  #' 
  #' @param accession: a character.
  #' 
  #' @return a scRNA-seq raw count matrix where columns are cells, and 
  #' all cells are associated to a unique id, with the following structure: {label}_{n}.
  #' The label corresponds to the ground truth of the authors.
  #' 
  dataset <- queryTME(geo_accession = accession)[[1]]
  scRNAseq_matrix <- get_scRNAseq_matrix.dataset(dataset)
  
  #_____________________________________________________________________Darmanis
  if (accession=="GSE84465") {scRNAseq_matrix <- head(scRNAseq_matrix, -5)}
  # drop the rows 'no_feature', 'ambiguous', 'too_low_aQual', 
  # 'not_aligned' and 'alignment_not_unique'.
  
  return(scRNAseq_matrix)
}

get_real_dataset <- function(scRNAseq_label) {
  #' Get a standard scRNA-seq matrix from a label, e.g. 'Darmanis_HumGBM'.
  #' 
  #' @param scRNAseq_label: a character.
  #' 
  #' @return a scRNA-seq raw count matrix where columns are cells, and 
  #' all cells are associated to a unique id, with the following structure: {label}_{n}.
  #' The label corresponds to the ground truth of the authors.
  #' 
  accessions <- list()
  
  # Li (2017)___________________________________________________________________
  # accession: GSE81861
  # cells: 364
  # genes: 57,241
  # clusters: 7
  # sequencing: Fluidigm C1
  # doi: 10.1038/ng.3818
  accessions[["Li_HumCRC_b"]] <- "GSE81861"
  
  # Jerby-Arnon (2018)__________________________________________________________
  # accession: GSE115978
  # cells: 6,879
  # genes: 23,686
  # clusters: 9
  # sequencing: SMART-Seq2
  # doi: 10.1016/j.cell.2018.09.006
  accessions[["JerbyArnon_HumMLM"]] <- "GSE115978"
  
  # Van Galen (2018)___________________________________________________________
  # accession: GSE116256
  # cells: 22,600
  # genes: 27,899
  # clusters: 17
  # sequencing: Seq-Well
  # doi: 10.1016/j.cell.2019.01.031
  accessions[["VanGalen_HumAML"]] <- "GSE116256"
  
  # Lambrechts (2018)___________________________________________________________
  # accession: E-MTAB-6149, E-MTAB-6653
  # cells: 51,775
  # genes: 22,180
  # clusters: 17
  # sequencing: 10x Genomics
  # doi: 10.1038/s41591-018-0096-5 
  accessions[["Lambrechts_HumNSCLC"]] <- "E-MTAB-6149,   E-MTAB-6653,"
  
  # Peng (2019)_________________________________________________________________
  # accession: CRA001160
  # cells: 57,530
  # genes: 24,005
  # clusters: 10
  # sequencing: 10x Genomics
  # doi: 10.1038/s41422-019-0195-y
  accessions[["Peng_HumPDAC"]] <- "CRA001160"
  
  # Darmanis (2017)_____________________________________________________________
  # accession: GSE84465
  # cells: 3,589
  # genes: 23,460
  # clusters: 7
  # sequencing: SMART-Seq2
  # doi: 10.1016/j.celrep.2017.10.030
  accessions[["Darmanis_HumGBM"]] <- "GSE84465"
  
  if (scRNAseq_label %in% names(accessions)) {
    expression <- get_scRNAseq_matrix.accession(accessions[[scRNAseq_label]])}
  else {expression <- read.csv(glue("./data/{scRNAseq_label}.csv"), header=TRUE, row.names=1)}
  rownames(expression) <- gsub("_", "--", rownames(expression))
  # _ in feature names can lead to errors downstream
  return(expression)
}


#______________________________________________________________________synthetic

get_cell_distribution <- function(n_populations, distribution, size) {
  #' Get a cell distribution. It corresponds to the number of cells in each population.
  #' 
  #' @param n_populations: a numeric. The number of cell populations.
  #' @param distribution: a character. The distribution of cell types: 'uniform' or 'geometric'.
  #' @param size: a numeric. The total number of cells.
  #' 
  #' @return a vector of numeric.
  #' 
  if (distribution == "uniform") {cell_distribution <- rep(size / n_populations, n_populations)}
  if (distribution=="geometric") {
    cell_distribution <- (1/2) ^ (1:n_populations) * size
    cell_distribution <- as.integer(cell_distribution)
    cell_distribution[n_populations] <- cell_distribution[n_populations - 1]
  }
  return(cell_distribution)
}

get_parameters.population <- function(n_cells, n_degs, parameters.init, label) {
  #' Generate parameters specific to one cell population from existing parameters.
  #' The parameters are generated by multiplying the existing parameters to
  #' artificially generate genes with differential expression (DEGs).
  #' 
  #' @param n_cells: a numeric. The number of cells in the population.
  #' @param parameters.init: a named list of vectors with intensity, variability and lib_size.
  #' @param label: a character. The label of the population.
  #' @param n_degs: a numeric. It corresponds to the the number of DEGs expected,
  #' w.r.t. the initial parameters input.
  #' 
  #' @return a named list of vectors with intensity, variability and lib_size.
  #' 
  underexpressed_genes_multipliers <- runif(n=n_degs/2, min=0.0001, max=0.25)
  overexpressed_genes_multipliers <- runif(n=n_degs/2, min=4, max=100)
  deg_multipliers <- c(underexpressed_genes_multipliers, overexpressed_genes_multipliers)
  
  n_genes <- length(parameters.init$intensity)
  genes_multipliers <- rep(1, n_genes - n_degs)
  multipliers <- c(genes_multipliers, deg_multipliers)
  
  fold_change_multipliers <- sample(multipliers)
  parameters <- SPARSim_create_DE_genes_parameter(sim_param=parameters.init,
                                                  fc_multiplier=fold_change_multipliers,
                                                  N_cells = n_cells,
                                                  condition_name = as.character(label))
  return(parameters)
}

get_synthetic_dataset <- function(n_populations, distribution, random_state, parameters.init,
                                  n_degs=500, size=10000) {
  #' Get a scRNA-seq synthetic dataset.
  #' This function follows the vignette of the SPARSim simulator.
  #' 
  #' @param n_populations: a numeric. The number of cell populations.
  #' @param distribution: a character. The distribution of cell types: 'uniform' or 'geometric'.
  #' @param random_state: a numeric.
  #' @param parameters.init: a named list of vectors with intensity, variability and lib_size.
  #' @param n_degs: a numeric. It corresponds to the the number of DEGs expected,
  #' w.r.t. the initial parameters input.
  #' @param size: a numeric. The total number of cells.
  #'
  #' @return a named list of parameters. The names are cell types.
  #' 
  set.seed(random_state)
  cell_distribution <- get_cell_distribution(n_populations, distribution, size)
  get_parameters.population.wrapper <- function(i) {
    get_parameters.population(cell_distribution[i], n_degs, parameters.init, i)
  }
  parameters <- lapply(X=1:n_populations, FUN = get_parameters.population.wrapper)
  dataset <- SPARSim_simulation(parameters)$count_matrix
  shuffled_cells <- sample(colnames(dataset), ncol(dataset))
  dataset <- dataset[, shuffled_cells]
  dataset <- as.data.frame(dataset)
  return(dataset)
}
