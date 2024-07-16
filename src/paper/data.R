"Functions used to get the scRNA-seq datasets used in the scEVE paper.

	2024/05/23 @yanisaspic"

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(SPARSim)
  library(TMExplorer) # packageVersion("TMExplorer")==1.10.0
})

#___________________________________________________________________________real

get_metadata <- function() {
  #' Get the metadata regarding the real datasets with counts data.
  #' 
  #' @return a data.frame with seven columns: 'sequencing', 'cells', 'clusters',
  #' 'genes', 'year', 'accession' and 'doi'.
  #' 
  metadata <- c(
    # Raw counts datasets still available (June 2024)
    # on TMExplorer 
    "SMARTer (Fluidigm C1)", 364, 7, 57241, 2017, "GSE81861", "10.1038/ng.3818", 
    "SMART-Seq2", 3589, 7, 23460, 2017, "GSE84465", "10.1016/j.celrep.2017.10.030",
    "SMART-Seq2", 6879, 9, 23686, 2018, "GSE115978", "10.1016/j.cell.2018.09.006",
    "10x Genomics", 18456, 18, 23580, 2020, "GSE125969", "10.1016/j.celrep.2020.108023",
    "Seq-Well", 22600, 17, 27899, 2018, "GSE116256", "10.1016/j.cell.2019.01.031",
    "10x Genomics", 51775, 17, 22180, 2018, "E-MTAB-6149, E-MTAB-6653", "10.1038/s41591-018-0096-5",
    "10x Genomics", 57530, 10, 24005, 2019, "CRA001160", "10.1038/s41422-019-0195-y",
    
    # Raw counts datasets used in the scEFSC paper still available (June 2024)
    # on https://hemberg-lab.github.io/scRNA.seq.datasets/
    "inDrop", 1937, 14, 20125, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
    "inDrop", 1724, 14, 20125, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
    "inDrop", 3605, 14, 20125, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
    "inDrop", 1303, 14, 20125, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
    "inDrop", 822, 13, 14878, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
    "inDrop", 1064, 13, 14878, 2016, "GSE84133", "10.1016/j.cels.2016.08.011",
    "SMARTer (Fluidigm C1)", 561, 9, 55186, 2017, "GSE81861", "10.1038/ng.3818",
    "SMARTer", 1679, 18, 24057, 2016, "GSE71585", "10.1038/nn.4216"
    )
  
  datasets <- c(
    # Raw counts datasets still available (June 2024)
    # on TMExplorer 
    "Li_HumCRC_b", "Darmanis_HumGBM", "JerbyArnon_HumMLM", "Gillen_HumEPDM",
    "VanGalen_HumAML", "Lambrechts_HumNSCLC", "Peng_HumPDAC",
    
    # Raw counts datasets used in the scEFSC paper still available (June 2024)
    # on https://hemberg-lab.github.io/scRNA.seq.datasets/
    "Baron_HumPan_1", "Baron_HumPan_2", "Baron_HumPan_3", "Baron_HumPan_4",
    "Baron_MouPan_1", "Baron_MouPan_2", "Li_HumCRC_a", "Tasic_MouBra"
    )
  characteristics <- c('protocol', 'cells', 'clusters', 'genes', 'year', 'accession', 'doi')
  
  metadata <- matrix(metadata, nrow=length(datasets), ncol=length(characteristics), byrow=TRUE)
  metadata <- as.data.frame(metadata, row.names=datasets)
  colnames(metadata) <- characteristics
  for (col in c("cells", "clusters", "genes", "year")) {metadata[, col] <- as.numeric(metadata[, col])}
  
  return(metadata)
}

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
  accessions <- list(
    Li_HumCRC_b="GSE81861",
    Darmanis_HumGBM="GSE84465",
    JerbyArnon_HumMLM="GSE115978",
    Gillen_HumEPDM="GSE125969",
    VanGalen_HumAML="GSE116256",
    Lambrechts_HumNSCLC="E-MTAB-6149,   E-MTAB-6653,",
    Peng_HumPDAC="CRA001160"
  )
  
  if (scRNAseq_label %in% names(accessions)) {
    expression <- get_scRNAseq_matrix.accession(accessions[[scRNAseq_label]])}
  else {expression <- read.csv(glue("./data/{scRNAseq_label}.csv"), header=TRUE, row.names=1)}
  rownames(expression) <- gsub("_", "--", rownames(expression))
  # _ in feature names can lead to errors downstream with Seurat
  return(expression)
}


#______________________________________________________________________synthetic

get_population_sizes <- function(n_populations, total_size, balanced=TRUE) {
  #' Get the sizes of cell populations in a synthetic dataset.
  #' 
  #' @param n_populations: a numeric. The number of cell populations.
  #' @param balanced: a boolean. If TRUE, every population has the same size. 
  #' If FALSE, populations follow a geometric distribution.
  #' @param total_size: a numeric. The total number of cells across all populations.
  #' 
  #' @return a vector of numeric.
  #' 
  if (balanced) {population_sizes <- rep(total_size / n_populations, n_populations)}
  else {
    if (n_populations==1) {population_sizes <- c(total_size)}
    else {
      population_sizes <- (1/2) ^ (1:n_populations) * total_size
      population_sizes[n_populations] <- population_sizes[n_populations - 1] 
    }
  }
  population_sizes <- as.integer(population_sizes)
  return(population_sizes)
}

get_fc_multipliers <- function(parameters.init, n_degs=500) {
  #' Get fold-change (fc) multiplier values to introduce DEGs to a cell population parameter set.
  #' 
  #' @param parameters.init: a named list of vectors with intensity, variability and lib_size.
  #' @param n_degs: a numeric. The number of DEGs.
  #' 
  #' @return a vector of numeric.
  #'
  underexpressed_genes_multipliers <- runif(n=n_degs/2, min=0.0001, max=0.25)
  overexpressed_genes_multipliers <- runif(n=n_degs/2, min=4, max=100)
  deg_multipliers <- c(underexpressed_genes_multipliers, overexpressed_genes_multipliers)

  n_genes <- length(parameters.init$variability)
  genes_multipliers <- rep(1, n_genes - n_degs)
  multipliers <- c(genes_multipliers, deg_multipliers)
  
  fold_change_multipliers <- sample(multipliers)
  return(fold_change_multipliers)
}

get_parameters.population <- function(n_cells, parameters.init, label) {
  #' Generate parameters specific to one cell population from existing parameters.
  #' The parameters are generated by both:
  #' - introducing genes with differential expression (DEGs).
  #' - changing the gene variability.
  #' These approaches are described in the SPARSim vignette.
  #' 
  #' @param n_cells: a numeric. The number of cells in the population.
  #' @param parameters.init: a named list of vectors with intensity, variability and lib_size.
  #' @param label: a character. The label of the population.
  #' 
  #' @return a named list of vectors with intensity, variability and lib_size.
  #' 
  fc_multipliers <- get_fc_multipliers(parameters.init)
  parameters <- SPARSim_create_DE_genes_parameter(sim_param=parameters.init,
                                                  fc_multiplier=fc_multipliers,
                                                  N_cells=n_cells,
                                                  condition_name=as.character(label))
  parameters$variability <- parameters$variability * runif(n=length(parameters$variability),
                                                           min=1, max=1.5)
  return(parameters)
}

get_parameters.related <- function(population_sizes, balanced, related, 
                                  parameters.init, total_size) {
  #' Get the parameters for each cell population of a related scRNA-seq dataset.
  #' 
  #' @param population_sizes: a vector of numeric.
  #' @param balanced: a boolean. If TRUE, every population has the same size. 
  #' If FALSE, populations follow a geometric distribution.  
  #' @param related: a boolean. If TRUE, the template parameter set is used recursively
  #' to generate new parameters: to represent related cell populations (e.g. T cells, T CD4+ cells).
  #' If false, a single template parameter set is used to generate new parameters: 
  #' to represent independent cell populations (e.g. T cells, neurons)
  #' @param parameters.init: a named list of vectors with intensity, variability and lib_size.
  #' @param total_size: a numeric. The total number of cells.
  #' 
  #' @return a named list of parameters.
  #'
  parameters.template <- parameters.init
  parameters <- list()
  for (i in 1:length(population_sizes)) {
    parameters[[i]] <- get_parameters.population(population_sizes[i],
                                                 parameters.template,
                                                 label=i)
    parameters.template <- parameters[[i]]
  }
  return(parameters)
}

get_synthetic_dataset <- function(n_populations, balanced, related, random_state,
                                  parameters.init, total_size=10000) {
  #' Get a scRNA-seq synthetic dataset.
  #' This function follows the vignette of the SPARSim simulator.
  #' 
  #' @param n_populations: a numeric. The number of cell populations.
  #' @param balanced: a boolean. If TRUE, every population has the same size. 
  #' If FALSE, populations follow a geometric distribution.  
  #' @param related: a boolean. If TRUE, the template parameter set is used recursively
  #' to generate new parameters: to represent related cell populations (e.g. T cells, T CD4+ cells).
  #' If false, a single template parameter set is used to generate new parameters: 
  #' to represent independent cell populations (e.g. T cells, neurons)
  #' @param random_state: a numeric.
  #' @param parameters.init: a named list of vectors with intensity, variability and lib_size.
  #' @param total_size: a numeric. The total number of cells.
  #'
  #' @return a named list of parameters. The names are cell types.
  #' 
  set.seed(random_state)
  get_parameters.population.wrapper <- function(i) {
    get_parameters.population(population_sizes[i], parameters.init, i)
  }
  population_sizes <- get_population_sizes(n_populations, total_size, balanced)
  
  if (related) {parameters <- get_parameters.related(population_sizes, balanced,
                                                   related, parameters.init, total_size)}
  else {parameters <- lapply(X=1:n_populations, FUN = get_parameters.population.wrapper)}
  dataset <- SPARSim_simulation(parameters)$count_matrix
  dataset <- as.data.frame(dataset)
  
  shuffled_cells <- sample(colnames(dataset), ncol(dataset))
  dataset <- dataset[, shuffled_cells]
  return(dataset)
}