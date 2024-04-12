"Run this script to generate the cancer scRNA-seq datasets required.

	2024/04/08 @yanisaspic"

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(TMExplorer) # packageVersion("TMExplorer")==1.10.0
})

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
  unknown_labels <- c("?", "NA", "")
  cell_labels <- colData(dataset)
  cell_labels <- cell_labels[!cell_labels$label %in% unknown_labels, , drop=FALSE]
  # remove cells with missing ground truth labels
  
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

get_scRNAseq_matrix <- function(scRNAseq_label) {
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
  accessions[["Jerby-Arnon_HumMLM"]] <- "GSE115978"
  
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
    return(get_scRNAseq_matrix.accession(accessions[[scRNAseq_label]]))
  }
  return(glue("./data/{scRNAseq_label}.csv"))
}