"Run this script to measure the results of a scROB analysis.

	2023/12/12 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(arules)
})

get_final_labels <- function(records.cells) {
  #' For every cell, get the final label.
  #' 
  #' @param records.cells: a binary data.frame with cluster labels as cols and cell names as rows.
  #' 
  #' @return a named factor.
  #' 
  get_final_label.cell <- function(u){
    is_in_class <- u>0
    final_index <- which.max(cumsum(is_in_class))
    final_label <- names(u)[final_index]
    return(final_label)
  }
  final_labels <- apply(X=records.cells,
                        MARGIN=1,
                        FUN=get_final_label.cell)
  final_labels <- factor(final_labels)
  return(final_labels)
}

test <- SEUROBJ.INIT
final_labels <- get_final_labels(RECORDS$cells)
test@active.ident <- final_labels
do_DimPlot(test) + do_DimPlot(SEUROBJ.INIT)
