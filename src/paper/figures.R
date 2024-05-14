"Functions used to generate the figures in the scEVE paper.

	2024/05/13 @yanisaspic"

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})

get_results <- function(path="./results") {
  #' Merge the results of each scRNA-seq dataset in a single dataframe.
  #' 
  #' @param path: a character.
  #' 
  #' @return a data.frame with six columns: 'peakRAM', 'time', 'ARI', 'NMI',
  #' 'method' and 'dataset'.
  #' 
  individual_filenames <- list.files(path, full.names=TRUE)
  individual_results <- lapply(individual_filenames, read.csv)
  results <- do.call(rbind, individual_results)
  for (col in c("peakRAM", "time", "ARI", "NMI")) {results[,col] <- as.numeric(results[,col])}
  return(results)
}

get_summary_barplot <- function(results, metric, log10=FALSE) {
  #' Generate a barplot summarizing the performances on a target metric across
  #' all methods and datasets.
  #' 
  #' @param results: a data.frame with six columns: 'peakRAM', 'time', 'ARI', 'NMI',
  #' 'method' and 'dataset'
  #' @param metric: a character. One of 'peakRAM', 'time', 'ARI', 'NMI.
  #' @param log10: a boolean. If TRUE, the value of the target metric is log10-transformed.
  #' 
  #' @return a ggplot2 barplot.
  #' 
  if (log10) {results[, metric] <- log10(results[, metric])}
  fill_values <- brewer.pal(n=length(unique(results$method)), name="Dark2")
    # Dark2 is color-blind friendly
  
  plot <- ggplot(results, aes(x=method, fill=method, y=ARI), color="black", size=2) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values=fill_values) +
    scale_y_continuous(expand=expansion(mult=0),
                       limits=c(min(0, results[, metric]), max(1, results[, metric]))
                       ) +
    facet_wrap(~ dataset)
      # bars start on the x-axis instead of floating around.
  return(plot)
}
