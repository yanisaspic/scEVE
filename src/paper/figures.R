"Functions used to generate the figures in the scEVE paper.

	2024/05/27 @yanisaspic"

suppressPackageStartupMessages({
  library(egg)
  library(ggplot2)
  library(ggplotify)
  library(glue)
  library(RColorBrewer)
  library(scales)
})

get_results.prior <- function() {
  #' Get the prior knowledge on the results regarding the figures of the scEVE paper.
  #' 
  #' @return a named list with: "real_datasets", "ensemble_results" and "algorithms".
  #' "real_datasets" is a vector of real scRNA-seq datasets labels.
  #' Note that the datasets are sorted by number of cells.
  #' "ensemble" is a data.frame with the performances of other ensemble algorithms.
  #' "algorithms" is an ordered list of algorithms benchmarked.
  #' 
  results.prior <- list(
    real_datasets = c("Li_HumCRC_b",
             "Li_HumCRC_a",
             "Baron_MouPan_1",
             "Baron_MouPan_2",
             "Baron_HumPan_4",
             "Tasic_MouBra",
             "Baron_HumPan_2",
             "Baron_HumPan_1",
             "Darmanis_HumGBM",
             "Baron_HumPan_3",
             "JerbyArnon_HumMLM",
             "VanGalen_HumAML",
             "Lambrechts_HumNSCLC",
             "Peng_HumPDAC"),
    ensemble_results = c(),
    algorithms = c("densityCut", "monocle3", "Seurat", "SHARP", "scEVE")
  )
  return(results.prior)
}

get_results <- function(path="./results") {
  #' Merge the results of each scRNA-seq dataset in a single dataframe.
  #' 
  #' @param path: a character. The path where results files are stored.
  #' 
  #' @return a data.frame with eight columns: 'peakRAM', 'time', 'ARI', 'NMI',
  #' 'method', 'dataset', 'real' and 'ARI*'.
  #' The 'ARI*' is a modified ARI ranging from 0 to 1 instead of -1 to 1.
  #' 
  individual_filenames <- list.files(path, full.names=TRUE)
  individual_results <- lapply(individual_filenames, read.csv)
  results <- do.call(rbind, individual_results)
  
  results.prior <- get_results.prior()
  results$method <- factor(results$method, levels=results.prior$algorithms)
  results$dataset <- factor(results$dataset, levels=results.prior$real_datasets)
  
  is_real <- function(row) {row["dataset"] %in% results.prior$real}
  results[,"real"] <- apply(X=results, MARGIN=1, FUN=is_real)
  for (col in c("peakRAM", "time", "ARI", "NMI")) {results[,col] <- as.numeric(results[,col])}
  for (col in c("peakRAM", "time")) {results[,glue("{col}**")] <- log10(results[,col])}
  return(results)
}

get_heatmap.real <- function(results.real, metric) {
  #' Get a heatmap where rows are datasets and cols methods.
  #' The value of the cells correspond to an input metric.
  #' Note that the datasets are sorted by number of cells.
  #' 
  #' @param results.real: a data.frame with three columns: 'dataset', 'method' and a metric.
  #' @param metric: a character. The metric of interest.
  #' 
  #' @return a plot.
  #' 
  plot <- ggplot(results.real, aes(x=method, y=dataset, fill=.data[[metric]])) +
    geom_tile() +
    geom_text(aes(label=round(.data[[metric]], 2))) +
    coord_fixed()
  
  minimum <- min(results.real[, metric])
  maximum <- max(results.real[, metric])
  midpoint <- median(results.real[, metric])
  colors <- brewer.pal(n=3, name="RdBu") # colorblind-friendly
  plot <- plot +
    scale_fill_gradient2(low=colors[1], mid=colors[2], high=colors[3],
                         midpoint=midpoint, limits=c(minimum,maximum)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=30, vjust=0.9, hjust=1),
          legend.key.width=unit(0.5, "lines"), legend.position="left",
          axis.title.x=element_blank(), axis.title.y=element_blank())
  return(plot)
}

get_boxplot.datasets <- function(results.real, metric) {
  #' Get a boxplot where the x-axis is datasets and the y-axis is a performance metric.
  #' 
  #' @param results.real: a data.frame with three columns: 'dataset', 'method' and a metric.
  #' @param metric: a character. The metric of interest.
  #' 
  #' @return a plot.
  #' 
  plot <- ggplot(results.real, aes(x=dataset, y=.data[[metric]])) +
    geom_boxplot(outlier.shape=NA) +
    coord_flip() +
    geom_point(aes(fill=method), color="black", size=3, pch=21)
  
  colors <- brewer.pal(n=length(unique(results.real$method)), name="Dark2") # colorblind-friendly
  plot <- plot + scale_fill_manual(values=colors) +
    theme_classic() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.text.x=element_text(hjust=1), panel.grid.major=element_line(linewidth=0.5),
          axis.line.y=element_blank(), axis.ticks.y=element_blank()) +
    guides(fill=guide_legend(ncol=2))
  return(plot)
}

get_boxplot.methods <- function(results.real, metric) {
  #' Get a boxplot where the x-axis is methods and the y-axis is a performance metric.
  #' 
  #' @param results.real: a data.frame with three columns: 'dataset', 'method' and a metric.
  #' @param metric: a character. The metric of interest.
  #' 
  #' @return a plot.
  #' 
  plot <- ggplot(results.real, aes(x=method, y=.data[[metric]])) +
    geom_boxplot(aes(fill=method)) +
    geom_point()
  
  colors <- brewer.pal(n=length(unique(results.real$method)), name="Dark2") # colorblind-friendly
  plot <- plot + scale_fill_manual(values=colors) +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          legend.position="left", axis.text.y=element_text(vjust=1),
          panel.grid.major=element_line(linewidth=0.5),
          axis.line.x=element_blank(), axis.ticks.x=element_blank()) +
    guides(fill="none")
  return(plot)
}

get_plot.real <- function(results.real, metric) {
  #' Generate a composite plot of the benchmark on real scRNA-seq datasets.
  #' It is composed of a heatmap and two boxplots w.r.t. the datasets and methods used respectively.
  #' 
  #' @param results.real: a data.frame with three columns: 'dataset', 'method' and a metric.
  #' @param metric: a character. The metric of interest.
  #' 
  #' @return a plot.
  #' 
  heatmap.real <- get_heatmap.real(results.real, metric)
  boxplot.methods <- get_boxplot.methods(results.real, metric)
  boxplot.datasets <- get_boxplot.datasets(results.real, metric)
  legend <- cowplot::get_legend(boxplot.datasets)
  legend <- as.ggplot(legend)
  boxplot.datasets <- boxplot.datasets + guides(fill="none")
  
  plot.real <- ggarrange(boxplot.methods, legend, heatmap.real, boxplot.datasets,
                         nrow=2, ncol=2, widths=c(6,4), heights=c(4,6))
  return(plot.real)
}
