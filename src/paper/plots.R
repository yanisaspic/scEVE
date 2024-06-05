"Functions used to generate the figures in the scEVE paper.

	2024/05/27 @yanisaspic"

suppressPackageStartupMessages({
  library(dplyr)
  library(egg)
  library(ggpattern)
  library(ggplot2)
  library(ggplotify)
  library(glue)
  library(RColorBrewer)
  library(readxl)
  library(reshape2)
  library(scales)
})

get_prior <- function() {
  #' Get the prior knowledge on the results regarding the figures of the scEVE paper.
  #' 
  #' @return a named list with: "real_datasets", "ensemble_results" and "algorithms".
  #' "algorithms" is an ordered list of algorithms benchmarked.
  #' "real_datasets" is a vector of real scRNA-seq datasets labels.
  #' Note that the datasets are sorted by number of cells.
  #' "ensemble" is a data.frame with the performances of other ensemble algorithms.
  #' They are reported from the scEFSC paper (Bian et al. 2022)
  #' 
  algorithms <- c("densityCut", "monocle3", "Seurat", "SHARP", "scEVE")
  real_datasets <- c(
    "Li_HumCRC_b",         # 364c
    "Li_HumCRC_a",         # 561c
    "Baron_MouPan_1",      # 822c
    "Baron_MouPan_2",      # 1,064c
    "Baron_HumPan_4",      # 1,303c
    "Tasic_MouBra",        # 1,679c
    "Baron_HumPan_2",      # 1,724c
    "Baron_HumPan_1",      # 1,937c
    "Darmanis_HumGBM",     # 3,589c
    "Baron_HumPan_3",      # 3,605c
    "JerbyArnon_HumMLM",   # 6,879c
    "VanGalen_HumAML",     # 27,899c
    "Lambrechts_HumNSCLC", # 51,775c
    "Peng_HumPDAC"         # 57,530c
    )

  ARI <- c(# Li_HumCRC_a
    0.75, 0.92, 0.79, 0.99,
    # Tasic_MouBra
    0.30, 0.78, 0.65, 0.84,
    # Baron_HumPan
    0.15, 0.51, 0.54, 0.54)
  NMI <- c(# Li_HumCRC_a
    0.8, 0.96, 0.86, 0.98,
    # Tasic_MouBra
    0.61, 0.84, 0.80, 0.87,
    # Baron_HumPan
    0.51, 0.72, 0.75, 0.75)
  methods <- c("RSEC", "SAME", "scCCESS", "scEFSC")
  datasets <- c("Li_HumCRC_a", "Tasic_MouBra", "Baron_HumPan")
  configurations <- expand.grid(list(method=methods, dataset=datasets))
  for (metric in c("ARI", "NMI")) {configurations[, metric] <- get(metric)}
  
  colormap <- list(densityCut="#e31a1c",
                   monocle3="#33a02c",
                   Seurat="#b15928",
                   SHARP="#ff7f00",
                   scEVE="#1f78b4",
                   RSEC="#f7f7f7",
                   SAME="#cccccc",
                   scCCESS="#969696",
                   scEFSC="#525252")
  
  prior <- list(algorithms=algorithms,
                real_datasets=real_datasets,
                colormap=colormap,
                ensemble_results=configurations)
  return(prior)
}

get_results <- function(path="./results") {
  #' Merge the results of each scRNA-seq dataset in a single dataframe.
  #' 
  #' @param path: a character. The path where results files are stored.
  #' 
  #' @return a data.frame with ten columns: 'peakRAM', 'time', 'ARI', 'NMI',
  #' 'method', 'dataset', 'real', 'log10(s)', 'log10(Mb)' and 'method_type'.
  #' 
  individual_filenames <- list.files(path, full.names=TRUE)
  individual_results <- lapply(individual_filenames, read.csv)
  results <- do.call(rbind, individual_results)
  
  prior <- get_prior()
  results$method <- factor(results$method, levels=prior$algorithms)

  is_real <- function(row) {row["dataset"] %in% prior$real}
  results[,"real"] <- apply(X=results, MARGIN=1, FUN=is_real)
  for (col in c("peakRAM", "time", "ARI", "NMI")) {results[,col] <- as.numeric(results[,col])}
  results[, "log10(s)"] <- log10(results[, "time"])
  results[, "log10(Mb)"] <- log10(results[, "peakRAM"])
  
  results[, "method_type"] <- "individual"
  results[results$method=="scEVE", "method_type"] <- "scEVE"
  return(results)
}

get_records <- function(path="./analysis/records") {
  #' Merge the meta records of each scRNA-seq dataset analyzed with scEVE.
  #' Only the leaf rows are kept for each dataset.
  #'
  #' @param path: a character. The path where record files are stored.
  #' 
  #' @return a data.frame with two columns: 'consensus' and 'dataset'.
  #' 
  individual_filenames <- list.files(path)
  get_dataset_label <- function(filename) {strsplit(filename, split=".xlsx")[1]}
  get_dataset_record <- function(filename) {
    record <- read_excel(glue("{path}/{filename}"), sheet="meta")
    record[, "dataset"] <- get_dataset_label(filename)
    record <- record[, c("consensus", "dataset")]
    return(record)
  }
  individual_records <- lapply(individual_filenames, FUN=get_dataset_record)
  records <- do.call(rbind, individual_records)
  records <- as.data.frame(records)
  return(records)
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
    geom_tile(color="white", lwd=1, linetype=1) +
    geom_text(aes(label=round(.data[[metric]], 2)))
  
  minimum <- min(results.real[, metric])
  maximum <- max(results.real[, metric])
  midpoint <- median(results.real[, metric])
  colors <- brewer.pal(n=3, name="RdBu") # colorblind-friendly
  plot <- plot +
    scale_fill_gradient2(low=colors[1], mid=colors[2], high=colors[3],
                         midpoint=midpoint, limits=c(minimum,maximum), 
                         breaks=c(minimum, midpoint, maximum),
                         labels=function(x) sprintf("%.2f", x)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=30, vjust=0.8, hjust=0.8),
          legend.key.height=unit(0.5, "lines"), legend.title.position="top",
          legend.position="bottom", axis.title.x=element_blank(), 
          axis.title.y=element_blank(), legend.title=element_text(hjust=0.6),
          strip.text=element_blank()) +
    facet_grid(~method_type, scales = "free_x", space = "free_x")
  
  arrow <- ifelse(metric %in% c("ARI", "NMI"), "\u2b08", "\u2b0a")
  plot$labels$fill <- paste(metric, arrow, sep=" ")
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
  
  plot <- plot + scale_fill_manual(values=get_prior()[["colormap"]]) +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          legend.position="left", panel.grid.major=element_line(linewidth=0.5),
          axis.line.x=element_blank(), axis.ticks.x=element_blank(),
          strip.text=element_blank()) +
    facet_grid(~method_type, scales = "free_x", space = "free_x")
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
  results.real$dataset <- factor(results.real$dataset, levels=get_prior()$real_datasets)
  heatmap.real <- get_heatmap.real(results.real, metric)
  boxplot.methods <- get_boxplot.methods(results.real, metric)
  plot.real <- ggarrange(boxplot.methods, heatmap.real, nrow=2, ncol=1, widths=1,
                         heights=c(3,7),byrow=TRUE)
  return(plot.real)
}

get_baron_results.scEVE <- function(results.real) {
  #' Get the median, minimum and maximum performances of scEVE on the Baron_HumPan dataset.
  #' 
  #' @param results.real: a data.frame with 'dataset', 'method', 'ARI' and 'NMI'.
  #' - datasets included are 'Baron_HumPan_X'.
  #' - method included is 'scEVE'.
  #' 
  #' @return a list with three rows: 'median', 'maximum' and 'minimum'.
  #' 
  is_baron_humpan <- grepl("Baron_HumPan", results.real$dataset)
  results.real$dataset[is_baron_humpan] <- "Baron_HumPan"
  data <- results.real[(results.real$dataset=="Baron_HumPan") & (results.real$method=="scEVE"),
                       c("method", "dataset", "ARI", "NMI")]
  
  baron_results <- list()
  for (metric in c("ARI", "NMI")) {
    row <- data.frame(method="scEVE", dataset="Baron_HumPan", metric=metric,
                      value=median(data[, metric]), ymin=min(data[, metric]),
                      ymax=max(data[, metric]))
    baron_results[[metric]] <- row
  }
  baron_results <- do.call(rbind, baron_results)
  return(baron_results)
}

get_barplots.ensemble <- function(results.real) {
  #' Generate barplots w.r.t. the performances of different ensemble clustering algorithms.
  #' 
  #' @param results.real: a data.frame with 'dataset', 'method', 'ARI' and 'NMI'.
  #' - datasets included are 'Li_HumCRC_a', 'Tasic_MouBra' and 'Baron_HumPan_X'.
  #' - method included is 'scEVE'.
  #' 
  #' @return a plot.
  #' 
  ensemble_results <- get_prior()$ensemble_results
  ensemble_datasets <- c("Li_HumCRC_a", "Tasic_MouBra", "Baron_HumPan")
  is_comparable <- (results.real$method == 'scEVE') & (results.real$dataset %in% ensemble_datasets)
  comparable_results <- results.real[is_comparable, c("ARI", "NMI", "method", "dataset")]
  data <- rbind(ensemble_results, comparable_results)
  
  data <- melt(data, variable.name="metric")
  for (col in c("ymin", "ymax")) {data[, col] <- NA}
  baron_results.scEVE <- get_baron_results.scEVE(results.real)
  data <- rbind(data, baron_results.scEVE)
  
  plot <- ggplot(data, aes(x=dataset, y=value, fill=method)) +
    geom_col(position="dodge", color="black") +
    geom_errorbar(aes(ymin=ymin, ymax=ymax), position=position_dodge(.9), width=.33) +
    facet_grid(~metric) +
    scale_y_continuous(expand=expansion(mult=0), limits=c(0,1)) +
    scale_fill_manual(values=get_prior()$colormap)
  
  plot <- plot +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.y=element_text(vjust=1),
          panel.grid.major=element_line(linewidth=0.5),
          strip.background=element_blank(), 
          strip.text=element_text(face="bold", size=13))
  
  return(plot)
}

get_consensus_plot <- function(results.real, records) {
  #' Draw a lineplot with 3 lines: 'ARI', 'NMI' and 'consensus' w.r.t. the dataset.
  #' 
  #' @param results.real: a data.frame with 'dataset', 'method', 'ARI' and 'NMI'.
  #' - method included is 'scEVE'.
  #' @param records: a data.frame with two columns: 'consensus' and 'records'.
  #' 
  #' @return a plot.
  #' 
  data <- results.real[results.real$method=="scEVE",]
  
  records$dataset <- factor(records$dataset, levels=get_prior()$real_datasets)
  plot <- ggplot(data=records, aes(x=dataset, y=consensus)) +
    geom_boxplot(fill="grey") +
    geom_point(position="jitter")
  
  plot <- plot +
    theme_classic() +
    theme(axis.text.x = element_text(angle=30, vjust=0.8, hjust=0.8),
          axis.title.x=element_blank())
  
  return(plot)
}
