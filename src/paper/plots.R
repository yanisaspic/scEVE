"Functions used to generate the figures in the scEVE paper.

	2024/05/27 @yanisaspic"

suppressPackageStartupMessages({
  library(dplyr)
  library(egg)
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

get_fontfaces <- function(data, rowsize, ascending=TRUE) {
  #' Get a vector of fontfaces, so that the best values are in bold in a heatmap.
  #' 
  #' @param data: a vector of numeric values.
  #' @param rowsize: an integer. It corresponds to the number of cells by row in the heatmap.
  #' @param ascending: a boolean. If TRUE, the best values are the highest ones.
  #' 
  #' @return a vector of character.
  #' 
  if (!ascending) {data <- data * -1}
  table <- matrix(data, ncol=rowsize, byrow=TRUE)
  maxima <- apply(X=table, MARGIN=1, FUN=which.max)
  fontfaces <- rep("plain", length(data))
  for (i in 1:nrow(table)) {fontfaces[maxima[i] + (i-1) * rowsize] <- "bold"}
  return(fontfaces)
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
  fontfaces <- get_fontfaces(data=results.real[, metric],
                             rowsize=length(levels(results.real$method)),
                             ascending=ifelse(metric %in% c("ARI", "NMI"), TRUE, FALSE))
    # bold best values.
  
  plot <- ggplot(results.real, aes(x=method, y=dataset, fill=.data[[metric]])) +
    geom_tile(color="white", lwd=.5, linetype=1) +
    geom_text(aes(label=round(.data[[metric]], 2)), fontface=fontfaces)
  
  minimum <- min(results.real[, metric])
  maximum <- max(results.real[, metric])
  midpoint <- median(results.real[, metric])
  colors <- brewer.pal(n=3, name="RdBu") # colorblind-friendly
  if (metric %in% c("ARI", "NMI")) {gradient <- c(colors[1], colors[2], colors[3])}
  else {gradient <- c(colors[3], colors[2], colors[1])}
  
  plot <- plot +
    scale_fill_gradient2(low=gradient[1], mid=gradient[2], high=gradient[3],
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
  
  arrow <- ifelse(metric %in% c("ARI", "NMI"), "[\u2b08]", "[\u2b0a]")
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
    facet_grid(~metric, labeller=as_labeller(c("ARI"="ARI [\u2b08]", "NMI"="NMI [\u2b08]"))) +
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

setup_results.synthetic <- function(results.synthetic) {}

get_scales <- function(results.synthetic) {
  #' Get a fill and a color scales for the synthetic plot.
  #' 
  #' @param results.synthetic: a data.frame with five columns: 'dataset', 'balanced',
  #' 'nested', 'populations' and a metric.
  #' 
  #' @return a named list with two elements: 'fill' and 'color'.
  #' 
  prior <- get_prior()
  get_base <- function(method) {prior[["colormap"]][[method]]}
  base <- sapply(X=results.synthetic$method, FUN=get_base)
  
  color_scale <- base
  color_scale[results.synthetic$balanced] <- "black"
  fill_scale <- base
  fill_scale[!results.synthetic$balanced] <- "white"
  scales <- list(color=color_scale, fill=fill_scale)
  return(scales)
}

get_plot.synthetic <- function(results.synthetic, metric) {
  #' Get a composite plot of the performances of the method under different conditions
  #' with synthetic scRNA-seq datasets.
  #' 
  #' @param results.synthetic: a data.frame with five columns: 'dataset', 'balanced',
  #' 'nested', 'populations' and a metric.
  #' @param metric: a character. The metric of interest.
  #' 
  #' @return a plot.
  #' 
  results.synthetic$populations <- as.character(results.synthetic$populations)
  results.synthetic$nested <- ifelse(results.synthetic$nested, "independent", "nested")
  
  scales <- get_scales(results.synthetic)
  
  plot <- ggplot(data=results.synthetic, aes(group=balanced)) +
    geom_bar(stat="identity", position="dodge", aes(x=populations, y=.data[[metric]]),
             fill=scales$fill, color=scales$color) +
    facet_grid(nested~method, switch="y")

  plot <- plot +
    scale_y_continuous(position="right", expand=expansion(mult=0))

  return(plot)
}
