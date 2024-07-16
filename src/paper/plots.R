"Functions used to generate the figures in the scEVE paper.

	2024/05/27 @yanisaspic"

suppressPackageStartupMessages({
  library(dplyr)
  library(egg)
  library(ggplot2)
  library(ggplotify)
  library(ggtree)
  library(glue)
  library(patchwork)
  library(RColorBrewer)
  library(readxl)
  library(reshape2)
  library(scales)
  library(tidytree)
})

get_prior <- function() {
  #' Get the prior knowledge on the benchmark regarding the figures of the scEVE paper.
  #' 
  #' @return a named list with: "real_datasets", "ensemble_benchmark" and "algorithms".
  #' "algorithms" is an ordered vector of algorithms benchmarked.
  #' "ensemble_algorithms" is an ordered vector of ensemble_algorithms compared.
  #' "algorithm_challenges" is a vector of algorithms addressing single-cell challenges (cf. Lähnemann et al.)
  #' "real_datasets" is a vector of real scRNA-seq datasets labels.
  #' Note that the datasets are sorted by number of cells.
  #' "ensemble_benchmark" is a data.frame with the performances of other ensemble algorithms.
  #' "colormap" associates each method to a specific color.
  #' They are reported from the scEFSC paper (Bian et al. 2022)
  #' 
  algorithms <- c("densityCut", "monocle3", "Seurat", "SHARP", "scEVE")
  ensemble_algorithms <- c("EC.PGMGR", "GRACE", "RSEC", "SAFE", "SAME", "scEFSC", "scEVE")
  challengers <- c("scEVE", "RSEC")
  
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
    "Gillen_HumEPDM",      # 18,456c
    "VanGalen_HumAML",     # 27,899c
    "Lambrechts_HumNSCLC", # 51,775c
    "Peng_HumPDAC"         # 57,530c
  )
  
  ensemble_performances.ARI <- c(
    # SAFE
    0.281, 0.8052, 0.595, 0.75, 0.7222, 0.1849, NA, NA,
    # SAME
    0.63, 0.8099, 0.73, 0.8621, 0.7461, 0.4103, 0.92, 0.78,
    # EC-PGMGR
    0.6987, 0.8976, 0.8454, 0.8903, 0.781, 0.9145, NA, NA,
    # scEFSC
    rep(NA, 4), NA, NA, 0.99, 0.84,
    # RSEC
    rep(NA, 4), NA, NA, 0.75, 0.30,
    # GRACE
    0.60, 0.86, 0.97, 0.93, 0.75, 0.50, NA, NA
  )
  ensemble_performances.NMI <- c(
    # SAFE
    0.6636, 0.7927, 0.7365, 0.7603, 0.7759, 0.5317, NA, NA,
    # SAME
    0.8103, 0.8059, 0.828, 0.8664, 0.802, 0.6253, 0.96, 0.84,
    # EC-PGMGR
    0.825, 0.8629, 0.8529, 0.876, 0.8331, 0.8306, NA, NA,
    # scEFSC
    rep(NA, 4), NA, NA, 0.98, 0.87, 0.75,
    # RSEC
    rep(NA, 4), NA, NA, 0.80, 0.61,
    # GRACE
    rep(NA, 8)
  )
  ensemble_datasets <- c(
    "Baron_HumPan_1", "Baron_HumPan_2", "Baron_HumPan_3", "Baron_HumPan_4",
    "Baron_MouPan_1", "Baron_MouPan_2", "Li_HumCRC_a", "Tasic_MouBra")
  ensemble_methods <- c("SAFE", "SAME", "EC.PGMGR", "scEFSC", "RSEC", "GRACE")
  configurations <- expand.grid(list(dataset=ensemble_datasets, method=ensemble_methods))
  for (metric in c("ARI", "NMI")) {
    configurations[, metric] <- get(glue("ensemble_performances.{metric}"))}
  
  baron_humpan.ARI <- c(median=0.54, min=0.54, max=0.54, # scEFSC
                        median=0.15, min=0.15, max=0.15) # RSEC
  baron_humpan.NMI <- c(median=0.75, min=0.75, max=0.75, # scEFSC
                        median=0.51, min=0.51, max=0.51) # RSEC
  
  colormap <- list(densityCut="#e31a1c",
                   monocle3="#33a02c",
                   Seurat="#b15928",
                   SHARP="#ff7f00",
                   EC.PGMGR="#7570b3",
                   GRACE="#e6ab02",
                   RSEC="#666666",
                   SAFE="#1b9e77",
                   SAME="#d95f02",
                   scEFSC="#e7298a",
                   scEVE="#1f78b4")
  
  prior <- list(algorithms=algorithms,
                ensemble_algorithms=ensemble_algorithms,
                challengers=challengers,
                real_datasets=real_datasets,
                colormap=colormap,
                ensemble_performances=configurations)
  return(prior)
}

get_benchmark <- function(path="./benchmark") {
  #' Merge the benchmark of each scRNA-seq dataset in a single dataframe.
  #' 
  #' @param path: a character. The path where benchmark files are stored.
  #' 
  #' @return a data.frame with ten columns: 'peakRAM', 'time', 'ARI', 'NMI',
  #' 'method', 'dataset', 'real', 'log10(s)', 'log10(Mb)' and 'method_type'.
  #' 
  individual_filenames <- list.files(path, full.names=TRUE)
  individual_benchmark <- lapply(individual_filenames, read.csv)
  benchmark <- do.call(rbind, individual_benchmark)
  
  prior <- get_prior()
  benchmark$method <- factor(benchmark$method, levels=prior$algorithms)
  
  is_real <- function(row) {row["dataset"] %in% prior$real}
  benchmark[,"real"] <- apply(X=benchmark, MARGIN=1, FUN=is_real)
  for (col in c("peakRAM", "time", "ARI", "NMI")) {benchmark[,col] <- as.numeric(benchmark[,col])}
  benchmark[, "log10(s)"] <- log10(benchmark[, "time"])
  benchmark[, "log10(Mb)"] <- log10(benchmark[, "peakRAM"])
  
  benchmark[, "method_type"] <- "individual"
  benchmark[benchmark$method=="scEVE", "method_type"] <- "scEVE"
  return(benchmark)
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

get_heatmap.real <- function(benchmark.real, metric) {
  #' Get a heatmap where rows are datasets and cols methods.
  #' The value of the cells correspond to an input metric.
  #' Note that the datasets are sorted by number of cells.
  #' 
  #' @param benchmark.real: a data.frame with three columns: 'dataset', 'method' and a metric.
  #' @param metric: a character. The metric of interest.
  #' 
  #' @return a plot.
  #' 
  fontfaces <- get_fontfaces(data=benchmark.real[, metric],
                             rowsize=length(levels(benchmark.real$method)),
                             ascending=ifelse(metric %in% c("ARI", "NMI"), TRUE, FALSE))
  # bold best values.
  
  plot <- ggplot(benchmark.real, aes(x=method, y=dataset, fill=.data[[metric]])) +
    geom_tile(color="white", lwd=.5, linetype=1) +
    geom_text(aes(label=round(.data[[metric]], 2)), fontface=fontfaces, size=4)
  
  minimum <- min(benchmark.real[, metric])
  maximum <- max(benchmark.real[, metric])
  midpoint <- median(benchmark.real[, metric])
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

get_boxplot.methods <- function(benchmark.real, metric) {
  #' Get a boxplot where the x-axis is methods and the y-axis is a performance metric.
  #' 
  #' @param benchmark.real: a data.frame with three columns: 'dataset', 'method' and a metric.
  #' @param metric: a character. The metric of interest.
  #' 
  #' @return a plot.
  #' 
  plot <- ggplot(benchmark.real, aes(x=method, y=.data[[metric]])) +
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

get_plot.real <- function(benchmark.real, metric) {
  #' Generate a composite plot of the benchmark on real scRNA-seq datasets.
  #' It is composed of a heatmap and two boxplots w.r.t. the datasets and methods used respectively.
  #' 
  #' @param benchmark.real: a data.frame with three columns: 'dataset', 'method' and a metric.
  #' @param metric: a character. The metric of interest.
  #' 
  #' @return a plot.
  #' 
  benchmark.real$dataset <- factor(benchmark.real$dataset, levels=get_prior()$real_datasets)
  heatmap.real <- get_heatmap.real(benchmark.real, metric)
  boxplot.methods <- get_boxplot.methods(benchmark.real, metric)
  plot.real <- ggarrange(boxplot.methods, heatmap.real, nrow=2, ncol=1, widths=1,
                         heights=c(3,7),byrow=TRUE)
  return(plot.real)
}

setup_benchmark.ensemble <- function(benchmark.real) {
  #' Merge the prior regarding other ensemble algorithms with the performances benchmarked of scEVE.
  #' 
  #' @param benchmark.real: a data.frame with four columns: 'dataset', 'method', 'ARI' and 'NMI'.
  #' 
  #' @return a data.frame with five columns: 'dataset', 'method', 'ARI', 'NMI' and 'challenger'.
  #' 
  prior <- get_prior()
  benchmark.ensemble <- benchmark.real[
    (benchmark.real$method=="scEVE") & 
      (benchmark.real$dataset %in% prior$ensemble_performances$dataset),
    c("ARI", "NMI", "method", "dataset")
  ]
  benchmark.ensemble <- rbind(benchmark.ensemble, prior$ensemble_performances)
  benchmark.ensemble$method <- factor(benchmark.ensemble$method, levels=prior$ensemble_algorithms)
  benchmark.ensemble$challenger <- benchmark.ensemble$method %in% prior$challengers
  return(benchmark.ensemble)
}

setup_baron_humpan.data <- function(benchmark.ensemble, metric) {
  #' Set-up the data to plot the Baron_HumPan performances. Baron_HumPan is split
  #' into four datasets. In the Baron_HumPan plot, the bar indicates the median
  #' performance, and the errorbar indicates the minimum and maximum performances.
  #' 
  #' @param benchmark.ensemble: a data.frame with four columns: 'dataset', 'method', 'ARI' and 'NMI'.
  #' @param metric: a character. One of 'ARI' or 'NMI'.
  #' 
  #' @return a data.frame with four columns: 'method', 'median', 'minimum' and 'maximum'.
  #' 
  is_baron_humpan <- grepl("Baron_HumPan", benchmark.ensemble$dataset)
  baron_humpan.data <- benchmark.ensemble[is_baron_humpan, c("method", metric)]
  
  get_data.method <- function(method) {
    data.method <- baron_humpan.data[baron_humpan.data$method==method,]
    data.method <- list(method=method, y=median(data.method[, metric]),
                        ymin=min(data.method[, metric]), ymax=max(data.method[, metric]))
    return(data.method)}
  
  baron_humpan.data <- lapply(X=unique(baron_humpan.data$method), FUN=get_data.method)
  baron_humpan.data <- do.call(rbind, baron_humpan.data)
  return(baron_humpan.data)
}

get_barplot.baron_humpan <- function(baron_humpan.data) {
  #' Get a barplot for the Baron_HumPan dataset, where the bar indicates the
  #' median performance, and the errorbar indicates the minimum and maximum performances.
  #' 
  #' @param baron_humpan.data: 
}

get_barplot.ensemble <- function(benchmark.ensemble, dataset, metric, zoom_aes=FALSE) {
  #' Get a barplot where the x-axis represents methods and the y-axis represents
  #' a clustering performance. The plot is titled after the dataset of interest.
  #' 
  #' @param benchmark.ensemble: a data.frame with four columns: 'method', 'dataset', 'ARI' and 'NMI'.
  #' @param dataset: a character.sss
  #' @param metric: a character. One of 'ARI' and 'NMI'.
  #' @param zoom_aes: a boolean. If TRUE, a different aesthetic is applied.
  #' 
  #' @return a barplot.
  #' 
  prior <- get_prior()
  data <- benchmark.ensemble[benchmark.ensemble$dataset==dataset,]
  get_color <- function(method) {ifelse(method %in% prior$challengers, "black", "na.value")}
  colors <- lapply(X=prior$ensemble_algorithms, FUN=get_color)
  names(colors) <- prior$ ensemble_algorithms
  print(colors)
  
  plot <- ggplot(data=data) +
    geom_col(aes(x=method, y=.data[[metric]], fill=method, color=method)) +
    scale_y_continuous(expand=expansion(mult=0), limits=c(0,1)) +
    scale_fill_manual(values=get_prior()$colormap) +
    scale_color_manual(values=colors) +
    ggtitle(dataset)

  plot <- plot +
    theme_classic() +
    theme(axis.title.x=element_blank(), 
          panel.grid.major=element_line(linewidth=0.5),
          strip.text=element_blank(), legend.position="none") +
    facet_grid(~challenger, scales = "free_x", space = "free_x")
  
  if (zoom_aes) {plot <- plot +
    theme(panel.background=element_rect(fill="#ebebeb"),
          panel.grid.major=element_line(colour="white"))}
  
  return(plot)
}

get_barplots.ensemble <- function(benchmark.real, metric) {
  #' Get a composite plot comparing scEVE to other ensemble algorithms in the literature.
  #' 
  #' @param benchmark.real: a data.frame with four columns: 'dataset', 'method', 'ARI' and 'NMI'.
  #' @param metric: a character. One of 'ARI' and 'NMI'.
  #' 
  #' @return a composite plot.
  #' 
  benchmark.ensemble <- setup_benchmark.ensemble(benchmark.real)
  baron_humpan.data <- setup_baron_humpan.data(benchmark.ensemble, metric)
}

get_barplots.ensemble <- function(benchmark.real) {
  #' Generate barplots w.r.t. the performances of different ensemble clustering algorithms.
  #' 
  #' @param benchmark.real: a data.frame with 'dataset', 'method', 'ARI' and 'NMI'.
  #' - datasets included are 'Li_HumCRC_a', 'Tasic_MouBra' and 'Baron_HumPan_X'.
  #' - method included is 'scEVE'.
  #' 
  #' @return a plot.
  #' 
  ensemble_benchmark <- get_prior()$ensemble_benchmark
  ensemble_datasets <- c("Li_HumCRC_a", "Tasic_MouBra", "Baron_HumPan")
  is_comparable <- (benchmark.real$method == 'scEVE') & (benchmark.real$dataset %in% ensemble_datasets)
  comparable_benchmark <- benchmark.real[is_comparable, c("ARI", "NMI", "method", "dataset")]
  data <- rbind(ensemble_benchmark, comparable_benchmark)
  data$method <- factor(data$method, levels=c("scEFSC", "SAME", "scEVE", "RSEC"))
  
  data <- melt(data, variable.name="metric")
  for (col in c("ymin", "ymax")) {data[, col] <- NA}
  baron_benchmark.scEVE <- get_baron_benchmark.scEVE(benchmark.real)
  data <- rbind(data, baron_benchmark.scEVE)
  
  plot <- ggplot(data, aes(x=dataset, y=value, fill=method)) +
    geom_col(position="dodge", color="black") +
    geom_errorbar(aes(ymin=ymin, ymax=ymax), position=position_dodge(.9), linewidth=1.5,
                  width=0) +
    facet_wrap(~metric, ncol=1,
               labeller=as_labeller(c("ARI"="ARI [\u2b08]", "NMI"="NMI [\u2b08]"))) +
    scale_y_continuous(expand=expansion(mult=0), limits=c(0,1)) +
    scale_fill_manual(values=get_prior()$colormap)
  
  plot <- plot +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
          axis.text.y=element_text(vjust=1),
          panel.grid.major=element_line(linewidth=0.5),
          strip.background=element_blank(), 
          strip.text=element_text(face="bold", size=13),
          legend.position="bottom")
  
  return(plot)
}

parse_dataset <- function(dataset) {
  #' Parse a dataset label to get a 1-row data.frame. The data.frame has 3 columns:
  #' 'populations' (a character), 'nested' (a boolean) and 'balanced' (a boolean).
  #' 
  #' @param dataset: a character.
  #' 
  #' @return a data.frame with 3 columns: 'populations', 'nested' and 'balanced'.
  #' 
  elements <- strsplit(dataset, split="_")[[1]]
  populations <- substr(elements[1], start=2, stop=nchar(elements[1]))
  balanced <- substr(elements[2], 2, 2)
  nested <- substr(elements[3], 2, 2)
  row <- data.frame(populations=as.character(populations),
                    balanced=as.logical(balanced),
                    nested=as.logical(nested))
  return(row)
}

setup_benchmark.synthetic <- function(benchmark.synthetic) {
  #' Get a slim data.frame usable by get_plot.synthetic.
  #' 
  #' @param benchmark.synthetic: a data.frame with 6 columns: 'dataset',
  #' 'method', 'ARI', 'NMI', 'log10(Mb)', 'log10(s)'.
  #' 
  #' @return a data.frame with 8 columns: 'populations', 'nested', 'balanced',
  #' 'ARI', 'NMI', 'log10(Mb)', 'log10(s)' and 'method'.
  #' 
  metadata <- lapply(X=benchmark.synthetic$dataset, FUN=parse_dataset)
  metadata <- do.call(rbind, metadata)
  benchmark.synthetic <- cbind(benchmark.synthetic, metadata)
  return(benchmark.synthetic)
}

get_plot.synthetic.method <- function(benchmark.synthetic, metric, method) {
  #' Get a composite plot of the performances of a method under different conditions
  #' with synthetic scRNA-seq datasets.
  #' 
  #' @param benchmark.synthetic: a data.frame with four columns: 'balanced',
  #' 'nested', 'populations' and a metric.
  #' @param metric: a character. The metric of interest.
  #' @param method: a character. The method of interest.
  #' 
  #' @return a plot.
  #' 
  benchmark.synthetic <- setup_benchmark.synthetic(benchmark.synthetic)
  data <- benchmark.synthetic[benchmark.synthetic$method==method,]
  
  data$nested <- factor(data$nested, levels=c(FALSE, TRUE))
  nested.labels <- c("independent", "nested")
  names(nested.labels) <- c(FALSE, TRUE)
  data$balanced <- factor(data$balanced, levels=c(TRUE, FALSE))
  balanced.labels <- c("balanced", "unbalanced")
  names(balanced.labels) <- c(TRUE, FALSE)
  
  plot <- ggplot(data=data) +
    geom_boxplot(aes(x=populations, group=populations, y=.data[[metric]]),
                 fill=get_prior()$colormap[[method]]) +
    facet_grid(nested~balanced,
               labeller=labeller(nested=nested.labels, balanced=balanced.labels))
  
  minimum <- min(benchmark.synthetic[, metric])
  plot <- plot +
    scale_y_continuous(limits=c(minimum, 1)) +
    theme(panel.grid.major=element_line(linewidth=0.5))
  
  return(plot)
}

get_plot.synthetic.scenario <- function(benchmark.synthetic, metric, nested, balanced) {
  #' Get a composite plot showing the performances of scEVE on a specific scenario
  #' w.r.t. the performances of its individual methods.
  #' 
  #' @param benchmark.synthetic: a data.frame with four columns: 'balanced',
  #' 'nested', 'populations' and a metric.
  #' @param metric: a character. The metric of interest.
  #' @param nested: a boolean.
  #' @param balanced: a boolean.
  #' 
  #' @return a plot.
  #' 
  prior <- get_prior()
  benchmark.synthetic <- setup_benchmark.synthetic(benchmark.synthetic)
  data <- benchmark.synthetic[(benchmark.synthetic$nested==nested) &
                                (benchmark.synthetic$balanced==balanced), ]
  maximum <- max(data[!is.na(data[, metric]), metric])
  minimum <- min(data[!is.na(data[, metric]), metric])
  data.sub <- data[data$method != "scEVE",]
  data.main <- data[data$method == "scEVE",]
  
  get_subplot <- function(method) {
    subplot <- ggplot(data=data.sub[data.sub$method==method,]) +
      geom_boxplot(aes(y=populations, group=populations, x=.data[[metric]]),
                   fill=prior$colormap[[method]]) +
      scale_x_continuous(limits=c(minimum, maximum)) +
      theme_classic() +
      theme(panel.background=element_rect(fill="#ebebeb"),
            panel.grid.major=element_line(colour="white", linewidth=0.5))
  }
  subplots <- lapply(X=unique(data.sub$method), FUN=get_subplot)
  subplots <- ggarrange(plots=subplots, draw=FALSE)
  
  main_plot <- ggplot(data=data.main) +
    geom_boxplot(aes(y=populations, group=populations, x=.data[[metric]]),
                 fill=prior$colormap$scEVE) +
    theme_classic() +
    scale_x_continuous(limits=c(minimum, maximum)) +
    theme(panel.grid.major=element_line(linewidth=0.5))
  
  plot.synthetic.scenario <- ggpubr::ggarrange(main_plot, subplots, nrow=2)
  
  return(plot.synthetic.scenario)
}

get_cluster_tree.data <- function(meta) {
  #' Get a tree representing the hierarchical clustering of scEVE.
  #' 
  #' @param meta: a data.frame with two columns: 'parent', 'consensus' and 'n'.
  #' 
  #' @return a phylo tree object.
  #' 
  data <- meta[-1, ]
  input <- data.frame(parent=as.character(data$parent), node=rownames(data))
  cluster_tree <- as.phylo(input)
  tree.coordinates <- fortify(cluster_tree)
  cluster_tree.data <- as_tibble(cluster_tree)
  
  for (col in c("consensus", "n")) {meta[, col] <- as.numeric(meta[, col])}
  rows_order <- rownames(meta)[order(meta$consensus, -meta$n)]
  # sort by consensus and ensure that C is the first row.
  cluster_tree.data <- cluster_tree.data[match(rows_order, cluster_tree.data$label),]
  tree.coordinates <- tree.coordinates[match(rows_order, tree.coordinates$label),]
  meta <- meta[rows_order,]
  
  for (col in c("consensus", "n")) {cluster_tree.data[, col] <- meta[, col]}
  for (col in c("isTip", "x", "y", "branch", "angle")) {
    cluster_tree.data[, col] <- tree.coordinates[, col]}
  
  cluster_tree.data$node_type <- "consensus.cluster"
  cluster_tree.data$node_type[cluster_tree.data$consensus==0] <- "root.leftover"
  
  return(cluster_tree.data)
}

get_cluster_tree <- function(cluster_tree.data) {
  #' Draw a classification tree summarizing the analysis of scEVE.
  #' 
  #' @param cluster_tree.data: a tibble object with 10 columns: 'parent', 'node',
  #' 'label', 'consensus', 'isTip', 'x', 'y', 'branch', 'angle' and 'node_type'.
  #' 
  #' @return a ggtree plot.
  #' 
  colors <- list(root.leftover="white", consensus.cluster="black")
  get_color <- function(node_type) {colors[[node_type]]}
  cluster_tree.data$colors <- sapply(X=cluster_tree.data$node_type, FUN=get_color)
  
  plot <- ggtree(cluster_tree.data) +
    geom_label(data=cluster_tree.data[-1,], aes(x=branch, label=round(consensus, 2)),
               label.size=NA, hjust="right") +
    geom_label(aes(label=label, fill=node_type, colour=node_type),
               hjust=ifelse(cluster_tree.data$isTip, "right", "middle")) +
    scale_colour_manual(values=cluster_tree.data$colors) +
    scale_fill_manual(values=list(root.leftover="#36393d", consensus.cluster="#e6e6e6")) +
    guides(fill="none", color="none")
  return(plot)
}

get_barplot.leaves <- function(distribution_leaves.data) {
  #' Get a horizontal barplot associating its leaf cluster to its cell type composition.
  #' 
  #' @param distribution_leaves.data: a data.frame with three columns: 'proportion',
  #' 'ground_truth' and 'population'.
  #'
  #' @return a plot.
  #' 
  plot <- ggplot(data=distribution_leaves.data) +
    geom_bar(aes(x=population, y=n, fill=ground_truth),
             position="fill", stat="identity", color="black") +
    theme_classic() +
    scale_y_continuous(expand=expansion(mult=0)) +
    coord_flip() + scale_fill_brewer(palette="Paired") +
    labs(y="proportion")
  return(plot)
}

get_counts.leaves <- function(distriubtion_leaves.data) {
  #' Get a horizontal barplot associating its leaf cluster to its number of cells.
  #' 
  #' @param distribution_leaves.data: a data.frame with three columns: 'proportion',
  #' 'ground_truth' and 'population'.
  #'
  #' @return a plot.
  #' 
  plot <- ggplot(data=distribution_leaves.data) +
    geom_bar(aes(x=population, y=n), stat="identity", fill="black") +
    theme_classic() +
    scale_y_continuous(expand=expansion(mult=0), trans="log10") +
    coord_flip() + scale_fill_brewer(palette="Paired") +
    labs(y="log10(n)")
  return(plot)
}

# need to rework tree to have proportion with n
# need to rework heatmap to center the colorbar
# need to check for synthetic datasets zooms

get_summary_tree <- function(cluster_tree.data, distribution_leaves.data) {
  #' Get a composite plot aligning a cluster tree to the barplots of its leaves.
  #' 
  #' @param cluster_tree.data: a tibble object with 10 columns: 'parent', 'node',
  #' 'label', 'consensus', 'isTip', 'x', 'y', 'branch', 'angle' and 'node_type'.
  #' @param distribution_leaves.data: a data.frame with three columns: 'proportion',
  #' 'ground_truth' and 'population'.
  #' 
  #' @return a composite plot.
  #' 
  cluster_tree <- get_cluster_tree(cluster_tree.data)
  cluster_tree <- cluster_tree +
    theme(plot.margin=unit(c(0, 0, 0, 0), "null"))
  barplot.leaves <- get_barplot.leaves(distribution_leaves.data)
  barplot.leaves <- barplot.leaves +
    theme(axis.line.y=element_blank(), axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "null"))
  counts.leaves <- get_counts.leaves(distribution_leaves.data)
  counts.leaves <- counts.leaves +
    theme(axis.line.y=element_blank(), axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "null"))
  summary_tree <- cluster_tree + counts.leaves + barplot.leaves +
    plot_layout(widths=c(8,2,2))
  return(summary_tree)
}
