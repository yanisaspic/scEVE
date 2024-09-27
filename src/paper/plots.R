"Functions used to generate the figures in the scEVE paper.

	2024/05/27 @yanisaspic"

suppressPackageStartupMessages({
  library(dplyr)
  library(egg)
  library(ggh4x)
  library(ggplot2)
  library(ggplotify)
  library(ggtree)
  library(ggVennDiagram)
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
  #' @return a named list.
  #' "algorithms" is an ordered vector of algorithms benchmarked.
  #' "ensemble_algorithms" is an ordered vector of ensemble_algorithms compared.
  #' "algorithm_challenges" is a vector of algorithms addressing single-cell challenges (cf. Lähnemann et al.)
  #' "real_datasets" is a vector of real scRNA-seq datasets labels; the datasets are sorted by size.
  #' "baron_humpan" is a data.frame with the performances of RSEC and scEFSC on the merged Baron dataset.
  #' "ensemble_benchmark" is a data.frame with the performances of other ensemble algorithms.
  #' "colormap" associates each method to a specific color.
  #' "bordermap" indicates if an ensemble method should be bordered or not.
  #' They are reported from the scEFSC paper (Bian et al. 2022)
  #' 
  algorithms <- c("densityCut", "monocle3", "Seurat", "SHARP", "scEVE")
  ensemble_algorithms <- c("RSEC", "scEVE", "EC.PGMGR", "GRACE", "SAFE", "SAME", "scEFSC")
  challengers <- c("RSEC", "scEVE")
  
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
    "Gillen_HumEPN",       # 18,456c
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
    rep(NA, 4), NA, NA, 0.98, 0.87,
    # RSEC
    rep(NA, 4), NA, NA, 0.80, 0.61,
    # GRACE
    rep(NA, 8)
  )
  ensemble_datasets <- c(
    "Baron_HumPan_1", "Baron_HumPan_2", "Baron_HumPan_3", "Baron_HumPan_4",
    "Baron_MouPan_1", "Baron_MouPan_2", "Li_HumCRC_a", "Tasic_MouBra")
  ensemble_algorithms.prior <- c("SAFE", "SAME", "EC.PGMGR", "scEFSC", "RSEC", "GRACE")
  ensemble_performances <- expand.grid(list(dataset=ensemble_datasets,
                                            method=ensemble_algorithms.prior))
  for (metric in c("ARI", "NMI")) {
    ensemble_performances[, metric] <- get(glue("ensemble_performances.{metric}"))}
  
  baron_humpan <- expand.grid(list(metric=c("ARI", "NMI"), method=c("scEFSC", "RSEC")))
  baron_humpan$y <- c(0.54, 0.75, 0.15, 0.51)
  for (col in c("ymin", "ymax")) {baron_humpan[, col] <- baron_humpan$median}
  
  colormap <- list(densityCut="#e31a1c",
                   monocle3="#33a02c",
                   Seurat="#b15928",
                   SHARP="#ff7f00",
                   scEVE="#1f78b4",
                   RSEC="#e6ab02",
                   EC.PGMGR="#7570b3",
                   GRACE="#666666",
                   SAFE="#1b9e77",
                   SAME="#d95f02",
                   scEFSC="#a6761d")
  
  bordermap <- c(rep("black", 2), rep("#00000000", 5))
  names(bordermap) <- ensemble_algorithms
  
  prior <- list(algorithms=algorithms,
                ensemble_algorithms=ensemble_algorithms,
                challengers=challengers,
                real_datasets=real_datasets,
                colormap=colormap,
                bordermap=bordermap,
                baron_humpan=baron_humpan,
                ensemble_performances=ensemble_performances)
  return(prior)
}

get_results <- function(path) {
  #' Merge the similarity table of each scRNA-seq dataset in a single data.frame.
  #' 
  #' @param path: a character. The path where the tables are stored.
  #' 
  #' @return a data.frame with five columns: 'method_1', 'method_2', 'dataset',
  #' 'ARI' and 'NMI'.
  #' 
  individual_filenames <- list.files(path, full.names=TRUE)
  individual_tables <- lapply(individual_filenames, read.csv)
  results <- do.call(rbind, individual_tables)
  return(results)
}

get_results.benchmark <- function(path="./results/benchmark") {
  #' Merge the benchmark of each scRNA-seq dataset in a single dataframe.
  #' 
  #' @param path: a character. The path where benchmark files are stored.
  #' 
  #' @return a data.frame with ten columns: 'peakRAM', 'time', 'ARI', 'NMI',
  #' 'method', 'dataset', 'real', 'log10(s)', 'log10(Mb)' and 'method_type'.
  #' 
  benchmark <- get_results(path)
  
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

parse_dataset.synthetic <- function(dataset) {
  #' Parse a dataset label to get a 1-row data.frame. The data.frame has 3 columns:
  #' 'populations' (a character), 'related' (a boolean) and 'balanced' (a boolean).
  #' 
  #' @param dataset: a character.
  #' 
  #' @return a data.frame with 3 columns: 'populations', 'related' and 'balanced'.
  #' 
  elements <- strsplit(dataset, split="_")[[1]]
  n_populations <- substr(elements[1], start=2, stop=nchar(elements[1]))
  balanced <- substr(elements[2], 2, 2)
  related <- substr(elements[3], 2, 2)
  row <- data.frame(n_populations=as.character(n_populations),
                    balanced=ifelse(as.logical(balanced), "yes", "no"),
                    related=ifelse(as.logical(related), "yes", "no"))
  return(row)
}

setup_data.synthetic <- function(data.synthetic) {
  #' Get a slim data.frame usable by get_plot.synthetic.
  #' 
  #' @param benchmark.synthetic: a data.frame with a column 'dataset' corresponding
  #' to synthetic datasets.
  #' 
  #' @return a data.frame with 3 extra columns: 'populations', 'related', 'balanced'.
  #' 
  metadata <- lapply(X=data.synthetic$dataset, FUN=parse_dataset.synthetic)
  metadata <- do.call(rbind, metadata)
  data.synthetic <- cbind(data.synthetic, metadata)
  return(data.synthetic)
}

#___________________________________________________________________________real
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

get_boxplot.methods <- function(benchmark.real, metric, legend=TRUE) {
  #' Get a boxplot where the x-axis is methods and the y-axis is a performance metric.
  #' 
  #' @param benchmark.real: a data.frame with three columns: 'dataset', 'method' and a metric.
  #' @param metric: a character. The metric of interest.
  #' @param legend: a boolean. If TRUE, the legend of the boxplot is plotted.
  #' 
  #' @return a plot.
  #' 
  plot <- ggplot(benchmark.real, aes(x=method, y=.data[[metric]])) +
    geom_boxplot(aes(fill=method)) +
    geom_point()
  
  plot <- plot + scale_fill_manual(values=get_prior()$colormap) +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
          legend.position="left", panel.grid.major=element_line(linewidth=0.5),
          axis.line.x=element_blank(), axis.ticks.x=element_blank(),
          legend.title=element_blank(), strip.text=element_blank()) +
    facet_grid(~method_type, scales = "free_x", space = "free_x")
  
  if (!legend) {plot <- plot + theme(legend.position="none")}
  return(plot)
}

get_plot.real <- function(benchmark.real, metric, legend=TRUE) {
  #' Generate a composite plot of the benchmark on real scRNA-seq datasets.
  #' It is composed of a heatmap and a boxplot w.r.t. the method used.
  #' 
  #' @param benchmark.real: a data.frame with three columns: 'dataset', 'method' and a metric.
  #' @param metric: a character. The metric of interest.
  #' @param legend: a boolean. If TRUE, the legend of the boxplot is drawn.
  #' 
  #' @return a plot.
  #' 
  benchmark.real$dataset <- factor(benchmark.real$dataset, levels=get_prior()$real_datasets)
  heatmap.real <- get_heatmap.real(benchmark.real, metric)
  boxplot.methods <- get_boxplot.methods(benchmark.real, metric, legend)
  plot.real <- ggarrange(boxplot.methods, heatmap.real, nrow=2, ncol=1, widths=1,
                         heights=c(3,7), byrow=TRUE)
  return(plot.real)
}

get_plot.real.summary <- function(benchmark.real) {
  #' Generate a composite plot of the benchmark on real scRNA-seq datasets.
  #' It includes plots for the 'NMI', the computation time and the memory usage.
  #' 
  #' @param benchmark.real: a data.frame with five columns: 'dataset', 'method', 'NMI',
  #' 'log10(s)' and 'log10(Mb)'.
  #' 
  #' @return a plot.
  #' 
  metrics <- c("NMI", "log10(s)", "log10(Mb)")
  get_plot.metric <- function(metric) {get_plot.real(benchmark.real, metric, legend=FALSE)}
  plots.metrics <- lapply(X=metrics, get_plot.metric)
  names(plots.metrics) <- metrics
  
  get_legend <- function() {
    plot <- get_boxplot.methods(benchmark.real, "NMI") + 
      theme(legend.position="top", legend.margin=margin(0, 0, 0, 100),
            legend.text=element_text(size=10))
    legend <- ggpubr::get_legend(plot)
    legend <- ggpubr::as_ggplot(legend)}
  
  blank <- ggplot() + theme_void()
  plots <- ggpubr::ggarrange(plots.metrics[["NMI"]], blank, plots.metrics[["log10(s)"]],
                             blank, plots.metrics[["log10(Mb)"]], nrow=1, ncol=5,
                             widths=c(4,.1,4,.1,4))
  plot.real <- ggpubr::ggarrange(get_legend(), plots, nrow=2, ncol=1,
                                 heights=c(1,10))
  return(plot.real)
}

#_______________________________________________________________________ensemble
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
  benchmark.ensemble$challenger <- "no"
  benchmark.ensemble$challenger[benchmark.ensemble$method %in% prior$challengers] <- "yes"
  benchmark.ensemble$challenger <- factor(benchmark.ensemble$challenger, levels=c("yes", "no"))
  return(benchmark.ensemble)
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
  data <- benchmark.ensemble[(benchmark.ensemble$dataset==dataset) &
                               (!is.na(benchmark.ensemble[, metric])),]
  prior <- get_prior()
  
  plot <- ggplot(data=data) +
    geom_col(aes(x=method, y=.data[[metric]], fill=method, color=method), linewidth=1) +
    scale_y_continuous(expand=expansion(mult=0), limits=c(0,1)) +
    scale_fill_manual(values=prior$colormap) +
    scale_colour_manual(values=prior$bordermap) +
    ggtitle(dataset)
  
  plot <- plot +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, vjust=0.8, hjust=0.8),
          panel.grid.major=element_line(linewidth=0.5),
          strip.text=element_blank(), legend.position="none",
          plot.title=element_text(hjust=0.5, size=12)) +
    facet_grid(~challenger, scales = "free_x", space = "free_x")
  
  if (zoom_aes) {plot <- plot +
    theme(panel.background=element_rect(fill="#ebebeb"),
          panel.grid.major=element_line(colour="white"))}
  
  return(plot)
}

setup_baron_humpan.data <- function(benchmark.ensemble, metric) {
  #' Set-up the data to plot the Baron_HumPan performances. Baron_HumPan is split
  #' into four datasets. In the Baron_HumPan plot, the bar indicates the median
  #' performance, and the errorbar indicates the minimum and maximum performances.
  #' 
  #' @param benchmark.ensemble: a data.frame with four columns: 'dataset', 'method', 'ARI' and 'NMI'.
  #' @param metric: a character. One of 'ARI' or 'NMI'.
  #' 
  #' @return a data.frame with five columns: 'method', 'ymin', 'ymax', 'challenger' and a metric.
  #' 
  is_baron_humpan <- grepl("Baron_HumPan", benchmark.ensemble$dataset)
  baron_humpan.data <- benchmark.ensemble[is_baron_humpan, c("method", "challenger", metric)]
  
  get_row.method <- function(method) {
    tmp <- baron_humpan.data[baron_humpan.data$method==method,]
    row.method <- list(method=method, ymin=min(tmp[, metric]), ymax=max(tmp[, metric]),
                       challenger=tmp[1, "challenger"])
    row.method[[metric]] <- median(tmp[, metric])
    row <- data.frame(row.method)
    return(row)}
  
  baron_humpan.data <- lapply(X=levels(baron_humpan.data$method), FUN=get_row.method)
  baron_humpan.data <- do.call(rbind, baron_humpan.data)
  
  # set up ymin and ymax for RSEC and scEFSC
  baron_humpan.prior <- get_prior()$baron_humpan
  for (method in c("scEFSC", "RSEC")) {
    y <- baron_humpan.prior$y[(baron_humpan.prior$method==method) &
                                (baron_humpan.prior$metric==metric)]
    baron_humpan.data[baron_humpan.data$method==method, metric] <- y
  }
  
  baron_humpan.data$method <- factor(baron_humpan.data$method,
                                     levels=get_prior()$ensemble_algorithms)
  return(baron_humpan.data)
}

get_barplot.baron_humpan <- function(baron_humpan.data, metric) {
  #' Get a barplot for the Baron_HumPan dataset, where the bar indicates the
  #' median performance, and the errorbar indicates the minimum and maximum performances.
  #' 
  #' @param baron_humpan.data: a data.frame with five columns: 'method', 'ymin', 'ymax',
  #' 'challenger' and a metric.
  #' @param metric: a character. One of 'ARI' or 'NMI'.
  #' 
  #' @return a plot.
  #' 
  data <- baron_humpan.data[!is.na(baron_humpan.data[, metric]), ]
  prior <- get_prior()
  
  plot <- ggplot(data=data, aes(x=method, y=.data[[metric]], fill=method)) +
    geom_col(aes(color=method), linewidth=1) +
    geom_errorbar(aes(ymin=ymin, ymax=ymax), position=position_dodge(.9), linewidth=1, width=0.1) +
    scale_y_continuous(expand=expansion(mult=0), limits=c(0,1)) +
    scale_fill_manual(values=prior$colormap) +
    scale_colour_manual(values=prior$bordermap) +
    ggtitle("Baron_HumPan")
  
  plot <- plot +
    theme_classic() +
    guides(fill=guide_legend(nrow=1, title.hjust=0)) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=30, vjust=0.8, hjust=0.8),
          panel.grid.major=element_line(colour="white", linewidth=0.5),
          strip.text=element_blank(), panel.background=element_rect(fill="#ebebeb"),
          plot.title=element_text(hjust=0.5)) +
    facet_grid(~challenger, scales = "free_x", space = "free_x")
  
  return(plot)
}

get_plot.ensemble <- function(benchmark.real, metric) {
  #' Get a composite plot representing the performances of multiple ensemble
  #' algorithms on 8 scRNA-seq datasets.
  #' 
  #' @param benchmark.real: a data.frame with four columns: 'dataset', 'method', 'ARI' and 'NMI'.
  #' @param metric: a character. One of 'ARI' or 'NMI'.
  #' 
  #' @return a composite plot.
  #' 
  benchmark.ensemble <- setup_benchmark.ensemble(benchmark.real)
  get_barplot.ensemble.wrapper <- function(dataset) {
    zoom <- ifelse(grepl("Baron_HumPan", dataset), TRUE, FALSE)
    barplot.ensemble <- get_barplot.ensemble(benchmark.ensemble, dataset, metric, zoom)
  }
  
  datasets <- unique(benchmark.ensemble$dataset)
  barplots <- lapply(X=datasets, FUN=get_barplot.ensemble.wrapper)
  names(barplots) <- datasets
  
  baron_humpan.data <- setup_baron_humpan.data(benchmark.ensemble, metric)
  barplot.baron_humpan <- get_barplot.baron_humpan(baron_humpan.data, metric) + 
    labs(fill="", color="")
  
  barplots.baron_humpan.subfigures <- ggpubr::ggarrange(
    barplots$Baron_HumPan_1, barplots$Baron_HumPan_2,
    barplots$Baron_HumPan_3, barplots$Baron_HumPan_4,
    nrow=2, ncol=2)
  
  barplots.upper <- ggpubr::ggarrange(barplots$Baron_MouPan_1, barplots$Baron_MouPan_2,
                                      barplots$Li_HumCRC_a, barplots$Tasic_MouBra, align="h",
                                      nrow=1, ncol=4)
  
  barplots.lower <- ggpubr::ggarrange(barplot.baron_humpan, barplots.baron_humpan.subfigures,
                                      nrow=1, ncol=2, common.legend=TRUE, legend="bottom")
  
  plot.ensemble <- ggpubr::ggarrange(barplots.upper, barplots.lower, nrow=2, heights=c(1,2))
  return(plot.ensemble)
}

#__________________________________________________________________________trees
get_resolution_tree.data <- function(meta) {
  #' Get a data.frame set up to generate a tree representing the 
  #' multi-resolution clustering of scEVE.
  #' 
  #' @param meta: a data.frame with two columns: 'parent', 'consensus' and 'n'.
  #' 
  #' @return a phylo tree object.
  #' 
  data <- meta[-1, ]
  input <- data.frame(parent=as.character(data$parent), node=rownames(data))
  resolution_tree <- as.phylo(input)
  resolution_tree.coordinates <- fortify(resolution_tree)
  resolution_tree.data <- as_tibble(resolution_tree)
  
  for (col in c("consensus", "n")) {meta[, col] <- as.numeric(meta[, col])}
  rows_order <- rownames(meta)[order(meta$consensus, -meta$n)]
  # sort by consensus and ensure that C is the first row.
  resolution_tree.data <- resolution_tree.data[match(rows_order, resolution_tree.data$label),]
  resolution_tree.coordinates <- resolution_tree.coordinates[match(rows_order, resolution_tree.coordinates$label),]
  meta <- meta[rows_order,]
  
  for (col in c("consensus", "n")) {resolution_tree.data[, col] <- meta[, col]}
  for (col in c("isTip", "x", "y", "branch", "angle")) {
    resolution_tree.data[, col] <- resolution_tree.coordinates[, col]}
  resolution_tree.data$node_type <- "consensus.cluster"
  resolution_tree.data$node_type[resolution_tree.data$consensus==0] <- "root.leftover"
  return(resolution_tree.data)
}

get_resolution_tree <- function(resolution_tree.data) {
  #' Draw a classification tree summarizing the analysis of scEVE.
  #' 
  #' @param resolution_tree.data: a tibble object with 10 columns: 'parent', 'node',
  #' 'label', 'consensus', 'isTip', 'x', 'y', 'branch', 'angle' and 'node_type'.
  #' 
  #' @return a ggtree plot.
  #' 
  colors <- list(root.leftover="white", consensus.cluster="black")
  get_color <- function(node_type) {colors[[node_type]]}
  resolution_tree.data$colors <- sapply(X=resolution_tree.data$node_type, FUN=get_color)
  
  plot <- ggtree(resolution_tree.data) +
    geom_label(data=resolution_tree.data[-1,], aes(x=branch, label=round(consensus, 2)),
               label.size=NA, hjust="right") +
    geom_label(aes(label=label, fill=node_type, colour=node_type),
               hjust=ifelse(resolution_tree.data$isTip, "right", "middle")) +
    scale_colour_manual(values=resolution_tree.data$colors) +
    scale_fill_manual(values=list(root.leftover="#36393d", consensus.cluster="#e6e6e6")) +
    guides(fill="none", color="none")
  return(plot)
}

get_barplot.tree.heterogeneity <- function(distributions.data) {
  #' Get a horizontal barplot associating its leaf cluster to its cell type composition.
  #' 
  #' @param distributions.data: a data.frame with three columns: 'n', 'ground_truth' and 'population'.
  #'
  #' @return a plot.
  #' 
  cell_types <- unique(distributions.data$ground_truth)
  nrows <- length(cell_types) %/% 4
  
  plot <- ggplot(data=distributions.data) +
    geom_bar(aes(x=population, y=n, fill=ground_truth),
             position="fill", stat="identity", color="black") +
    theme_classic() + 
    scale_y_continuous(expand=expansion(mult=0), labels=percent) +
    guides(fill=guide_legend(nrow=nrows, title="cell types")) +
    coord_flip() + scale_fill_brewer(palette="Dark2") +
    ylab("cell types")
  return(plot)
}

get_barplot.tree.counts <- function(distributions.data) {
  #' Get a horizontal barplot associating its leaf cluster to its number of cells.
  #' 
  #' @param distributions.data: a data.frame with three columns: 'n', 'ground_truth' and 'population'.
  #'
  #' @return a plot.
  #' 
  data <- distributions.data %>%
    group_by(population) %>%
    summarise(n=sum(n))
  
  plot <- ggplot(data=data) +
    geom_bar(aes(x=population, y=n), stat="identity", fill="black") +
    theme_classic() +
    scale_y_log10(expand=expansion(mult=0), guide="axis_logticks") +
    coord_flip() + scale_fill_brewer(palette="Paired") +
    theme(panel.grid.major.x=element_line(linewidth=0.5)) +
    ylab("# cells")
  return(plot)
}

get_plot.tree <- function(resolution_tree.data, distributions.data, dataset) {
  #' Get a composite plot aligning a resolution tree to barplots representing
  #' the proportion of cell types and the total count of cells in leaf populations.
  #' 
  #' @param resolution_tree.data: a tibble object with 10 columns: 'parent', 'node',
  #' 'label', 'consensus', 'isTip', 'x', 'y', 'branch', 'angle' and 'node_type'.
  #' @param distributions.data: a data.frame with three columns: 'proportion',
  #' 'ground_truth' and 'population'.
  #' @param dataset: a character.
  #' 
  #' @return a composite plot.
  #' 
  #'
  leaves <- resolution_tree.data$label[resolution_tree.data$isTip]
  distributions.data <- distributions.data[distributions.data$population %in% leaves, ]
  barplot.heterogeneity <- get_barplot.tree.heterogeneity(distributions.data)
  barplot.counts <- get_barplot.tree.counts(distributions.data)
  resolution_tree <- get_resolution_tree(resolution_tree.data)
  
  resolution_tree <- resolution_tree + theme(plot.margin=unit(c(0, 0, 0, 0), "null"))
  
  apply_aes <- function(barplot) {barplot +
      theme(axis.line.y=element_blank(), axis.title.y=element_blank(),
            axis.text.y=element_blank())}
  
  plot.tree <- resolution_tree +
    apply_aes(barplot.heterogeneity) +
    apply_aes(barplot.counts) +
    plot_layout(widths=c(8,3,3), guides="collect") &
    theme(legend.position="bottom")
  
  return(plot.tree)
}

#_____________________________________________________________________signatures
setup_signatures.data.population <- function(population, markers, cancer_signatures) {
  #' Get a data.frame associating the markers of a cell population to cancer signatures.
  #' 
  #' @param population: a character.
  #' @param markers: a data.frame where rows are markers and cols are cell populations.
  #' @param cancer_signatures: a data.frame with two columns: 'gene' and 'signature'.
  #' 
  #' @return a data.frame with three columns: 'population', 'signature' and 'n_markers'.
  #' 
  markers.population <- rownames(markers[markers[, population]==1, ])
  signatures <- unique(cancer_signatures$signature)
  
  get_n_markers.signature <- function(signature) {
    markers.signature <- cancer_signatures[cancer_signatures$signature==signature, "gene"]
    data <- intersect(markers.population, markers.signature)
    n_markers <- length(data)
    return(n_markers)
  }
  
  signatures.data.population <- lapply(X=signatures, FUN=get_n_markers.signature)
  signatures.data.population <- do.call(rbind, signatures.data.population)
  signatures.data.population <- as.data.frame(signatures.data.population)
  colnames(signatures.data.population) <- c("n_markers")
  
  signatures.data.population$signature <- signatures
  signatures.data.population$cell_population <- population
  colnames(signatures.data.population) <- c("n_markers", "signature", "population")
  return(signatures.data.population)
}

setup_signatures.data <- function(markers, cancer_signatures) {
  #' Get a data.frame set up to generate a barplot representing the cancer signatures.
  #' 
  #' @param markers: a data.frame where rows are markers and cols are cell populations.
  #' @param cancer_signatures: a data.frame with two columns: 'gene' and 'signature'.
  #' 
  #' @return a data.frame with three columns: 'population', 'signature' and 'n_markers'.
  #' 
  n_markers.signature <- table(cancer_signatures$signature)
  get_label.signature <- function(signature) {
    glue("{signature} ({n_markers.signature[[signature]]})")} 
  labels.signature <- sapply(X=cancer_signatures$signature, FUN=get_label.signature)
  cancer_signatures$signature <- labels.signature
  
  populations <- colnames(markers)
  signatures.data.population <- lapply(X=populations, FUN=setup_signatures.data.population,
                                       markers=markers, cancer_signatures=cancer_signatures)
  signatures.data <- do.call(rbind, signatures.data.population)
  return(signatures.data)
}

get_plot.signatures <- function(signatures.data) {
  #' Get a barplot representing the number of cancer markers w.r.t. the cell
  #' population and the cancer signature.
  #'
  #' @param signatures.data: a data.frame with three columns: 'population', 'signature'
  #' and 'n_markers'.
  #' 
  #' @return a plot.
  #' 
  data <- signatures.data[signatures.data$population %in% c("C.5.1", "C.5.3"), ]
  plot <- ggplot(data=data) +
    geom_bar(aes(x=signature, y=n_markers, fill=population), color="black", 
             stat="identity", position="dodge", width=0.75) +
    theme_classic() +
    scale_y_continuous(expand=expansion(mult=0)) +
    coord_flip() + scale_fill_manual(values=list(C.5.1="#ebebeb", C.5.3="#1F78B4")) +
    theme(panel.grid.major.x=element_line(linewidth=0.5),
          legend.title=element_blank(), axis.title.y=element_blank(),
          legend.position="top", legend.key.height=unit(0.4, "line")) + 
    ylab("# markers")
  return(plot)
}

#______________________________________________________________________synthetic
get_plot.synthetic.method <- function(benchmark.synthetic, metric, method) {
  #' Get a composite plot of the performances of a method under different conditions
  #' with synthetic scRNA-seq datasets.
  #' 
  #' @param benchmark.synthetic: a data.frame with five columns: 'balanced',
  #' 'related', 'populations', 'method' and a metric.
  #' @param metric: a character.
  #' @param method: a character. 
  #' 
  #' @return a plot.
  #' 
  data <- benchmark.synthetic[benchmark.synthetic$method==method, ]
  
  related.labels <- c(yes="related", no="independent")
  balanced.labels <- c(yes="balanced", no="imbalanced")
  color <- get_prior()$colormap[[method]]
  
  plot <- ggplot(data=data) +
    geom_boxplot(aes(x=n_populations, y=.data[[metric]]),
                 fill="white", color=color) +
    scale_y_continuous(limits=c(0, 1)) +
    facet_grid(related ~ balanced, #switch="y",
               labeller=labeller(related=related.labels, balanced=balanced.labels))
  
  plot <- plot + theme_classic() + xlab("# clusters") +
    theme(panel.background=element_rect(fill="#ebebeb"),
          panel.grid.major=element_line(colour="white", linewidth=0.5),
          panel.border=element_rect(colour="black", fill=NA, linewidth=1),
          axis.line=element_blank(), strip.text=element_text(color="white"),
          strip.background=element_rect(color=NA, fill=color))
  
  top_facet <- ggpubr::text_grob("cluster sizes", hjust=0.3)
  plot <- grid.arrange(plot, top=top_facet, right="cluster transcriptomes")
  return(plot)
}

get_boxplot.synthetic <- function(benchmark.synthetic, metric, method, related, balanced) {
  #' Get a boxplot showing the performances of a given method on synthetic datasets
  #' that can be balanced or related.
  #' 
  #' @param benchmark.synthetic: a data.frame with five columns: 'balanced',
  #' 'related', 'populations', 'method' and a metric.
  #' @param metric: a character. The metric of interest.
  #' @param method: a character. The method of interest.
  #' @param related: a character. One of 'yes' or 'no'.
  #' @param balanced: a character. One of 'yes' or 'no'.
  #' 
  #' @return a plot.
  #' 
  data <- benchmark.synthetic[(benchmark.synthetic$method==method) &
                                (benchmark.synthetic$related==related) &
                                (benchmark.synthetic$balanced==balanced),]
  
  plot <- ggplot(data=data) +
    geom_boxplot(aes(x=n_populations, y=.data[[metric]]),
                 fill="white", color=get_prior()$colormap[[method]]) +
    scale_y_continuous(limits=c(0, 1)) + xlab("# clusters")
  
  plot <- plot + facet_grid(. ~ method) +
    theme_classic() + 
    theme(panel.background=element_rect(fill="#ebebeb"),
          panel.grid.major=element_line(colour="white", linewidth=0.5),
          strip.background=element_rect(fill=get_prior()$colormap[[method]], color=NA),
          strip.text=element_text(face="bold", colour="white")
    )
  
  return(plot)
}

get_plot.synthetic.setting <- function(benchmark.synthetic, metric, related, balanced) {
  #' Get a composite plot showing the performances of scEVE on a specific data setting
  #' w.r.t. the performances of its individual methods.
  #' 
  #' @param benchmark.synthetic: a data.frame with four columns: 'balanced',
  #' 'related', 'populations' and a metric.
  #' @param metric: a character. The metric of interest.
  #' @param related: a character. One of 'yes' or 'no'.
  #' @param balanced: a character. One of 'yes' or 'no'.
  #' 
  #' @return a plot.
  #' 
  get_boxplot.synthetic.wrapper <- function(method) {
    get_boxplot.synthetic(benchmark.synthetic, metric, method, related, balanced)}
  methods <- unique(benchmark.synthetic$method)
  boxplots <- lapply(X=methods, FUN=get_boxplot.synthetic.wrapper)
  names(boxplots) <- methods

  blank <- ggplot() + theme_void()
  boxplots.individual <- ggpubr::ggarrange(boxplots$densityCut, boxplots$monocle3,
                                           blank, blank,
                                           boxplots$Seurat, boxplots$SHARP,
                                           nrow=3, ncol=2, heights=c(5,0.1,5))
  boxplot.scEVE <- ggpubr::ggarrange(blank, boxplots$scEVE, blank, nrow=1, ncol=3,
                                     widths=c(.5,1,.5))
  plot.synthetic.setting <- ggpubr::ggarrange(boxplot.scEVE, boxplots.individual,
                                              nrow=2, ncol=1, heights=c(1,2))
  return(plot.synthetic.setting)
}

get_plot.synthetic.summary <- function(benchmark.synthetic, metric) {
  #' Get a composite plot showing the performances of all the methods
  #' on the synthetic datasets.
  #' 
  #' @param benchmark.synthetic: a data.frame with four columns: 'balanced',
  #' 'related', 'populations' and a metric.
  #' @param metric: a character. The metric of interest.
  #' 
  #' @return a plot.
  #' 
  get_plot.synthetic.method.wrapper <- function(method) {
    plot <- get_plot.synthetic.method(benchmark.synthetic, metric, method)
  }
  
  methods <- unique(benchmark.synthetic$method)
  plots <- lapply(X=methods, FUN=get_plot.synthetic.method.wrapper)
  names(plots) <- methods
  
  blank <- ggplot() + theme_void()
  plot.scEVE <- ggpubr::ggarrange(blank, plots$scEVE, blank, nrow=1,
                                  ncol=3, widths=c(0.5, 1, 0.5))
  plots.methods <- ggpubr::ggarrange(plots$densityCut, blank, plots$monocle3,
                                     blank, blank, blank,
                                     plots$Seurat, blank, plots$SHARP, nrow=3, ncol=3,
                                     widths=c(5,0.5,5), heights=c(5,0.1,5))
  plot.synthetic.summary <- grid.arrange(plots.methods, plot.scEVE, nrow=2, ncol=1,
                                         heights=c(2,1))
  
  return(plot.synthetic.summary)
}
