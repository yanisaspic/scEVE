"Library with functions to draw patchwork figures in R.

	2023/12/04 @yanisaspic"

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(Seurat)
  library(SCpubr)
  library(tibble)
  library(ggpubr)
  library(reshape)
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
})


get_existing_pdfs <- function(population) {
  #' Get the names of the existing pdf files wrt the population of interest.
  #'
  #' @param population: the name of a population.
  #'
  #' @return a vector of filenames.
  #'
  files <- c()
  for (cat in c("trim", "cluster", "draw")) {
    filename <- glue("./results/figures/{POPULATION}_{cat}.pdf")
    if (file.exists(filename)) {
      files <- c(files, filename)
    }
  }
  return(files)
}


is_child_of <- function(group_label, population) {
  #' Test if a label corresponds to a child of the population label.
  #' 
  #' @param group_label: a seed or a cluster label.
  #' @param parent: a population.
  #' 
  #' @return a boolean.
  #' 
  is_next_resolution <- nchar(group_label) == nchar(population) + 1
  is_related <- substr(group_label, 1, nchar(population)) == population
  return(is_next_resolution & is_related)
}


get_ground_truth <- function(cell_name) {
  #' Get the expected group of a cell.
  #' 
  #' @param cell_name: a cell name of the form '{ground_truth}_X'
  #' 
  #' @return a char.
  #' 
  return(strsplit(x = cell_name, split="_")[[1]][1])
}


get_classification <- function(cells) {
  #' Get the prediction and the ground truth of a subset of cells.
  #' 
  #' @param cells: a binary data.frame associating cells (rows) to their group (cols).
  #' 
  #' @return a data.frame with two columns: 'pred' and 'GT'.
  #' 
  data <- melt(as.matrix(cells), varnames=c("cell_name", "pred"))
  data[, "GT"] <- sapply(as.character(data$cell_name), get_ground_truth)
  data <- data[data$value == 1, ]
  rownames(data) <- data$cell_name
  
  data$pred <- as.character(data$pred)
  data$GT <- as.character(data$GT)
  return(data[, c("pred", "GT")])
}


get_composition_barplot <- function(classification) {
  #' Get a composition barplot.
  #' 
  #' @param classification: a data.frame of 2 columns: 'pred' and 'GT'
  #' 
  #' @return a barplot with x=pred, fill=GT and y=n_cells (i.e. number of cells in pred and GT).
  #' 
  data <- classification %>% group_by(pred, GT) %>% summarise(n_cells=n())
  data_wo.gt <- classification %>% group_by(pred) %>% summarise(n_cells=n())
  
  # functional parameters
  g <- ggplot(data=data, aes(fill=GT, y=n_cells, x=pred)) + 
    geom_bar(position="stack", stat="identity") +
    coord_flip()
  
  # aesthetic parameters
  g <- g +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor= element_blank(),
      panel.border = element_blank()
    ) +
    theme(axis.line = element_line(color = 'black')) +
    scale_y_continuous(
      expand = c(0,0), 
      limits = c(0, sum(data$n_cells)), 
      breaks = c(0, sum(data$n_cells)),
    ) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()
    )
  
  return(g)
}


get_specificity.markers <- function(expr.markers_inout) {
  #' Quantify the specificity of a set of markers with regard to a group.
  #' 
  #' @param expr.markers_inout: a list with two data.frames of scRNA-seq, named 'inside' and 'outside': genes (rows) and cells (cols)
  #' 
  #' @return a data.frame with two columns: 'marker' and 'specificity'.
  #' 
  get_binary <- function(v){as.logical(v)}
  get_representation <- function(v){(sum(v)+1)/(length(v)+1)}
  get_representation.gene <- function(expr.gene){get_representation(get_binary(expr.gene))}
  get_representation.markers <- function(expr.markers){apply(X=expr.markers, MARGIN = 1, FUN = get_representation.gene)}
  representation.markers <- lapply(expr.markers_inout, FUN=get_representation.markers)
  
  specificity.markers <- representation.markers$inside / representation.markers$outside
  specificity.markers <- enframe(specificity.markers, name="marker", value="specificity")
  return(specificity.markers)
}
  

split_cells_inout <- function(classification, group_label) {
  #' Split a cell population into an 'inside' and an 'outside' group.
  #'
  #' @param classification: a data.frame with two columns: 'pred' and 'GT'.
  #' @param group_label: the name of a group.
  #' 
  #' @return a list with two vectors of cell names, named 'inside' and 'outside'.
  #' 
  is_in_group <- classification$pred==group_label
  cells <- list(
    inside=rownames(classification[is_in_group,]),
    outside=rownames(classification[!is_in_group,])
  )
  return(cells)
}


get_markers_of_group <- function(markers, group_label) {
  #' Get a vector of markers associated to a group.
  #' 
  #' @param markers: a binary data.frame associating markers (rows) to their groups (cols).
  #' @param group_label: the name of a group.
  #' 
  #' @return a vector of gene names
  #' 
  data <- markers %>% 
    select(all_of(group_label)) %>%
    filter(markers[group_label]==1)
  return(rownames(data))
}


get_top_threshold <- function(v) {
  #' Get the threshold value to split top individuals from their distribution.
  #' It corresponds to Q3 + 1.5 * IQR.
  #' 
  #' @param v: a vector of numerics
  #' 
  #' @return a numeric
  #'
  parameters <- summary(v)
  q1 <- parameters[["1st Qu."]]
  q3 <- parameters[["3rd Qu."]]
  iqr <- q3 - q1
  threshold <- q3 + 1.5 * iqr
  return(threshold)
}


add_top_column <- function(specificity.group) {
  #' Test if a marker is very specific. True if a marker is an outlier or the single most specific.
  #' 
  #' @param specificity.group: a data.frame with three columns: 'marker', 'pred', and 'specificity'.
  #' 
  #' @return a data.frame with four columns: 'marker', 'pred', 'specificity', 'top'.
  #' 
  max_specificity <- max(specificity.group$specificity)
  top_specificity_threshold <- get_top_threshold(specificity.group$specificity)
  is_top.gene <- function(x){x >= top_specificity_threshold | x == max_specificity}
  specificity.group[, "top"] <- sapply(X=specificity.group$specificity, FUN=is_top.gene)
  return(specificity.group)
}


get_specificity.group <- function(group_label, markers, classification, expression) {
  #' Quantify the specificity of the marker genes of a single group of cells.
  #' 
  #' @param group_label: the name of a group.
  #' @param markers: a binary data.frame associating markers (rows) to their groups (cols).
  #' @param classification: a data.frame with two columns: 'pred' and 'GT'.
  #' @param expression: a data.frame of scRNA-seq: genes (rows) and cells (cols).
  #'
  #' @return a data.frame with three columns: 'marker', 'pred', and 'specificity'.
  #' 
  cells_inout <- split_cells_inout(classification, group_label)
  markers_of_group <- get_markers_of_group(markers, group_label)
  expr_inout <- list(
    inside=expression[markers_of_group, cells_inout$inside],
    outside=expression[markers_of_group, cells_inout$outside])
  specificity.group <- get_specificity.markers(expr_inout)
  specificity.group[, "pred"] <- group_label
  return(specificity.group)
}


get_specificity.frame <- function(expression, markers, classification) {
  #' Quantify the specificity of the marker genes to their respective group of cells.
  #' 
  #' @param expression: a data.frame of scRNA-seq: genes (rows) and cells (cols).
  #' @param markers: a binary data.frame associating markers (rows) to their groups (cols).
  #' @param classification: a data.frame with two columns: 'pred' and 'GT'.
  #' 
  #' @return a data.frame with four columns: 'marker', 'pred', 'specificity', 'top'.
  #' 
  groups <- as.character(unique(classification$pred))
  specificity.groups <- lapply(
    X=groups, 
    FUN=get_specificity.group, 
    markers=markers, 
    classification=classification, 
    expression=expression)
  specificity.groups <- lapply(X=specificity.groups, FUN=add_top_column)
  specificity.frame <- merge_all(specificity.groups)
  specificity.frame <- specificity.frame %>%
    arrange(pred, desc(specificity))
  
  specificity.frame$pred <- as.character(specificity.frame$pred)
  specificity.frame$marker <- as.character(specificity.frame$marker)
  specificity.frame$specificity <- as.numeric(specificity.frame$specificity)
  specificity.frame$top <- as.logical(specificity.frame$top)
  return(specificity.frame)
}


get_widths <- function(specificity.frame) {
  #' Get the widths of the bars in a markers barplot.
  #' They are relative to the number of specific genes.
  #' 
  #' @param specificity.frame: a data.frame with four columns: 'marker', 'pred', 'specificity' and 'top'.
  #' 
  #' @return a data.frame with three columns: 'start', 'middle' and 'end'.
  #' 
  data <- specificity.frame[specificity.frame$top, ] %>%
    arrange(desc(specificity)) %>%
    distinct(marker, .keep_all = TRUE) %>%
    arrange(pred, desc(specificity))
  widths <- data %>%
    group_by(pred) %>%
    summarise(rel=n()/nrow(data)*100)
  
  gap <- 1
  widths[, "end"] <- cumsum(widths$rel) - gap
  widths[, "start"] <- widths$end - widths$rel + gap
  widths[, "middle"] <- widths$start + (widths$end - widths$start) / 2
  widths[widths$start==0, "start"] = gap
  return(widths)
}


get_markers_barplot <- function(specificity.frame) {
  #' Get a markers barplot.
  #' 
  #' @param specificity.frame: a data.frame with four columns: 'marker', 'pred', 'specificity' and 'top'.
  #' 
  #' @return a barplot with x=pred, width=n_top and y=n_markers 
  #' (i.e. number of very specific markers and number of total markers respectively).
  #'
  data <- specificity.frame %>% 
    group_by(pred) %>% 
    summarise(n_markers=n())
  data[, "n_top"] <- specificity.frame[specificity.frame$top, ] %>% 
    group_by(pred) %>% 
    summarise(n_markers=n()) %>% 
    select(n_markers)
  
  widths <- get_widths(specificity.frame)
  data <- inner_join(data, widths, by="pred")

  # functional parameters
  g <- ggplot(data=data, aes(y=n_markers, label=pred)) +
    geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=n_markers), position="identity", stat="identity")
  
  # aesthetic parameters
  g <- g +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor= element_blank(),
      panel.border = element_blank()
    ) +
    theme(axis.line = element_line(color = 'black')) +
    scale_y_continuous(
      expand = expansion(mult=c(0,0.1)), 
      limits = c(0, NA),
      breaks = c(0, data$n_markers)
    ) +
    scale_x_continuous(
      expand = expansion(mult=c(0,0)), 
      limits = c(0, 100)
    ) +
    theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank()) +
    geom_label(aes(x=middle, y=0, vjust=0))
  
  return(g)
}


get_children_consensus <- function(meta) {
  #' Get a data.frame associating groups to their consensus.
  #' 
  #' @param meta: a data.frame with a row 'consensus' and groups (cols).
  #' 
  #' @return a data.frame with two columns: 'pred' and 'consensus'.
  #' 
  data <- meta["consensus",]
  data <- melt(data, id.vars=NULL)
  colnames(data) <- c("pred", "consensus")
  
  data$pred <- as.character(data$pred)
  data$consensus <- as.numeric(data$consensus)
  return(data)
}


get_consensus_barplot <- function(children_consensus, parent_consensus) {
  #' Get a consensus barplot.
  #' 
  #' @param children_consensus: a data.frame with two columns: 'pred' and 'consensus'.
  #' @param parent_consensus: a numeric.
  #'
  #' @return a barplot with x=pred and y=consensus. The red line represents the parent consensus. 
  #' 
  
  # functional parameters
  g <- ggplot(data=children_consensus, aes(y=consensus, x=pred, label=pred)) +
    geom_bar(position="identity", stat="identity") +
    scale_x_discrete(name = "", position = "top") +
    coord_flip()
  
  # aesthetic parameters
  g <- g +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor= element_blank(),
      panel.border = element_blank()
    ) +
    theme(axis.line = element_line(color = 'black')) +
    geom_hline(yintercept = as.numeric(parent_consensus), color="red", linetype=2) +
    scale_y_reverse(
      expand = expansion(mult=c(0.1,0)),
      limits = c(1, 0),
      breaks = c(0, 1)
    ) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    geom_label(aes(hjust=1, y=0))
  

  return(g)
}


get_top_dotplot <- function(expression, classification, specificity.frame) {
  #' Get a dotplot of the top marker genes.
  #'
  #' @param expression: a data.frame of scRNA-seq: genes (rows) and cells (cols).
  #' @param classification: a data.frame with two columns: 'pred' and 'GT'.
  #' @param specificity.frame: a data.frame with four columns: 'marker', 'pred', 'specificity' and 'top'.
  #' 
  #' @return a dotplot with x=pred and y=top_markers.
  #' 
  data <- specificity.frame[specificity.frame$top, ]
  data <- data %>%
    arrange(desc(specificity)) %>%
    distinct(marker, .keep_all = TRUE) %>%
    arrange(pred, desc(specificity))
  specific_markers <- data$marker
  
  # set-up Seurat Object
  so <- CreateSeuratObject(expression)
  preds <- factor(classification$pred)
  names(preds) <- rownames(classification)
  so@active.ident <- preds

  # functional parameters
  g <- do_DotPlot(
    sample = so, 
    features = specific_markers, 
    scale = TRUE,
    dot.scale = 8,
    legend.type = "normal",
    use_viridis = TRUE)
  
  # aesthetic parameters
  g <- g +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  return(g)
}


get_summary_plot <- function(expression, classification, markers, meta, parent_consensus) {
  #' Get a summary plot with: 
  #'  - a consensus barplot, 
  #'  - a markers barplot, 
  #'  - a composition barplot,
  #'  - and a top markers dotplot.
  #'  
  #' @param expression: a data.frame of scRNA-seq: genes (rows) and cells (cols).
  #' @param classification: a data.frame with two columns: 'pred' and 'GT'.
  #' @param markers: a binary data.frame associating markers (rows) to their groups (cols).
  #' @param meta: a data.frame with a row 'consensus' and groups (cols).
  #' @param parent_consensus: a numeric.
  #'
  #' @return a composite plot
  #' 
  consensus_barplot <- get_consensus_barplot(get_children_consensus(meta), parent_consensus)
  composition_barplot <- get_composition_barplot(classification) +
    theme(legend.key.size = unit(.2, "cm")) +
    guides(fill = guide_legend(ncol = 2))
  composition_barplot_legend <- get_legend(composition_barplot)
  params <- composition_barplot$theme
  
  specificity.frame <- get_specificity.frame(expression, markers, classification)
  markers_barplot <- get_markers_barplot(specificity.frame)
  top_dotplot <- get_top_dotplot(expression, classification, specificity.frame)
  top_dotplot_legend <- get_legend(top_dotplot)

  composition_barplot <- composition_barplot +
    theme(plot.margin = margin(0, 3, 15, 0),
          legend.position="none")
  consensus_barplot <- consensus_barplot +
    theme(plot.margin = margin(0, 0, 15, 3))
  markers_barplot <- markers_barplot +
    theme(plot.margin = margin(0, -3, 3, -30))
  top_dotplot <- top_dotplot +
    theme(plot.margin = margin(0, 0, 0, -3),
          legend.position="none",
          text = params$text,
          axis.title = params$axis.title,
          axis.text = params$axis.text)
  
  filler <- ggplot() + theme_void()

  WIDTHS <- c(0.15, 0.65, 0.20)
  upper_g <- ggarrange(
    plotlist=list(filler, markers_barplot, composition_barplot_legend),
    ncol=3,
    widths=WIDTHS)

  lower_g <- ggarrange(
    plotlist=list(consensus_barplot, top_dotplot, composition_barplot),
    ncol=3,
    widths=WIDTHS
  )

  
  g <- ggarrange(plotlist=list(upper_g, lower_g), nrow=2, heights=c(0.2, 0.8))

  g <- g +
    theme(plot.margin = margin(10, 10, 10, 10))

  return(g)
}


draw_seeds <- function(SeurObj, seeds) {
  #' Get a summary plot from records.
  #' 
  #' @param expression: a data.frame of scRNA-seq: genes (rows) and cells (cols).
  #' @param records: a list of three data.frames: 'meta', 'cells', 'markers'
  #' @param parent_consensus: a numeric.
  #' 
  #' @return a composite plot
  #' 
  markers <- records[["markers"]][, order(colnames(records[["markers"]]))]
  cells <- records[["cells"]][, order(colnames(records[["cells"]]))]
  meta <- records[["meta"]][, order(colnames(records[["meta"]]))]

  classification <- get_classification(cells)
  g <- get_summary_plot(expression, classification, markers, meta, parent_consensus)

  return(g)
}