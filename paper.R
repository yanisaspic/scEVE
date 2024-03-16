"Run this script to generate the comparisons required for the scEVE paper.

	2024/02/28 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(SummarizedBenchmark)
})
source("./src/paper/methods.R")
source("./src/paper/metrics.R")

#___________________________________________________________________________init
EXPRESSION.INIT <- read.csv("./data/datasets/Li_HumCRC.csv", header=TRUE, row.names=1)
GROUND_TRUTH <- get_ground_truth(EXPRESSION.INIT)


DATA <- list(expression.init=EXPRESSION.INIT, ground_truth=GROUND_TRUTH, n_HVGs=5000, random_state=0)

bd <- BenchDesign(data=DATA)
bd <- bd %>%
  add.individual_method(bd, "Seurat") %>%
  add.individual_method(bd, "monocle3") %>%
  add.individual_method(bd, "SHARP") %>%
  add.individual_method(bd, "densityCut")

#____________________________________________________________________________run
sb <- buildBench(bd=bd, truthCols=c(labels="ground", time="time", peakRAM="peakRAM"))

# the iterative pre-processing steps of scEVE
# are incompatible with SummarizedBenchmark
# -> we benchmark it separately below
###
params <- get_default_hyperparameters()
benchmark_scEVE.strategy <- function(strategy) {
  params$leftovers_strategy <- strategy
  results <- benchmark_scEVE(expression.init=expression.init,
                             params=params,
                             random_state=RANDOM_STATE)
  return(results)
}
leftovers_strategy <- c("default", "assign", "assign_weighted")
scEVE.results <- lapply(X=leftovers_strategy, FUN=benchmark_scEVE.with_leftovers_strategy)

#############
# add metrics
#############
sb <- sb %>%
  addPerformanceMetric(assay="labels", evalMetric="ARI", evalFunction=get_ARI) %>%
  addPerformanceMetric(assay="labels", evalMetric="AMI", evalFunction=get_AMI) %>%
  addPerformanceMetric(assay="labels", evalMetric="NMI", evalFunction=get_NMI) %>%
  addPerformanceMetric(assay="time", evalMetric="time", evalFunction=get_time) %>%
  addPerformanceMetric(assay="peakRAM", evalMetric="peakRAM", evalFunction=get_peakRAM)
metrics <- estimatePerformanceMetrics(sb, tidy=TRUE)[, c("label", "value", "performanceMetric")]
metrics$highlight <- "no"
metrics <- add_scEVE_metrics(metrics, scEVE_results, ground_truth)

# draw metrics
##############
metrics$label <- as.factor(metrics$label)
metrics$performanceMetric <- as.factor(metrics$performanceMetric)
metrics$value <- as.numeric(metrics$value)

pdf(file="benchmark.pdf")
quality_plot <- ggplot() + 
  geom_col(data=metrics[metrics$performanceMetric %in% c("AMI", "ARI", "NMI"),], 
           aes(x=reorder(label, -value), y=value, fill=label, color=highlight, linetype=highlight)) +
  scale_color_manual(values=c("yes"="black"), guide=FALSE) +
  scale_linetype_manual(values=c("yes"="solid", "no"="blank"), guide=FALSE) +
  facet_grid(performanceMetric ~ ., scale="free_y")
computation_plot <- ggplot() +
  geom_col(data=metrics[metrics$performanceMetric %in% c("time", "peakRAM"),],
           aes(x=reorder(label, -value), y=value, fill=label, color=highlight, linetype=highlight)) +
  scale_color_manual(values=c("yes"="black"), guide=FALSE) +
  scale_linetype_manual(values=c("yes"="solid", "no"="blank"), guide=FALSE) +
  facet_grid(performanceMetric ~ ., scale="free_y") +
  ylab("log10(value)")
print(quality_plot)
print(computation_plot)
dev.off()
