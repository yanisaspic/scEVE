"Run this script to generate the benchmark of scEVE.

	2023/01/20 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(SummarizedBenchmark)
})
source("./scEVE.R")
source("./src/benchmark/metrics.R")
source("./src/benchmark/methods.R")

get_benchmark.scEVE <- function(expression.init, params, random_state) {
  #' Get the results of scEVE by applying it with a given set of parameters on a scRNA-seq raw count matrix.
  #' These results include:
  #' - peakRAM (Mo): the maximum memory usage of the method.
  #' - time (s): the computation time in seconds.
  #' - preds: a named factor, where names are cells and values are cluster labels.
  #' 
  #' @param expression.init: a scRNA-seq raw-count matrix.
  #' @param params: a list of parameters.
  #' @param random_state: a numeric.
  #' 
  #' @return a list of three elements: 'peakRAM', 'time' and 'preds'.
  #' 
  time.before <- Sys.time()
  memory_summary <- gc(reset=TRUE)
  peakRAM.before <- memory_summary[11] + memory_summary[12]
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  results <- do_scEVE(expression.init, params, figures=FALSE, random_state=random_state)
  preds <- results$preds
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  time.after <- Sys.time()
  memory_summary <- gc()
  peakRAM.after <- memory_summary[11] + memory_summary[12]
  
  results <- list(peakRAM = peakRAM.after - peakRAM.before,
                  time = as.numeric(time.after - time.before, units="secs"),
                  preds = preds)
  return(results)
}

N_HVGs_FOR_INDIVIDUAL_METHODS=5000
RANDOM_STATE=0

####################
# init the benchmark
####################
expression.init <- read.csv("./data/Jerby-Arnon_MLM.csv", header=TRUE, row.names=1)
ground_truth <- sapply(strsplit(colnames(expression.init), split="_"), "[", 1)
expression.count <- get_expression.count(expression.init, N_HVGs)
SeurObj.count <- get_SeurObj.count(expression.count)

data <- list(expr=expression.count, seurobj=SeurObj.count, 
             ground=ground_truth, time=0, peakRAM=0)
bd <- BenchDesign(data=data)

########################
# add individual methods
########################
bd <- bd %>%
  addMethod(label="Seurat",
            func=benchmark_Seurat,
            post=list(labels=get_labels, time=save_time, peakRAM=save_peakRAM),
            params=rlang::quos(SeurObj.count=seurobj,
                               random_state=RANDOM_STATE)) %>%
  addMethod(label="densityCut",
            func=benchmark_densityCut,
            post=list(labels=get_labels, time=save_time, peakRAM=save_peakRAM),
            params=rlang::quos(expression.count=expr,
                               random_state=RANDOM_STATE)) %>%
  addMethod(label="monocle3",
            func=benchmark_monocle3,
            post=list(labels=get_labels, time=save_time, peakRAM=save_peakRAM),
            params=rlang::quos(SeurObj.count=seurobj,
                               random_state=RANDOM_STATE)) %>%
  addMethod(label="SHARP",
            func=benchmark_SHARP,
            post=list(labels=get_labels, time=save_time, peakRAM=save_peakRAM),
            params=rlang::quos(expression.count=expr,
                               random_state=RANDOM_STATE))

######################
# add ensemble methods
######################
bd <- bd %>%
  addMethod(label="SAME",
            func=benchmark_SAME,
            post=list(labels=get_labels, time=save_time, peakRAM=save_peakRAM),
            params=rlang::quos(expression.count=expr,
                               SeurObj.count=seurobj,
                               clustering_methods=get_default_hyperparameters()$clustering_methods,
                               random_state=RANDOM_STATE,
                               criterion="BIC"))

# run the benchmark
###################
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
