"Run this script to generate the benchmark of scEVE.

	2023/01/20 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(SummarizedBenchmark)
})

source("./scEVE.R")
source("./src/benchmark/metrics.R")
source("./src/benchmark/methods.R")

N_HVGs=5000
RANDOM_STATE=0

####################
# init the benchmark
####################
expression.init <- read.csv("./data/Darmanis_GBM.csv", header=TRUE, row.names=1)
ground_truth <- sapply(strsplit(colnames(expression.init), split="_"), "[", 1)
expression.count <- get_expression.count(expression.init, N_HVGs)
SeurObj.count <- get_SeurObj.count(expression.count)

data <- list(expr=expression.count, seurobj=SeurObj.count, 
             ground=ground_truth, time=0, memory=0)
bd <- BenchDesign(data=data)

########################
# add individual methods
########################
bd <- bd %>%
  addMethod(label="Seurat",
            func=benchmark_Seurat,
            post=list(labels=get_labels, time=save_time, memory=save_memory),
            params=rlang::quos(SeurObj.count=seurobj,
                               random_state=RANDOM_STATE)) %>%
  addMethod(label="densityCut",
            func=benchmark_densityCut,
            post=list(labels=get_labels, time=save_time, memory=save_memory),
            params=rlang::quos(expression.count=expr,
                               random_state=RANDOM_STATE)) %>%
  addMethod(label="monocle3",
            func=benchmark_monocle3,
            post=list(labels=get_labels, time=save_time, memory=save_memory),
            params=rlang::quos(SeurObj.count=seurobj,
                               random_state=RANDOM_STATE)) %>%
  addMethod(label="SHARP",
            func=benchmark_SHARP,
            post=list(labels=get_labels, time=save_time, memory=save_memory),
            params=rlang::quos(expression.count=expr,
                               random_state=RANDOM_STATE))

# run the benchmark
###################
sb <- buildBench(bd=bd, truthCols=c(labels="ground", time="time", memory="memory"))

#############
# add metrics
#############
sb <- sb %>%
  addPerformanceMetric(assay="labels", evalMetric="ARI", evalFunction=get_ARI) %>%
  addPerformanceMetric(assay="labels", evalMetric="NMI", evalFunction=get_NMI) %>%
  addPerformanceMetric(assay="time", evalMetric="time", evalFunction=get_time) %>%
  addPerformanceMetric(assay="memory", evalMetric="memory", evalFunction=get_memory)
metrics <- estimatePerformanceMetrics(sb, tidy=TRUE)
