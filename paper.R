"Run this script to generate the results of the paper submitted to JOBIM2024.

	2024/04/03 @yanisaspic"

source("./scEVE.R")
source("./src/paper/methods.R")
source("./src/paper/metrics.R")

#_______________________________________________________________________showcase
expression.init <- read.csv("./data/Baron_HumPan.csv", header=TRUE, row.names=1)
ground_truth <- get_ground_truth(expression.init)
output <- do_scEVE(expression.init)

draw_tree(output$records, ground_truth)
figure_2.heatmap <- draw_heatmap(output$preds, ground_truth)
print(figure_2.heatmap)
# figure_2: tree is manually annotated
# figure_3: see './figures/C5.pdf'
# figure_4: markers extracted from './records.xlsx'

#______________________________________________________________________benchmark
get_benchmark.dataset <- function(dataset) {
  expression.init <- read.csv(dataset, header=TRUE, row.names=1)
  ground_truth <- get_ground_truth(expression.init)
  results <- get_benchmark.scEVE(expression.init, params=get_default_hyperparameters(), random_state=0)
  metrics <- list(ARI=get_metric(results$preds, ground_truth, metric="ARI"),
                  NMI=get_metric(results$preds, ground_truth, metric="NMI"),
                  time=results$time,
                  peakRAM=results$peakRAM)
  return(metrics)
}

datasets <- list(Baron="./data/Baron_HumPan.csv",
                 Li="./data/Li_HumCRC.csv",
                 Darmanis="./data/Darmanis_HumGBM.csv")
benchmark <- lapply(X=datasets, FUN=get_benchmark.dataset)
benchmark.table <- do.call(rbind, benchmark)
print(benchmark.table)
