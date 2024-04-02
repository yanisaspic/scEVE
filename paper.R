"Run this script to generate the results of the paper submitted to JOBIM2024.

	2024/04/02 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})
source("./scEVE.R")
source("./src/paper/methods.R")
source("./src/paper/metrics.R")

#_______________________________________________________________________showcase
expression.init <- read.csv("./data/datasets/Darmanis_HumGBM.csv", header=TRUE, row.names=1)
ground_truth <- get_ground_truth(expression.init)
output <- do_scEVE(expression.init)
figure_2.tree <- draw_tree(output$records, ground_truth)
figure_2.heatmap <- draw_heatmap(output$preds, ground_truth)

#______________________________________________________________________benchmark
datasets <- c("./data/datasets/Baron_HumPan.csv",
              "./data/datasets/Li_HumCRC.csv")
for (ds in datasets) {
  expression.init <- read.csv(ds, header=TRUE, row.names=1)
  ground_truth <- get_ground_truth(expression.init)
  results <- get_benchmark.scEVE(expression.init, params=get_default_hyperparameters(), random_state=0)
  print("ARI:")
  print(get_metric(preds=results$preds, ground_truth=ground_truth, metric="ARI"))
  print("NMI:")
  print(get_metric(preds=results$preds, ground_truth=ground_truth, metric="NMI"))
  print("time:")
  print(results$time)
  print("peakRAM")
  print(results$peakRAM)
}
