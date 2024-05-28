"Run this script to generate the figures used in the scEVE paper.

	2024/05/27 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})

source("./src/paper/figures.R")

results <- get_results()
results.real <- results[results$real,]
plot.real <- get_plot.real(results.real, "ARI")
for (metric in c("ARI", "NMI", "peakRAM**", "time**")) {
  plot.real <- get_plot.real(results.real, metric)
  ggsave(glue("./figures/{metric}.png"), plot.real)
}