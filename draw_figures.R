"Run this script to generate the figures used in the scEVE paper.

	2024/05/27 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})

source("./src/paper/figures.R")

results <- get_results()
results.real <- results[results$real,]
results.synthetic <- results[!results$real,]

# performances of individual methods and scEVE on real datasets
for (metric in c("ARI", "NMI", "peakRAM**", "time**")) {
  plot.real <- get_plot.real(results.real, metric)
  ggsave(glue("./figures/{metric}.png"), plot.real)
}

# performances of scEVE w.r.t. previous ensemble clustering algorithms
barplots.ensemble <- get_barplots.ensemble(results.real)
ggsave("./figures/ensemble.png", barplots.ensemble, width=7, height=2.3)