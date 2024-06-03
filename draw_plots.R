"Run this script to generate the plots used in the scEVE paper.

	2024/06/03 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})

source("./src/paper/plots.R")

results <- get_results()
results.real <- results[results$real,]
results.synthetic <- results[!results$real,]

# performances of individual methods and scEVE on real datasets
for (metric in c("ARI", "NMI", "peakRAM**", "time**")) {
  plot.real <- get_plot.real(results.real, metric)
  plot.synthetic <- get_plot.synthetic(results.synthetic, metric)
  ggsave(glue("./plots/real/{metric}.png"), plot.real, width=4, height=7)
  ggsave(glue("./plots/synthetic/{metric}.png"), plot.synthetic)
}

# performances of scEVE w.r.t. previous ensemble clustering algorithms
barplots.ensemble <- get_barplots.ensemble(results.real)
ggsave("./plots/ensemble.png", barplots.ensemble, width=7, height=2.3)