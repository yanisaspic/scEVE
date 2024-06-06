"Run this script to generate the plots used in the scEVE paper.

	2024/06/03 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})

source("./src/paper/plots.R")

results <- get_results()
results.real <- results[results$real,]
results.synthetic <- results[!results$real,]
records <- get_records()

# performances of individual methods and scEVE on real datasets
for (metric in c("ARI", "NMI", "log10(s)", "log10(Mb)")) {
  plot.real <- get_plot.real(results.real, metric)
  ggsave(glue("./plots/real/{metric}.png"), plot.real, width=4.5, height=7)
}

# performances of scEVE w.r.t. previous ensemble clustering algorithms
barplots.ensemble <- get_barplots.ensemble(results.real)
ggsave("./plots/ensemble.png", barplots.ensemble, width=8.5, height=2)