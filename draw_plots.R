"Run this script to generate the plots used in the scEVE paper.

	2024/06/03 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(openxlsx) # packageVersion("openxlsx")==4.2.5.2
})

source("./src/paper/plots.R")
source("./src/paper/metrics.R")

benchmark <- get_benchmark()

# performances of individual methods and scEVE on real datasets.
for (metric in c("ARI", "NMI", "log10(s)", "log10(Mb)")) {
  plot.real <- get_plot.real(benchmark[benchmark$real,], metric)
  ggsave(glue("./plots/real/{metric}.png"), plot.real, width=4.5, height=7)
}

# performances of scEVE w.r.t. previous ensemble clustering algorithms.
for (metric in c("ARI", "NMI")) {
  plot.ensemble <- get_plot.ensemble(benchmark[benchmark$real,], metric)
  ggsave(glue("./plots/ensemble/{metric}.png"), plot.ensemble, width=10, height=7)
}

# trees of real scRNA-seq datasets clustering.
for (dataset in get_prior()$real_datasets) {
  records <- get_records(glue("./results/records/{dataset}.xlsx"))
  resolution_tree.data <- get_resolution_tree.data(records$meta)
  distributions.data <- get_distributions.data(records$cells)
  tree <- get_plot.tree(resolution_tree.data, distributions.data, dataset)
  ggsave(glue("./plots/trees/{dataset}.png"), tree, width=13, height=5.5)
}

# performances of individual methods and scEVE on synthetic datasets.
for (method in get_prior()$algorithm) {
  for (metric in c("ARI", "NMI")) {
    plot.synthetic <- get_plot.synthetic.method(benchmark.synthetic, metric, method)
    ggsave(glue("./plots/synthetic/{method}_{metric}.png"), plot.synthetic, width=4, height=4)
  }
}

# focus on the nested-balanced synthetic datasets.
for (metric in c("ARI", "NMI")) {
  plot.synthetic <- get_plot.synthetic.scenario(benchmark.synthetic, metric, TRUE, TRUE)
  ggsave(glue("./plots/synthetic/NB_{metric}.png"), plot.synthetic, width=4.5, height=7)  
}

# disruption and consensi of clusters predicted by scEVE on real datasets.
# real_datasets <- get_prior()$real_datasets
# disruptions <- get_disruption(real_datasets)
# test <- disruptions %>% mutate(bins=cut(consensus, breaks=10), )
# # plot <- ggplot(data=test) +
# #   geom_boxplot(aes(x=bins, y=disruption, group=bins, fill=bins)) +
# #   geom_point(aes(x=bins, y=disruption, fill=bins), position=position_jitterdodge(), alpha=.5) 
# 
# print(cons)
# ground_truth <- get_ground_truth(rownames(sheet.cells))
# distribution.population <- get_distribution.population("C.3", sheet.cells, ground_truth)