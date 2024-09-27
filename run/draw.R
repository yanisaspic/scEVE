"Run this script to generate the plots used in the scEVE paper.

	2024/06/03 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
  library(openxlsx) # packageVersion("openxlsx")==4.2.5.2
})

source("./src/paper/plots.R")
source("./src/paper/metrics.R")

benchmark <- get_results.benchmark()
benchmark.synthetic <- setup_data.synthetic(benchmark[!benchmark$real, ])
similarities <- get_results("./results/similarity")
similarities <- setup_data.synthetic(similarities)

# performances of individual methods and scEVE on real datasets.
plot.real.summary <- get_plot.real.summary(benchmark[benchmark$real,])
ggsave("./plots/real/summary.png", plot.real.summary, width=12.5, height=8)
for (metric in c("ARI", "NMI", "log10(s)", "log10(Mb)")) {
  plot.real.metric <- get_plot.real(benchmark[benchmark$real,], metric)
  ggsave(glue("./plots/real/{metric}.png"), plot.real.metric, width=4.5, height=7)
}

# performances of scEVE w.r.t. previous ensemble clustering algorithms.
for (metric in c("ARI", "NMI")) {
  plot.ensemble <- get_plot.ensemble(benchmark[benchmark$real,], metric)
  ggsave(glue("./plots/ensemble/{metric}.png"), plot.ensemble, width=10, height=6)
}

# trees of real scRNA-seq datasets clustering.
for (dataset in get_prior()$real_datasets) {
  records <- get_records(glue("./results/records/{dataset}.xlsx"))
  resolution_tree.data <- get_resolution_tree.data(records$meta)
  distributions.data <- get_distributions.data(records$cells)
  tree <- get_plot.tree(resolution_tree.data, distributions.data, dataset)
  ggsave(glue("./plots/trees/{dataset}.png"), tree, width=13, height=5.5)
}

# barplot of cancer signatures in the Darmanis dataset
records <- get_records(glue("./results/records/Darmanis_HumGBM.xlsx"))
cancer_signatures <- read.csv("./data/cancer_signatures.csv")
signatures.data <- setup_signatures.data(records$markers, cancer_signatures)
plot.signatures <- get_plot.signatures(signatures.data)
ggsave("./plots/signatures.png", plot.signatures, width=4, height=6)

# performances of scEVE on synthetic datasets, and explanation.
for (metric in c("ARI", "NMI")) {
  plot.synthetic.scEVE <- get_plot.synthetic.method(benchmark.synthetic, metric, method="scEVE")
  plot.synthetic.hard <- get_plot.synthetic.setting(benchmark.synthetic, metric,
                                                    related="yes", balanced="yes")
  ggsave(glue("./plots/synthetic/scEVE_{metric}.png"), plot.synthetic.scEVE, width=5, height=5)
  ggsave(glue("./plots/synthetic/hard_{metric}.png"), plot.synthetic.hard, width=4.5, height=6.5)
}

# performances of individual methods on synthetic datasets.
for (metric in c("ARI", "NMI")) {
  plot.synthetic.summary <- get_plot.synthetic.summary(benchmark.synthetic, metric)
  ggsave(glue("./plots/synthetic/summary_{metric}.png"), plot.synthetic.summary,
         width=9.5, height=10.5)
}

# mean (+ std) performance of the clustering methods on real datasets
data <- benchmark[(benchmark$real) & (benchmark$method != "scEVE"), ]
performances <- data %>%
  group_by(dataset) %>%
  summarise(ARI=mean(ARI), NMI=mean(NMI))
print(performances)

# mean (+ std) similarity of the predictions of the clustering methods on synthetic datasets
similarities.summary <- similarities %>%
  group_by(n_populations) %>%
  summarise(mean_ARI=mean(ARI), std_ARI=sd(ARI), mean_NMI=mean(NMI), std_NMI=sd(NMI))
print(similarities.summary)