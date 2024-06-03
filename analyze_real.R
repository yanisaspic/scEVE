"Run this script to analyze the real scRNA-seq datasets with figures and records.

	2024/06/03 @yanisaspic"

suppressPackageStartupMessages({
  library(glue)
})

source("./scEVE.R")
source("./src/paper/data.R")

store_results <- function(dataset) {
  #' Store the results of a scEVE analysis in a directory named after the dataset analyzed.
  #' 
  #' @param dataset: a character.
  #' 
  dir.create(glue("./analysis/figures/{dataset}"))
  figures <- list.files("./figures")
  for (fig in figures) {
    file.rename(glue("./figures/{fig}"), glue("./analysis/figures/{dataset}/{fig}"))
  }
  file.rename("./records.xlsx", glue("./analysis/records/{dataset}.xlsx"))
}

real_datasets <- c(
  "Li_HumCRC_b",         # 364c
  "Li_HumCRC_a",         # 561c
  "Baron_MouPan_1",      # 822c
  "Baron_MouPan_2",      # 1,064c
  "Baron_HumPan_4",      # 1,303c
  "Tasic_MouBra",        # 1,679c
  "Baron_HumPan_2",      # 1,724c
  "Baron_HumPan_1",      # 1,937c
  "Darmanis_HumGBM",     # 3,589c
  "Baron_HumPan_3",      # 3,605c
  "JerbyArnon_HumMLM",   # 6,879c
  "VanGalen_HumAML",     # 27,899c
  "Lambrechts_HumNSCLC", # 51,775c
  "Peng_HumPDAC"         # 57,530c
)

params <- get_default_hyperparameters()
for (dataset in real_datasets) {
  expression.init <- get_real_dataset(dataset)
  results <- do_scEVE(expression.init=expression.init, params=params,
                      figures=TRUE, random_state=0, save=TRUE)
  store_results(dataset)
}
