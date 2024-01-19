"Generate a records .xlsx file.

  Expected global variables:
    - DATASET: a complete scRNA-seq dataset.
    - PARAMS: a named list of hyper-parameters, with:
      -> 'root_consensus': a numeric.

	2023/12/12 @yanisaspic"

source("./src/utils/misc.R")

init_records(DATASET, PARAMS)
