"Run this script to get the R requirements and their versions.

	2024/09/27 @yanisaspic"

source("./scEVE.R")
source("./src/paper/data.R")
source("./src/paper/methods.R")
source("./src/paper/metrics.R")
source("./src/paper/plots.R")

suppressPackageStartupMessages({
  library(sessioninfo)
})

packageVersion("sessioninfo")
sessioninfo::session_info(to_file="session.log")