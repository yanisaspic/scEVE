"Dependencies of the R scripts.

	2023/07/19 @yanisaspic"

CRAN <- "http://cran.us.r-project.org"
if (!require("devtools")) install.packages("devtools", dependencies=T)
# packageVersion("devtools")    2.4.5

if (!require("stringr")) install_version("stringr", version="1.5.0", repos=CRAN, dependencies=T)
if (!require("readxl")) install_version("readxl", version="1.4.3", repos=CRAN, dependencies=T)
if (!require("qpdf")) install_version("qpdf", version="1.3.2", repos=CRAN, dependencies=T)
if (!require("BiocManager")) install_version("BiocManager", version="1.30.20", repos=CRAN, dependencies=T)
if (!require("aricode")) install_version("aricode", version="1.0.3", repos=CRAN, dependencies=T)
if (!require("rjson")) install_version("rjson", version="0.2.21", repos=CRAN, dependencies=T)
if (!require("lubridate")) install_version("lubridate", version="1.9.2", repos=CRAN, dependencies=T)
if (!require("SCpubr")) install_version("SCpubr", version="2.0.2", repos=CRAN, dependencies=T)
if (!require("ggplotify")) install_version("ggplotify", version="0.1.2", repos=CRAN, dependencies=T)

if (!require("BiocGenerics")) install("BiocGenerics", dependencies=T, update=F)
if (!require("DelayedArray")) install("DelayedArray", dependencies=T, update=F) 
if (!require("DelayedMatrixStats")) install("DelayedMatrixStats", dependencies=T, update=F)
if (!require("limma")) install("limma", dependencies=T, update=F)
if (!require("lme4")) install("lme4", dependencies=T, update=F)
if (!require("S4Vectors")) install("S4Vectors", dependencies=T, update=F)
if (!require("SingleCellExperiment")) install("SingleCellExperiment", dependencies=T, update=F)
if (!require("batchelor")) install("batchelor", dependencies=T, update=F)
if (!require("HDF5Array")) install("HDF5Array", dependencies = T, update=F)
if (!require("ggrastr")) install("ggrastr", dependencies = T, update=F)
if (!require("spdep")) install("spdep", dependencies = T, update=F)
# packageVersion("BiocGenerics")          0.46.0
# packageVersion("DelayedArray")          0.26.6
# packageVersion("DelayedMatrixStats")    1.22.1
# packageVersion("limma")                 3.56.2
# packageVersion("lme4")                  1.1.34
# packageVersion("S4Vectors")             0.38.1
# packageVersion("SingleCellExperiment")  1.22.0
# packageVersion("batchelor")             1.16.0
# packageVersion("HDF5Array")             1.28.1
# packageVersion("ggrastr")               1.0.2
# packageVersion("spdep")                 1.2.8
if (!require("cluster")) install_version("cluster", version="2.1.4", repos=CRAN, dependencies = T)
if (!require("mclust")) install_version("mclust", version="6.0.0", repos=CRAN, dependencies = T)
if (!require("RMTstat")) install_version("RMTstat", version="0.3.1", repos=CRAN, dependencies = T)
if (!require("clusterCrit")) install_version("clusterCrit", version="1.2.8", dependencies = T)
if (!require("tensorflow")) {
  install_version("tensorflow", version="2.11.0", dependencies = T)
  install_tensorflow()
}

if (!require("scater")) install("scater", dependencies=T, update=F)
if (!require("Seurat")) install_version("Seurat", version="4.3.0.1", repos=CRAN, dependencies=T)
if (!require("SeuratWrappers")) install_github('satijalab/seurat-wrappers', dependencies=T)
if (!require("monocle3")) install_github('cole-trapnell-lab/monocle3', dependencies = T)
if (!require("SHARP")) install_github("shibiaowan/SHARP", dependencies=T)
if (!require("densitycut")) install_bitbucket("jerry00/densitycut_dev", dependencies=T)
if (!require("Wind")) install_github("haowulab/Wind", dependencies=T)
if (!require("scMayoMap")) install_github("chloelulu/scMayoMap", dependencies=T)
if (!require("ggVennDiagram")) install_github("gaospecial/ggVennDiagram", dependencies=T)
# packageVersion("scater")          1.28.0
# packageVersion("splatter")        1.22.1
# packageVersion("SeuratWrappers")  0.3.1
# packageVersion("monocle3")        1.3.1
# packageVersion("cidr")            0.1.5
# packageVersion("SIMLR")           1.26.1
# packageVersion("scCCESS")         0.3.2
# packageVersion("SHARP")           1.1.0
# packageVersion("densitycut")      0.0.1
# packageVersion("Wind")            0.9.1
# packageVersion("ggVennDiagram")   1.4.13