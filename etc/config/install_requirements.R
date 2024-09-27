"Run this script to install the R dependencies.

	2024/04/11 @yanisaspic"

if (!require("devtools")) install.packages("devtools", dependencies=T)
CRAN <- "http://cran.us.r-project.org"

#___________________________________________________________________________CRAN
if (!require("aricode")) install_version("aricode", version="1.0.3", repos=CRAN, dependencies=T)
if (!require("arules")) install_version("arules", version="1.7.7", repos=CRAN, dependencies=T)
if (!require("dbplyr")) install_version("dbplyr", version="2.5.0", repos=CRAN, dependencies=T)
if (!require("dplyr")) install_version("dplyr", version="1.1.4", repos=CRAN, dependencies=T)
if (!require("egg")) install_version("egg", version="0.4.5", repos=CRAN, dependencies=T)
if (!require("ggh4x")) install_version("ggh4x", version="0.2.8", repos=CRAN, dependencies=T)
if (!require("ggplot2")) install_version("ggplot2", version="3.5.1", repos=CRAN, dependencies=T)
if (!require("ggplotify")) install_version("ggplotify", version="0.1.2", repos=CRAN, dependencies=T)
if (!require("ggVennDiagram")) install_version("ggVennDiagram", version="1.5.2", repos=CRAN, dependencies=T)
if (!require("glue")) install_version("glue", version="1.7.0", repos=CRAN, dependencies=T)
if (!require("gridExtra")) install_version("gridExtra", version="2.3", repos=CRAN, dependencies=T)
if (!require("igraph")) install_version("igraph", version="2.0.3", repos=CRAN, dependencies=T)
if (!require("Matrix")) install_version("Matrix", version="1.6-5", repos=CRAN, dependencies=T)
if (!require("matrixStats")) install_version("matrixStats", version="1.4.1", repos=CRAN, dependencies=T)
if (!require("openxlsx")) install_version("openxlsx", version="4.2.5.2", repos=CRAN, dependencies=T)
if (!require("patchwork")) install_version("patchwork", version="1.3.0", repos=CRAN, dependencies=T)
if (!require("qpdf")) install_version("qpdf", version="1.3.3", repos=CRAN, dependencies=T)
if (!require("RColorBrewer")) install_version("RColorBrewer", version="1.1-3", repos=CRAN, dependencies=T)
if (!require("readxl")) install_version("readxl", version="1.4.3", repos=CRAN, dependencies=T)
if (!require("reshape2")) install_version("reshape2", version="1.4.4", repos=CRAN, dependencies=T)
if (!require("rlang")) install_version("rlang", version="1.1.4", repos=CRAN, dependencies=T)
if (!require("scales")) install_version("scales", version="1.3.0", repos=CRAN, dependencies=T)
if (!require("SCpubr")) install_version("SCpubr", version="2.0.2", repos=CRAN, dependencies=T)
if (!require("sessioninfo")) install_version("sessioninfo", version="1.2.2", repos=CRAN, dependencies=T)
if (!require("Seurat")) install_version("Seurat", version="5.1.0", repos=CRAN, dependencies=T)
if (!require("SeuratObject")) install_version("SeuratObject", version="5.0.2", repos=CRAN, dependencies=T)
if (!require("sp")) install_version("sp", version="2.1-4", repos=CRAN, dependencies=T)
if (!require("tidytree")) install_version("tidytree", version="0.4.6", repos=CRAN, dependencies=T)

#___________________________________________________________________Bioconductor
if (!require("BiocManager")) install_version("BiocManager", version="1.30.23", repos=CRAN, dependencies=T)
BiocManager::install(version="3.18")

if (!require("Biobase")) install("Biobase", dependencies=T, update=F)
# packageVersion("Biobase") = 2.62.0
if (!require("BiocFileCache")) install("BiocFileCache", dependencies=T, update=F)
# packageVersion("BiocFileCache") = 2.10.2
if (!require("BiocGenerics")) install("BiocGenerics", dependencies=T, update=F)
# packageVersion("BiocGenerics") = 0.48.1
if (!require("GenomeInfoDb")) install("GenomeInfoDb", dependencies=T, update=F)
# packageVersion("GenomeInfoDb") = 1.38.8
if (!require("GenomicRanges")) install("GenomicRanges", dependencies=T, update=F)
# packageVersion("GenomicRanges") = 1.54.1
if (!require("ggtree")) {install("YuLab-SMU/treedataverse", dependencies=T, update=F)}
# packageVersion("ggtree") = 3.10.1
if (!require("IRanges")) install("IRanges", dependencies=T, update=F)
# packageVersion("IRanges") = 2.36.0
if (!require("MatrixGenerics")) install("MatrixGenerics", dependencies=T, update=F)
# packageVersion("MatrixGenerics") = 1.14.0
if (!require("S4Vectors")) install("S4Vectors", dependencies=T, update=F)
# packageVersion("S4Vectors") = 0.40.2
if (!require("scater")) install("scater", dependencies=T, update=F)
# packageVersion("scater") = 1.30.1
if (!require("scuttle")) install("scuttle", dependencies=T, update=F)
# packageVersion("scuttle") = 1.12.0
if (!require("SingleCellExperiment")) install("SingleCellExperiment", dependencies=T, update=F)
# packageVersion("SingleCellExperiment") = 1.24.0
if (!require("SummarizedExperiment")) install("SummarizedExperiment", dependencies=T, update=F)
# packageVersion("SummarizedExperiment") = 1.32.0
if (!require("TMExplorer")) install("TMExplorer", dependencies=T, update=F)
# packageVersion("TMExplorer") = 1.12.0

#______________________________________________________Bitbucket, Github, Gitlab
if (!require("densitycut")) install_bitbucket("jerry00/densitycut_dev@3878b5f", dependencies=T)
# packageVersion("densityCut") = 0.0.1
if (!require("monocle3")) install_github("cole-trapnell-lab/monocle3@98402ed", dependencies=T)
# packageVersion("monocle3") = 1.3.7
if (!require("SeuratWrappers")) install_github("satijalab/seurat-wrappers@28f074c", dependencies=T)
# packageVersion("SeuratWrappers") = 0.3.4
if (!require("SHARP")) install_github("shibiaowan/SHARP@24a9866", dependencies=T)
# packageVersion("SHARP") = 1.1.0
if (!require("SPARSim")) install_gitlab("sysbiobig/sparsim@19f86cb", depencencies=T)
# packageVersion("SPARSim") = 0.9.5