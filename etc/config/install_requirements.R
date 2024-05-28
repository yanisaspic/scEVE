"Run this script to install the R dependencies.

	2024/04/11 @yanisaspic"

if (!require("devtools")) install.packages("devtools", dependencies=T)

#___________________________________________________________________________CRAN
CRAN <- "http://cran.us.r-project.org"
if (!require("aricode")) install_version("aricode", version="1.0.3", repos=CRAN, dependencies=T)
if (!require("arules")) install_version("arules", version="1.7.7", repos=CRAN, dependencies=T)
if (!require("BiocManager")) install_version("BiocManager", version="1.30.22", repos=CRAN, dependencies=T)
if (!require("dplyr")) install_version("dplyr", version="1.1.4", repos=CRAN, dependencies=T)
if (!require("egg")) install_version("egg", version="0.4.5", repos=CRAN, dependencies=T)
if (!require("ggplot2")) install_version("ggplot2", version="3.4.4", repos=CRAN, dependencies=T)
if (!require("ggplotify")) install_version("ggplotify", version="0.1.2", repos=CRAN, dependencies=T)
if (!require("ggVennDiagram")) install_version("ggVennDiagram", version="1.5.0", repos=CRAN, dependencies=T)
if (!require("glue")) install_version("glue", version="1.7.0", repos=CRAN, dependencies=T)
if (!require("gridExtra")) install_version("gridExtra", version="2.3", repos=CRAN, dependencies=T)
if (!require("igraph")) install_version("igraph", version="1.6.0", repos=CRAN, dependencies=T)
if (!require("Matrix")) install_version("Matrix", version="1.6.5", repos=CRAN, dependencies=T)
if (!require("openxlsx")) install_version("openxlsx", version="4.2.5.2", repos=CRAN, dependencies=T)
if (!require("qpdf")) install_version("qpdf", version="1.3.2", repos=CRAN, dependencies=T)
if (!require("RColorBrewer")) install_version("RColorBrewer", version="1.1.3", repos=CRAN, dependencies=T)
if (!require("reshape2")) install_version("reshape2", version="1.4.4", repos=CRAN, dependencies=T)
if (!require("rlang")) install_version("rlang", version="1.1.3", repos=CRAN, dependencies=T)
if (!require("SCpubr")) install_version("SCpubr", version="2.0.2", repos=CRAN, dependencies=T)
if (!require("Seurat")) install_version("Seurat", version="5.0.1", repos=CRAN, dependencies=T)
if (!require("SeuratObject")) install_version("SeuratObject", version="5.0.1", repos=CRAN, dependencies=T)
if (!require("sp")) install_version("sp", version="2.1.2", repos=CRAN, dependencies=T)


#___________________________________________________________________Bioconductor
if (!require("BiocGenerics")) install("BiocGenerics", dependencies=T, update=F)
if (!require("GenomeInfoDb")) install("GenomeInfoDb", dependencies=T, update=F)
if (!require("GenomicRanges")) install("GenomicRanges", dependencies=T, update=F)
if (!require("IRanges")) install("IRanges", dependencies=T, update=F)
if (!require("MatrixGenerics")) install("MatrixGenerics", dependencies=T, update=F)
if (!require("S4Vectors")) install("S4Vectors", dependencies=T, update=F)
if (!require("scater")) install("scater", dependencies=T, update=F)
if (!require("scuttle")) install("scuttle", dependencies=T, update=F)
if (!require("SingleCellExperiment")) install("SingleCellExperiment", dependencies=T, update=F)
if (!require("SummarizedExperiment")) install("SummarizedExperiment", dependencies=T, update=F)
if (!require("TMExplorer")) install("TMExplorer", dependencies=T, update=F)


#______________________________________________________Bitbucket, Github, Gitlab
if (!require("densitycut")) install_bitbucket("jerry00/densitycut_dev", dependencies=T)
if (!require("monocle3")) install_github('cole-trapnell-lab/monocle3', dependencies = T)
if (!require("SeuratWrappers")) install_github('satijalab/seurat-wrappers', dependencies = T)
if (!require("SHARP")) install_github("shibiaowan/SHARP", dependencies = T)
if (!require("SPARSim")) install_gitlab("sysbiobig/sparsim", depencencies = T)