# scEVE: a scRNA-seq ensemble clustering algorithm

[![R](https://img.shields.io/badge/R-4.3.1-blue.svg)](https://cran.r-project.org/) [![Ubuntu](https://img.shields.io/badge/20.04.6-black.svg)](https://ubuntu.com/) [![bash](https://img.shields.io/badge/5.0.17-black.svg)](https://www.gnu.org/software/bash/) [![Python](https://img.shields.io/badge/Python-3.9-blue.svg)](https://www.python.org/) 

![The figure is composed of three parts: a multi-resolution cluster tree (left), a cell type plot (middle), and an abundance plot (right).\\
The multi-resolution cluster tree (left) represents the cell populations predicted by scEVE as nodes. Edges between nodes connect sub-populations to their parents, and their weight indicates the robustness of the sub-populations. The nodes representing the original cell pool and the leftover clusters are black. The cell type plot (middle) associates the highest-resolution populations to bars. If a population is homogeneous (\textit{i.e.} composed of a unique cell type), its bar is monochrome. Otherwise, its bar is polychrome. The abundance plot (right) also associates highest-resolution populations to bars. The bars are log-scaled, and they indicate the size of the cell populations.](/plots/Darmanis_HumGBM.png?raw=true "Clustering results of scEVE")

___
*Single-cell RNA-sequencing is used to measure the individual transcriptomes of cells within a tissue. In a decade only, it motivated the development of hundreds of clustering methods that attempt to identify cell populations with similar transcriptomes. Because each method relies on its specific hypotheses, the results of a clustering analysis can vary drastically based on the method used, which makes employing a single method ill-advised. To address that issue, ensemble algorithms integrated multiple clustering methods by minimizing the differences in their results. While that approach is sensible, it does not address some methodological challenges in single-cell data science: namely, the need to generate clustering results with multiple cellular resolutions and explicit uncertainty values. In this work, we present a novel approach to ensemble clustering that addresses these challenges by leveraging the differences between clustering results. We present our algorithm scEVE, and we evaluate it on 15 experimental and up to 600 synthetic datasets. Our results highlight how beneficial these challenges are to the biological downstream analysis, and they show that scEVE is the best ensemble algorithm to address them, yet. Finally, we argue that future ensemble algorithm would profit from combining both approaches.*
___
The repository includes the code of the scEVE algorithm, and the codes used to generate the datasets and the results in the related paper "**The differences of results between clustering methods are informative in single-cell RNA-seq ensemble clustering analyses.**" 

## Architecture

# 1. Dependencies
The dependencies required to run the scEVE algorithm and to generate our results are summarized at the end of the README.
The summary is generated with `session_info()` from the R package `sessioninfo` (see https://github.com/r-lib/sessioninfo).

The dependencies are also listed in the `etc/config/` directory. Run the script `etc/config/install_requirements.R` to install them.
```bash
Rscript ./etc/config/install_requirements.R
```

# 2. Data
Run the script `setup.sh` to download and set-up the data used in the paper. The generated files are stored in the `data/` directory.
```bash
chmod +x setup.sh
./setup.sh
```

## Summary of the dependencies
The summary is generated with `session_info()` from the R package `sessioninfo` (see https://github.com/r-lib/sessioninfo).
```
─ Session info ────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.3.3 (2024-02-29)
 os       Ubuntu 22.04.4 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language (EN)
 collate  fr_FR.UTF-8
 ctype    fr_FR.UTF-8
 tz       Europe/Paris
 date     2024-09-06
 rstudio  2022.07.2+576 Spotted Wakerobin (desktop)
 pandoc   2.9.2.1 @ /usr/bin/pandoc

─ Packages ────────────────────────────────────────────────────────────────
 package              * version   date (UTC) lib source
 abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
 ape                    5.7-1     2023-03-13 [1] CRAN (R 4.3.0)
 aplot                  0.2.2     2023-10-06 [1] CRAN (R 4.3.0)
 aricode              * 1.0.3     2023-10-20 [1] CRAN (R 4.3.0)
 arules               * 1.7-7     2023-11-29 [1] CRAN (R 4.3.0)
 askpass                1.2.0     2023-09-03 [1] CRAN (R 4.3.0)
 backports              1.5.0     2024-05-23 [1] CRAN (R 4.3.3)
 beachmat               2.16.0    2023-04-25 [1] Bioconductor
 beeswarm               0.4.0     2021-06-01 [1] CRAN (R 4.3.0)
 Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
 BiocFileCache        * 2.8.0     2023-04-25 [1] Bioconductor
 BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
 BiocManager            1.30.23   2024-05-04 [1] CRAN (R 4.3.3)
 BiocNeighbors          1.18.0    2023-04-25 [1] Bioconductor
 BiocParallel           1.36.0    2023-10-24 [1] Bioconductor
 BiocSingular           1.16.0    2023-04-25 [1] Bioconductor
 bit                    4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
 bit64                  4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
 bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
 blob                   1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
 boot                   1.3-30    2024-02-26 [4] CRAN (R 4.3.3)
 broom                  1.0.5     2023-06-09 [1] CRAN (R 4.3.0)
 cachem                 1.1.0     2024-05-16 [1] CRAN (R 4.3.3)
 car                    3.1-2     2023-03-30 [1] CRAN (R 4.3.0)
 carData                3.0-5     2022-01-06 [1] CRAN (R 4.3.0)
 caTools                1.18.2    2021-03-28 [1] CRAN (R 4.3.0)
 cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.3.0)
 cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.0)
 clipr                  0.8.0     2022-02-22 [1] CRAN (R 4.3.0)
 clues                  0.6.2.2   2019-12-03 [1] url (https://cran.r-project.org/src/contrib/Archive/clues/clues_0.6.2.2.tar.gz)
 cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.0)
 clusterCrit            1.3.0     2023-11-23 [1] CRAN (R 4.3.0)
 codetools              0.2-20    2024-03-31 [1] CRAN (R 4.3.3)
 colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
 cowplot                1.1.2     2023-12-15 [1] CRAN (R 4.3.0)
 crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
 curl                   5.2.1     2024-03-01 [1] CRAN (R 4.3.3)
 data.table             1.15.4    2024-03-30 [1] CRAN (R 4.3.3)
 DBI                    1.2.3     2024-06-02 [1] CRAN (R 4.3.3)
 dbplyr               * 2.3.4     2023-09-26 [1] CRAN (R 4.3.2)
 DelayedArray           0.26.7    2023-07-28 [1] Bioconductor
 DelayedMatrixStats     1.22.6    2023-08-28 [1] Bioconductor
 deldir                 2.0-4     2024-02-28 [1] CRAN (R 4.3.3)
 densitycut           * 0.0.1     2023-07-19 [1] bitbucket (jerry00/densitycut_dev@3878b5f)
 digest                 0.6.35    2024-03-11 [1] CRAN (R 4.3.3)
 doParallel             1.0.17    2022-02-07 [1] CRAN (R 4.3.0)
 dotCall64              1.1-1     2023-11-28 [1] CRAN (R 4.3.0)
 dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.0)
 egg                  * 0.4.5     2019-07-13 [1] CRAN (R 4.3.3)
 ellipsis               0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
 entropy              * 1.3.1     2021-10-02 [1] CRAN (R 4.3.3)
 fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.0)
 fastDummies            1.7.3     2023-07-06 [1] CRAN (R 4.3.0)
 fastmap                1.2.0     2024-05-15 [1] CRAN (R 4.3.3)
 filelock               1.0.3     2023-12-11 [1] CRAN (R 4.3.0)
 fitdistrplus           1.1-11    2023-04-25 [1] CRAN (R 4.3.0)
 flashClust             1.01-2    2012-08-21 [1] CRAN (R 4.3.0)
 foreach                1.5.2     2022-02-02 [1] CRAN (R 4.3.0)
 fs                     1.6.4     2024-04-25 [1] CRAN (R 4.3.3)
 future                 1.33.1    2023-12-22 [1] CRAN (R 4.3.0)
 future.apply           1.11.1    2023-12-21 [1] CRAN (R 4.3.0)
 generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
 GenomeInfoDb         * 1.36.4    2023-10-02 [1] Bioconductor
 GenomeInfoDbData       1.2.11    2024-06-19 [1] Bioconductor
 GenomicRanges        * 1.52.1    2023-10-08 [1] Bioconductor
 ggbeeswarm             0.7.2     2023-04-29 [1] CRAN (R 4.3.0)
 ggfun                  0.1.5     2024-05-28 [1] CRAN (R 4.3.3)
 ggh4x                * 0.2.8     2024-01-23 [1] CRAN (R 4.3.3)
 ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.3.3)
 ggplotify            * 0.1.2     2023-08-09 [1] CRAN (R 4.3.0)
 ggpubr                 0.6.0     2023-02-10 [1] CRAN (R 4.3.0)
 ggrepel                0.9.5     2024-01-10 [1] CRAN (R 4.3.0)
 ggridges               0.5.6     2024-01-23 [1] CRAN (R 4.3.3)
 ggsignif               0.6.4     2022-10-13 [1] CRAN (R 4.3.0)
 ggtree               * 3.10.1    2024-02-25 [1] Bioconductor 3.18 (R 4.3.3)
 ggVennDiagram        * 1.5.0     2024-01-13 [1] CRAN (R 4.3.0)
 globals                0.16.3    2024-03-08 [1] CRAN (R 4.3.3)
 glue                 * 1.7.0     2024-01-09 [1] CRAN (R 4.3.0)
 goftest                1.2-3     2021-10-07 [1] CRAN (R 4.3.0)
 gplots                 3.1.3.1   2024-02-02 [1] CRAN (R 4.3.3)
 gridExtra            * 2.3       2017-09-09 [1] CRAN (R 4.3.0)
 gridGraphics           0.5-1     2020-12-13 [1] CRAN (R 4.3.0)
 gtable                 0.3.5     2024-04-22 [1] CRAN (R 4.3.3)
 gtools                 3.9.5     2023-11-20 [1] CRAN (R 4.3.0)
 htmltools              0.5.8.1   2024-04-04 [1] CRAN (R 4.3.3)
 htmlwidgets            1.6.4     2023-12-06 [1] CRAN (R 4.3.0)
 httpuv                 1.6.13    2023-12-06 [1] CRAN (R 4.3.0)
 httr                   1.4.7     2023-08-15 [1] CRAN (R 4.3.0)
 ica                    1.0-3     2022-07-08 [1] CRAN (R 4.3.0)
 igraph               * 1.6.0     2023-12-11 [1] CRAN (R 4.3.3)
 IRanges              * 2.34.1    2023-06-22 [1] Bioconductor
 irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.0)
 iterators              1.0.14    2022-02-05 [1] CRAN (R 4.3.0)
 jsonlite               1.8.8     2023-12-04 [1] CRAN (R 4.3.0)
 KernSmooth             2.23-24   2024-05-17 [1] CRAN (R 4.3.3)
 later                  1.3.2     2023-12-06 [1] CRAN (R 4.3.0)
 lattice                0.22-6    2024-03-20 [1] CRAN (R 4.3.3)
 lazyeval               0.2.2     2019-03-15 [1] CRAN (R 4.3.0)
 leiden                 0.4.3.1   2023-11-17 [1] CRAN (R 4.3.0)
 lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.0)
 listenv                0.9.1     2024-01-29 [1] CRAN (R 4.3.3)
 lme4                   1.1-35.1  2023-11-05 [1] CRAN (R 4.3.0)
 lmtest                 0.9-40    2022-03-21 [1] CRAN (R 4.3.0)
 magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
 MASS                   7.3-60    2023-05-04 [4] CRAN (R 4.3.1)
 Matrix               * 1.6-5     2024-01-11 [1] CRAN (R 4.3.0)
 MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
 matrixStats          * 1.3.0     2024-04-11 [1] CRAN (R 4.3.3)
 memoise                2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
 mime                   0.12      2021-09-28 [1] CRAN (R 4.3.0)
 miniUI                 0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
 minqa                  1.2.7     2024-05-20 [1] CRAN (R 4.3.3)
 monocle3             * 1.3.4     2023-08-21 [1] Github (cole-trapnell-lab/monocle3@2b17745)
 munsell                0.5.1     2024-04-01 [1] CRAN (R 4.3.3)
 mvtnorm                1.2-5     2024-05-21 [1] CRAN (R 4.3.3)
 nlme                   3.1-163   2023-08-09 [4] CRAN (R 4.3.1)
 nloptr                 2.0.3     2022-05-26 [1] CRAN (R 4.3.0)
 openxlsx             * 4.2.5.2   2023-02-06 [1] CRAN (R 4.3.0)
 parallelly             1.37.1    2024-02-29 [1] CRAN (R 4.3.3)
 patchwork            * 1.2.0     2024-01-08 [1] CRAN (R 4.3.0)
 pbapply                1.7-2     2023-06-27 [1] CRAN (R 4.3.0)
 pheatmap               1.0.12    2019-01-04 [1] CRAN (R 4.3.0)
 pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
 pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
 plotly                 4.10.4    2024-01-13 [1] CRAN (R 4.3.0)
 plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.0)
 png                    0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
 polyclip               1.10-6    2023-09-27 [1] CRAN (R 4.3.0)
 progressr              0.14.0    2023-08-10 [1] CRAN (R 4.3.0)
 promises               1.3.0     2024-04-05 [1] CRAN (R 4.3.3)
 purrr                  1.0.2     2023-08-10 [1] CRAN (R 4.3.0)
 qpdf                 * 1.3.3     2024-03-25 [1] CRAN (R 4.3.3)
 R.methodsS3            1.8.2     2022-06-13 [1] CRAN (R 4.3.0)
 R.oo                   1.26.0    2024-01-24 [1] CRAN (R 4.3.3)
 R.utils                2.12.3    2023-11-18 [1] CRAN (R 4.3.0)
 R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
 RANN                   2.6.1     2019-01-08 [1] CRAN (R 4.3.0)
 RColorBrewer         * 1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
 Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.0)
 RcppAnnoy              0.0.22    2024-01-23 [1] CRAN (R 4.3.3)
 RcppHNSW               0.6.0     2024-02-04 [1] CRAN (R 4.3.3)
 RCurl                  1.98-1.14 2024-01-09 [1] CRAN (R 4.3.0)
 readxl               * 1.4.3     2023-07-06 [1] CRAN (R 4.3.0)
 remotes                2.5.0     2024-03-17 [1] CRAN (R 4.3.3)
 reshape2             * 1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
 reticulate             1.37.0    2024-05-21 [1] CRAN (R 4.3.3)
 rlang                * 1.1.4     2024-06-04 [1] CRAN (R 4.3.3)
 ROCR                   1.0-11    2020-05-02 [1] CRAN (R 4.3.0)
 RSpectra               0.16-1    2022-04-24 [1] CRAN (R 4.3.0)
 RSQLite                2.3.7     2024-05-27 [1] CRAN (R 4.3.3)
 rstatix                0.7.2     2023-02-01 [1] CRAN (R 4.3.0)
 rstudioapi             0.16.0    2024-03-24 [1] CRAN (R 4.3.3)
 rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.3.0)
 Rtsne                  0.17      2023-12-07 [1] CRAN (R 4.3.0)
 S4Arrays               1.2.0     2023-10-24 [1] Bioconductor
 S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.3)
 ScaledMatrix           1.8.1     2023-05-03 [1] Bioconductor
 scales               * 1.3.0     2023-11-28 [1] CRAN (R 4.3.0)
 scater               * 1.28.0    2023-04-25 [1] Bioconductor
 scattermore            1.2       2023-06-12 [1] CRAN (R 4.3.0)
 SCpubr               * 2.0.2     2023-10-11 [1] CRAN (R 4.3.0)
 sctransform            0.4.1     2023-10-19 [1] CRAN (R 4.3.0)
 scuttle              * 1.10.3    2023-10-10 [1] Bioconductor
 sessioninfo            1.2.2     2021-12-06 [1] CRAN (R 4.3.3)
 Seurat               * 5.0.1     2023-11-17 [1] CRAN (R 4.3.0)
 SeuratObject         * 5.0.1     2023-11-17 [1] CRAN (R 4.3.0)
 SeuratWrappers       * 0.3.1     2023-08-21 [1] Github (satijalab/seurat-wrappers@d28512f)
 SHARP                * 1.1.0     2023-07-19 [1] Github (shibiaowan/SHARP@b1eddac)
 shiny                  1.8.0     2023-11-17 [1] CRAN (R 4.3.0)
 SingleCellExperiment * 1.22.0    2023-04-25 [1] Bioconductor
 sp                   * 2.1-2     2023-11-26 [1] CRAN (R 4.3.0)
 spam                   2.10-0    2023-10-23 [1] CRAN (R 4.3.0)
 sparseMatrixStats      1.12.2    2023-07-02 [1] Bioconductor
 SPARSim              * 0.9.5     2024-01-16 [1] gitlab (sysbiobig/sparsim@19f86cb)
 spatstat.data          3.0-4     2024-01-15 [1] CRAN (R 4.3.0)
 spatstat.explore       3.2-5     2023-10-22 [1] CRAN (R 4.3.0)
 spatstat.geom          3.2-7     2023-10-20 [1] CRAN (R 4.3.0)
 spatstat.random        3.2-2     2023-11-29 [1] CRAN (R 4.3.0)
 spatstat.sparse        3.0-3     2023-10-24 [1] CRAN (R 4.3.0)
 spatstat.utils         3.0-5     2024-06-17 [1] CRAN (R 4.3.3)
 stringi                1.8.4     2024-05-06 [1] CRAN (R 4.3.3)
 stringr                1.5.1     2023-11-14 [1] CRAN (R 4.3.0)
 SummarizedExperiment * 1.30.2    2023-06-06 [1] Bioconductor
 survival               3.7-0     2024-06-05 [1] CRAN (R 4.3.3)
 tensor                 1.5       2012-05-05 [1] CRAN (R 4.3.0)
 terra                  1.7-78    2024-05-22 [1] CRAN (R 4.3.3)
 tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
 tidyr                  1.3.0     2023-01-24 [1] CRAN (R 4.3.0)
 tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
 tidytree             * 0.4.6     2023-12-12 [1] CRAN (R 4.3.3)
 TMExplorer           * 1.12.0    2023-10-26 [1] Bioconductor
 treeio                 1.26.0    2023-10-24 [1] Bioconductor
 utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.0)
 uwot                   0.1.16    2023-06-29 [1] CRAN (R 4.3.0)
 vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.0)
 vipor                  0.4.7     2023-12-18 [1] CRAN (R 4.3.0)
 viridis                0.6.5     2024-01-29 [1] CRAN (R 4.3.3)
 viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
 withr                  3.0.0     2024-01-16 [1] CRAN (R 4.3.3)
 xtable                 1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
 XVector                0.40.0    2023-04-25 [1] Bioconductor
 yulab.utils            0.1.3     2024-01-08 [1] CRAN (R 4.3.0)
 zip                    2.3.1     2024-01-27 [1] CRAN (R 4.3.3)
 zlibbioc               1.48.2    2024-03-13 [1] Bioconductor 3.18 (R 4.3.3)
 zoo                    1.8-12    2023-04-13 [1] CRAN (R 4.3.0)

 [1] /home/yasloudj/R/x86_64-pc-linux-gnu-library/4.3
 [2] /usr/local/lib/R/site-library
 [3] /usr/lib/R/site-library
 [4] /usr/lib/R/library

───────────────────────────────────────────────────────────────────────────
```
