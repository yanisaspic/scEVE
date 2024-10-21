# Capitalizing on the differences of prediction between multiple scRNA-seq clustering methods.

[![R](https://img.shields.io/badge/R-4.3.1-blue.svg)](https://cran.r-project.org/)
[![Linux](https://img.shields.io/badge/Ubuntu-20.04.6-white.svg)](https://ubuntu.com/)
[![bash](https://img.shields.io/badge/bash-5.0.17-white.svg)](https://www.gnu.org/software/bash/)
[![Python](https://img.shields.io/badge/Python-4.3.1-white.svg)](https://www.python.org/) 

![The figure is composed of three parts: a multi-resolution cluster tree (left), a cell type plot (middle), and an abundance plot (right).\\
The multi-resolution cluster tree (left) represents the cell populations predicted by scEVE as nodes. Edges between nodes connect sub-populations to their parents, and their weight indicates the robustness of the sub-populations. The nodes representing the original cell pool and the leftover clusters are black. The cell type plot (middle) associates the highest-resolution populations to bars. If a population is homogeneous (\textit{i.e.} composed of a unique cell type), its bar is monochrome. Otherwise, its bar is polychrome. The abundance plot (right) also associates highest-resolution populations to bars. The bars are log-scaled, and they indicate the size of the cell populations.](/plots/Darmanis_HumGBM.png?raw=true "Clustering results of scEVE")

___
*Single-cell RNA-sequencing measures individual cell transcriptomes in a tissue. In the past decade, that technology has motivated the development of hundreds of clustering methods. These methods attempt to group cells into populations by leveraging the similarity of their transcriptomes. Because each method relies on specific hypotheses, their predictions can vary drastically. To address that issue, ensemble algorithms detect cell populations by integrating multiple clustering methods, and minimizing the differences of their predictions. While that approach is sensible, it has yet to address some conceptual challenges in single-cell data science; namely, ensemble algorithms have yet to generate clustering results with uncertainty values and multiple resolutions. In this work, we present an original approach to ensemble clustering that addresses these challenges; by describing the differences between clustering results, instead of minimizing them. We present our algorithm scEVE, and we evaluate it on 15 experimental and up to 600 synthetic datasets. Our results reveal that scEVE performs decently, and is the first scRNA-seq ensemble algorithm to address both challenges. They also highlight how biological downstream analyses benefit from addressing these challenges. Overall, we expect that our work will provide an alternative direction for developing single-cell ensemble clustering algorithms.*
___
The repository includes the code of the scEVE algorithm, and the codes used to generate the datasets and the results in the related paper "**scEVE: a single-cell RNA-seq ensemble clustering algorithm that capitalizes on the differences of predictions between multiple clustering methods.**" 

# Architecture
The repository is composed of 7 directories:
- `./data/` some of the experimental datasets are stored here.
- `./etc/` there are scripts to run before conducing the clustering analyses.
  - `./etc/config/` to install the dependencies.
  - `./etc/source/` to download and format the experimental datasets.
- `./null/` empty, but used to parallelize multiple scEVE clustering analyses.
- `./plots/` the figures used in the paper are generated here.
- `./results/` the results of each clustering analysis.
    - `./results/benchmark/` the clustering and computational performances.
    - `./results/figures/` the plots generated during the clustering analysis.
    - `./results/records/` the files summarizing a clustering analysis.
    - `./results/similarity` the pairwise similarity between the clustering results.
- `./run/` the scripts stored here are ran to produce our results.
  - `./run/array/` this sub-directory includes SLURM scripts to parallelize the dataset analyses on computer clusters.
- `./src/` the functions imported by `scEVE.R` and the scripts in `./run/`.
  - `./src/scEVE` functions used in the scEVE algorithm.
  - `./src/paper` functions used exclusively to generate our results.

### Regarding the scripts in `./run/` 
- `./run/setup.sh` is ran to download and format the experimental datasets.
- `./run/get_real.R` is ran to analyze one of 15 experimental datasets.
- `./run/get_synthetic.R` is ran to analyze one of 600 synthetic datasets.
- `./run/get_similarity.R` is ran to measure the similarity between clustering results on 150 specific datasets.

# Producing our results
Because of the large number of datasets, we recommend the use of a computer cluster to produce our results.

## 1. Installing the dependencies.
The dependencies and their respective versions are summarized at the end of the README.
They can be installed automatically by running the script `install_requirements.R`.

```bash
Rscript ./etc/config/install_requirements.R
```
Note that the versions of the dependencies installed that way may differ from the versions in the summary.

## 2. Formatting the datasets.
The experimental datasets used in the paper are publicly available.
Some of them must be downloaded and formatted by running the script `setup.sh`.

```bash
chmod +x ./run/setup.sh
./run/setup.sh
```

## 3. Conducting the analyses.
If a computer cluster is used, the SLURM scripts stored in `./run/array/` can be used directly.

```bash
sbatch ./run/array/real.sbatch
sbatch ./run/array/synthetic.sbatch
sbatch ./run/array/similarity.sbatch
```

Note that hidden SLURM scripts are also available in order to analyze a specific dataset.
The scripts expect an argument corresponding to the dataset of interest, _i.e._ a label for experimental datasets, or an integer **n** for synthetic datasets. 

```bash
sbatch ./run/array/real.sbatch Darmanis_HumGBM
sbatch ./run/array/synthetic.sbatch 123
```
**n** ranges from 1 to 600 for `./run/array/synthetic.sbatch/`, and 1 to 150 for `./run/array/similarity.sbatch/`.

## 4. Generating our results.
The script `./run/draw.R` is ran to generate the figures and the contents of our tables. The results are directly available in the `./plots/` directory and in the terminal.

```bash
Rscript ./run/draw.R
```

## Summary of the dependencies
The summary is generated with `session_info()` from the R package `sessioninfo` (see https://github.com/r-lib/sessioninfo).
Required dependencies are indicated by an asterisk (*).

```
─ Session info ───────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.3.1 (2023-06-16)
 os       Ubuntu 20.04.6 LTS
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       Europe/Paris
 date     2024-09-27
 pandoc   3.2.1 @ /shared/ifbstor1/software/miniconda/envs/r-4.3.1/bin/pandoc

─ Packages ───────────────────────────────────────────────────────────────────
 package              * version    date (UTC) lib source
 abind                  1.4-8      2024-09-12 [1] CRAN (R 4.3.1)
 ape                    5.8        2024-04-11 [2] CRAN (R 4.3.1)
 aplot                  0.2.3      2024-06-17 [2] CRAN (R 4.3.1)
 aricode              * 1.0.3      2023-10-20 [2] CRAN (R 4.3.1)
 arules               * 1.7-7      2023-11-29 [1] CRAN (R 4.3.1)
 askpass                1.2.0      2023-09-03 [2] CRAN (R 4.3.1)
 beachmat               2.18.1     2024-02-14 [2] Bioconductor 3.18 (R 4.3.1)
 beeswarm               0.4.0      2021-06-01 [2] CRAN (R 4.3.1)
 Biobase              * 2.62.0     2024-03-21 [2] bioc_xgit (@8201fbb)
 BiocFileCache        * 2.10.2     2024-03-27 [2] Bioconductor 3.18 (R 4.3.1)
 BiocGenerics         * 0.48.1     2023-11-01 [2] Bioconductor
 BiocManager            1.30.25    2024-08-28 [1] CRAN (R 4.3.1)
 BiocNeighbors          1.20.2     2024-01-07 [2] Bioconductor 3.18 (R 4.3.1)
 BiocParallel           1.36.0     2023-10-24 [2] Bioconductor
 BiocSingular           1.18.0     2023-10-24 [2] Bioconductor
 bit                    4.5.0      2024-09-20 [1] CRAN (R 4.3.1)
 bit64                  4.5.2      2024-09-22 [1] CRAN (R 4.3.1)
 bitops                 1.0-8      2024-07-29 [1] CRAN (R 4.3.1)
 blob                   1.2.4      2023-03-17 [2] CRAN (R 4.3.1)
 boot                   1.3-30     2024-02-26 [2] CRAN (R 4.3.1)
 cachem                 1.1.0      2024-05-16 [2] CRAN (R 4.3.1)
 caTools                1.18.3     2024-09-04 [1] CRAN (R 4.3.1)
 cellranger             1.1.0      2016-07-27 [2] CRAN (R 4.3.1)
 cli                    3.6.3      2024-06-21 [2] CRAN (R 4.3.1)
 clues                  0.6.2.2    2019-12-03 [1] url (https://cran.r-project.org/src/contrib/Archive/clues/clues_0.6.2.2.tar.gz)
 cluster                2.1.6      2023-12-01 [2] CRAN (R 4.3.1)
 clusterCrit            1.3.0      2023-11-23 [1] CRAN (R 4.3.1)
 codetools              0.2-20     2024-03-31 [2] CRAN (R 4.3.1)
 colorspace             2.1-1      2024-07-26 [1] CRAN (R 4.3.1)
 cowplot                1.1.3      2024-01-22 [2] CRAN (R 4.3.1)
 crayon                 1.5.3      2024-06-20 [2] CRAN (R 4.3.1)
 curl                   5.2.3      2024-09-20 [1] CRAN (R 4.3.1)
 data.table             1.16.0     2024-08-27 [1] CRAN (R 4.3.1)
 DBI                    1.2.3      2024-06-02 [2] CRAN (R 4.3.1)
 dbplyr               * 2.5.0      2024-03-19 [2] CRAN (R 4.3.1)
 DelayedArray           0.28.0     2023-10-24 [2] Bioconductor
 DelayedMatrixStats     1.24.0     2023-10-24 [2] Bioconductor
 deldir                 2.0-4      2024-02-28 [2] CRAN (R 4.3.1)
 densitycut           * 0.0.1      2024-09-20 [1] bitbucket (jerry00/densitycut_dev@3878b5f)
 digest                 0.6.37     2024-08-19 [1] CRAN (R 4.3.1)
 doParallel             1.0.17     2022-02-07 [2] CRAN (R 4.3.1)
 dotCall64              1.1-1      2023-11-28 [2] CRAN (R 4.3.1)
 dplyr                * 1.1.4      2023-11-17 [2] CRAN (R 4.3.1)
 egg                  * 0.4.5      2019-07-13 [2] CRAN (R 4.3.1)
 fansi                  1.0.6      2023-12-08 [2] CRAN (R 4.3.1)
 farver                 2.1.2      2024-05-13 [2] CRAN (R 4.3.1)
 fastDummies            1.7.4      2024-08-16 [1] CRAN (R 4.3.1)
 fastmap                1.2.0      2024-05-15 [2] CRAN (R 4.3.1)
 filelock               1.0.3      2023-12-11 [2] CRAN (R 4.3.1)
 fitdistrplus           1.2-1      2024-07-12 [1] CRAN (R 4.3.1)
 flashClust             1.01-2     2012-08-21 [2] CRAN (R 4.3.1)
 foreach                1.5.2      2022-02-02 [2] CRAN (R 4.3.1)
 fs                     1.6.4      2024-04-25 [2] CRAN (R 4.3.1)
 future                 1.34.0     2024-07-29 [1] CRAN (R 4.3.1)
 future.apply           1.11.2     2024-03-28 [2] CRAN (R 4.3.1)
 generics               0.1.3      2022-07-05 [2] CRAN (R 4.3.1)
 GenomeInfoDb         * 1.38.8     2024-03-15 [2] Bioconductor 3.18 (R 4.3.1)
 GenomeInfoDbData       1.2.11     2023-11-03 [2] Bioconductor
 GenomicRanges        * 1.54.1     2023-10-29 [2] Bioconductor
 ggbeeswarm             0.7.2      2023-04-29 [2] CRAN (R 4.3.1)
 ggfun                  0.1.6      2024-08-28 [1] CRAN (R 4.3.1)
 ggh4x                * 0.2.8      2024-01-23 [1] CRAN (R 4.3.1)
 ggplot2              * 3.5.1      2024-04-23 [2] CRAN (R 4.3.1)
 ggplotify            * 0.1.2      2023-08-09 [2] CRAN (R 4.3.1)
 ggrepel                0.9.6      2024-09-07 [1] CRAN (R 4.3.1)
 ggridges               0.5.6      2024-01-23 [2] CRAN (R 4.3.1)
 ggtree               * 3.10.1     2024-02-25 [2] Bioconductor 3.18 (R 4.3.1)
 ggVennDiagram        * 1.5.2      2024-02-20 [2] CRAN (R 4.3.1)
 globals                0.16.3     2024-03-08 [2] CRAN (R 4.3.1)
 glue                 * 1.7.0      2024-01-09 [2] CRAN (R 4.3.1)
 goftest                1.2-3      2021-10-07 [2] CRAN (R 4.3.1)
 gplots                 3.1.3.1    2024-02-02 [2] CRAN (R 4.3.1)
 gridExtra            * 2.3        2017-09-09 [2] CRAN (R 4.3.1)
 gridGraphics           0.5-1      2020-12-13 [2] CRAN (R 4.3.1)
 gtable                 0.3.5      2024-04-22 [2] CRAN (R 4.3.1)
 gtools                 3.9.5      2023-11-20 [2] CRAN (R 4.3.1)
 htmltools              0.5.8.1    2024-04-04 [2] CRAN (R 4.3.1)
 htmlwidgets            1.6.4      2023-12-06 [2] CRAN (R 4.3.1)
 httpuv                 1.6.15     2024-03-26 [2] CRAN (R 4.3.1)
 httr                   1.4.7      2023-08-15 [2] CRAN (R 4.3.1)
 ica                    1.0-3      2022-07-08 [2] CRAN (R 4.3.1)
 igraph               * 2.0.3      2024-03-13 [2] CRAN (R 4.3.1)
 IRanges              * 2.36.0     2023-10-24 [2] Bioconductor
 irlba                  2.3.5.1    2022-10-03 [2] CRAN (R 4.3.1)
 iterators              1.0.14     2022-02-05 [2] CRAN (R 4.3.1)
 jsonlite               1.8.9      2024-09-20 [1] CRAN (R 4.3.1)
 KernSmooth             2.23-24    2024-05-17 [2] CRAN (R 4.3.1)
 later                  1.3.2      2023-12-06 [2] CRAN (R 4.3.1)
 lattice                0.22-6     2024-03-20 [2] CRAN (R 4.3.1)
 lazyeval               0.2.2      2019-03-15 [2] CRAN (R 4.3.1)
 leiden                 0.4.3.1    2023-11-17 [2] CRAN (R 4.3.1)
 lifecycle              1.0.4      2023-11-07 [2] CRAN (R 4.3.1)
 listenv                0.9.1      2024-01-29 [2] CRAN (R 4.3.1)
 lme4                   1.1-35.5   2024-07-03 [2] CRAN (R 4.3.1)
 lmtest                 0.9-40     2022-03-21 [2] CRAN (R 4.3.1)
 magrittr               2.0.3      2022-03-30 [2] CRAN (R 4.3.1)
 MASS                   7.3-60.0.1 2024-01-13 [2] CRAN (R 4.3.1)
 Matrix               * 1.6-5      2024-01-11 [2] CRAN (R 4.3.1)
 MatrixGenerics       * 1.14.0     2023-10-24 [2] Bioconductor
 matrixStats          * 1.4.1      2024-09-08 [1] CRAN (R 4.3.1)
 memoise                2.0.1      2021-11-26 [2] CRAN (R 4.3.1)
 mime                   0.12       2021-09-28 [2] CRAN (R 4.3.1)
 miniUI                 0.1.1.1    2018-05-18 [2] CRAN (R 4.3.1)
 minqa                  1.2.8      2024-08-17 [1] CRAN (R 4.3.1)
 monocle3             * 1.3.7      2024-09-27 [1] Github (cole-trapnell-lab/monocle3@98402ed)
 munsell                0.5.1      2024-04-01 [2] CRAN (R 4.3.1)
 mvtnorm                1.3-1      2024-09-03 [1] CRAN (R 4.3.1)
 nlme                   3.1-165    2024-06-06 [2] CRAN (R 4.3.1)
 nloptr                 2.1.1      2024-06-25 [2] CRAN (R 4.3.1)
 openxlsx             * 4.2.5.2    2023-02-06 [1] CRAN (R 4.3.1)
 parallelly             1.38.0     2024-07-27 [1] CRAN (R 4.3.1)
 patchwork            * 1.3.0      2024-09-16 [1] CRAN (R 4.3.1)
 pbapply                1.7-2      2023-06-27 [2] CRAN (R 4.3.1)
 pheatmap               1.0.12     2019-01-04 [2] CRAN (R 4.3.1)
 pillar                 1.9.0      2023-03-22 [2] CRAN (R 4.3.1)
 pkgconfig              2.0.3      2019-09-22 [2] CRAN (R 4.3.1)
 plotly                 4.10.4     2024-01-13 [2] CRAN (R 4.3.1)
 plyr                   1.8.9      2023-10-02 [2] CRAN (R 4.3.1)
 png                    0.1-8      2022-11-29 [2] CRAN (R 4.3.1)
 polyclip               1.10-7     2024-07-23 [1] CRAN (R 4.3.1)
 progressr              0.14.0     2023-08-10 [2] CRAN (R 4.3.1)
 promises               1.3.0      2024-04-05 [2] CRAN (R 4.3.1)
 purrr                  1.0.2      2023-08-10 [2] CRAN (R 4.3.1)
 qpdf                 * 1.3.3      2024-03-25 [1] CRAN (R 4.3.1)
 R.methodsS3            1.8.2      2022-06-13 [2] CRAN (R 4.3.1)
 R.oo                   1.26.0     2024-01-24 [2] CRAN (R 4.3.1)
 R.utils                2.12.3     2023-11-18 [2] CRAN (R 4.3.1)
 R6                     2.5.1      2021-08-19 [2] CRAN (R 4.3.1)
 RANN                   2.6.2      2024-08-25 [1] CRAN (R 4.3.1)
 RColorBrewer         * 1.1-3      2022-04-03 [2] CRAN (R 4.3.1)
 Rcpp                   1.0.13     2024-07-17 [1] CRAN (R 4.3.1)
 RcppAnnoy              0.0.22     2024-01-23 [2] CRAN (R 4.3.1)
 RcppHNSW               0.6.0      2024-02-04 [2] CRAN (R 4.3.1)
 RCurl                  1.98-1.16  2024-07-11 [1] CRAN (R 4.3.1)
 readxl               * 1.4.3      2023-07-06 [2] CRAN (R 4.3.1)
 remotes                2.5.0      2024-03-17 [2] CRAN (R 4.3.1)
 reshape2             * 1.4.4      2020-04-09 [2] CRAN (R 4.3.1)
 reticulate             1.39.0     2024-09-05 [1] CRAN (R 4.3.1)
 rlang                * 1.1.4      2024-06-04 [2] CRAN (R 4.3.1)
 ROCR                   1.0-11     2020-05-02 [2] CRAN (R 4.3.1)
 RSpectra               0.16-2     2024-07-18 [1] CRAN (R 4.3.1)
 RSQLite                2.3.7      2024-05-27 [2] CRAN (R 4.3.1)
 rsvd                   1.0.5      2021-04-16 [2] CRAN (R 4.3.1)
 Rtsne                  0.17       2023-12-07 [2] CRAN (R 4.3.1)
 S4Arrays               1.2.1      2024-03-04 [2] Bioconductor 3.18 (R 4.3.1)
 S4Vectors            * 0.40.2     2023-11-23 [2] Bioconductor 3.18 (R 4.3.1)
 ScaledMatrix           1.10.0     2023-10-24 [2] Bioconductor
 scales               * 1.3.0      2023-11-28 [2] CRAN (R 4.3.1)
 scater               * 1.30.1     2023-11-16 [2] Bioconductor
 scattermore            1.2        2023-06-12 [2] CRAN (R 4.3.1)
 SCpubr               * 2.0.2      2023-10-11 [1] CRAN (R 4.3.1)
 sctransform            0.4.1      2023-10-19 [2] CRAN (R 4.3.1)
 scuttle              * 1.12.0     2023-10-24 [2] Bioconductor
 sessioninfo          * 1.2.2      2021-12-06 [2] CRAN (R 4.3.1)
 Seurat               * 5.1.0      2024-05-10 [2] CRAN (R 4.3.1)
 SeuratObject         * 5.0.2      2024-05-08 [2] CRAN (R 4.3.1)
 SeuratWrappers       * 0.3.4      2024-09-27 [1] Github (satijalab/seurat-wrappers@28f074c)
 SHARP                * 1.1.0      2024-09-27 [1] Github (shibiaowan/SHARP@24a9866)
 shiny                  1.9.1      2024-08-01 [1] CRAN (R 4.3.1)
 SingleCellExperiment * 1.24.0     2023-10-24 [2] Bioconductor
 sp                   * 2.1-4      2024-04-30 [2] CRAN (R 4.3.1)
 spam                   2.10-0     2023-10-23 [2] CRAN (R 4.3.1)
 SparseArray            1.2.4      2024-02-11 [2] Bioconductor 3.18 (R 4.3.1)
 sparseMatrixStats      1.14.0     2023-10-24 [2] Bioconductor
 SPARSim              * 0.9.5      2024-09-23 [1] gitlab (sysbiobig/sparsim@19f86cb)
 spatstat.data          3.1-2      2024-06-21 [2] CRAN (R 4.3.1)
 spatstat.explore       3.3-2      2024-08-21 [1] CRAN (R 4.3.1)
 spatstat.geom          3.3-3      2024-09-18 [1] CRAN (R 4.3.1)
 spatstat.random        3.3-2      2024-09-18 [1] CRAN (R 4.3.1)
 spatstat.sparse        3.1-0      2024-06-21 [2] CRAN (R 4.3.1)
 spatstat.univar        3.0-1      2024-09-05 [1] CRAN (R 4.3.1)
 spatstat.utils         3.1-0      2024-08-17 [1] CRAN (R 4.3.1)
 stringi                1.8.4      2024-05-06 [2] CRAN (R 4.3.1)
 stringr                1.5.1      2023-11-14 [2] CRAN (R 4.3.1)
 SummarizedExperiment * 1.32.0     2023-10-24 [2] Bioconductor
 survival               3.7-0      2024-06-05 [2] CRAN (R 4.3.1)
 tensor                 1.5        2012-05-05 [2] CRAN (R 4.3.1)
 tibble                 3.2.1      2023-03-20 [2] CRAN (R 4.3.1)
 tidyr                  1.3.1      2024-01-24 [2] CRAN (R 4.3.1)
 tidyselect             1.2.1      2024-03-11 [2] CRAN (R 4.3.1)
 tidytree             * 0.4.6      2023-12-12 [2] CRAN (R 4.3.1)
 TMExplorer           * 1.12.0     2023-10-26 [1] Bioconductor
 treeio                 1.26.0     2023-10-24 [2] Bioconductor
 utf8                   1.2.4      2023-10-22 [2] CRAN (R 4.3.1)
 uwot                   0.2.2      2024-04-21 [2] CRAN (R 4.3.1)
 vctrs                  0.6.5      2023-12-01 [2] CRAN (R 4.3.1)
 vipor                  0.4.7      2023-12-18 [2] CRAN (R 4.3.1)
 viridis                0.6.5      2024-01-29 [2] CRAN (R 4.3.1)
 viridisLite            0.4.2      2023-05-02 [2] CRAN (R 4.3.1)
 withr                  3.0.1      2024-07-31 [1] CRAN (R 4.3.1)
 xtable                 1.8-4      2019-04-21 [2] CRAN (R 4.3.1)
 XVector                0.42.0     2023-10-24 [2] Bioconductor
 yulab.utils            0.1.7      2024-08-26 [1] CRAN (R 4.3.1)
 zip                    2.3.1      2024-01-27 [2] CRAN (R 4.3.1)
 zlibbioc               1.48.2     2024-03-13 [2] Bioconductor 3.18 (R 4.3.1)
 zoo                    1.8-12     2023-04-13 [2] CRAN (R 4.3.1)

 [1] /shared/ifbstor1/home/yasloudj/scEVE.paper
 [2] /shared/ifbstor1/software/miniconda/envs/r-4.3.1/lib/R/library
──────────────────────────────────────────────────────────────────────────────
```
