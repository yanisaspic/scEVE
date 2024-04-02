# scEVE: a scRNA-seq ensemble clustering algorithm that leverages the extrinsic variability to prevent over-clustering

[![R](https://img.shields.io/badge/R-4.3.3-ffffff.svg)](https://cran.r-project.org/)

[![Linux](https://img.shields.io/badge/GNU_bash-5.1.16-ffffff.svg)](https://www.gnu.org/software/bash/)

[![Python](https://img.shields.io/badge/Python-3.10.12-ffffff.svg)](https://www.python.org/)
[![Code style: black](https://img.shields.io/badge/code_style-black-000000.svg)](https://github.com/psf/black)

___
*Clustering analyses play a fundamental role in single-cell data science. Hundreds of methods  
have been developed to conduct this analysis, but they all generate different results. Benchmarks and reviews make this issue obvious, and they also show that no single clustering method outperforms all the others. Thus, to address this issue and to generate clustering results robust to the method used, scRNA-seq ensemble clustering algorithms have been developed. They usually tackle this issue by trying to minimize the differences across multiple clustering solutions. In this paper, we propose a novel approach to tackle it. We name ”extrinsic variability” the variations in the clustering solutions that are due to methodological choices. We hypothesize that this extrinsic variability is not an issue to be minimized, but rather an informative signal, and that it can be leveraged to prevent over-clustering. To verify our hypothesis, we have developed scEVE, an algorithm that embraces this approach. In this paper, we present the algorithm of scEVE. We apply it on a human glioblastoma scRNA-seq dataset, and we compare its performance to three state of the art ensemble clustering algorithms, on two different scRNA-seq datasets. We start by presenting scEVE, and how it effectively prevents over-clustering. Then, we showcase it on the public glioblastoma dataset, and we reveal the existence of a sub-cluster of cancer cells, that we characterize biologically. Finally, we show that scEVE performs well compared to existing methods, on top of addressing two main challenges in scRNA-seq clustering. Overall, our work shows that the extrinsic variability can be informative and we present scEVE, an algorithm that leverages it to generate a multi-resolution clustering with explicit consensus values.*
___
Generate scRNA-seq clustering results robust to the choice of methodology by integrating a pool of orthogonal methods iteratively.

The pool includes :
1. Seurat
2. monocle3
3. scCCESS-SIMLR
4. scCCESS-k-means

This package includes :
- interacting Python & R scripts to generate the iterative ensemblist clustering.

## References

1.
2.
3.
4.
5.