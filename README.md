# scEVE: a scRNA-seq ensemble clustering algorithm that leverages the extrinsic variability to prevent over-clustering

[![R](https://img.shields.io/badge/R-4.3.3-ffffff.svg)](https://cran.r-project.org/)

[![Linux](https://img.shields.io/badge/GNU_bash-5.1.16-ffffff.svg)](https://www.gnu.org/software/bash/)

[![Python](https://img.shields.io/badge/Python-3.10.12-ffffff.svg)](https://www.python.org/)
[![Code style: black](https://img.shields.io/badge/code_style-black-000000.svg)](https://github.com/psf/black)

This repository includes the code of the scEVE algorithm. 
It also includes the code required to generate the results of the JOBIM2024 paper.
___
*Clustering analyses play a fundamental role in single-cell data science. Hundreds of methods  
have been developed to conduct this analysis, but they all generate different results. Benchmarks and reviews make this issue obvious, and they also show that no single clustering method outperforms all the others. Thus, to address this issue and to generate clustering results robust to the method used, scRNA-seq ensemble clustering algorithms have been developed. They usually tackle this issue by trying to minimize the differences across multiple clustering solutions. In this paper, we propose a novel approach to tackle it. We name ”extrinsic variability” the variations in the clustering solutions that are due to methodological choices. We hypothesize that this extrinsic variability is not an issue to be minimized, but rather an informative signal, and that it can be leveraged to prevent over-clustering. To verify our hypothesis, we have developed scEVE, an algorithm that embraces this approach.*
___
***1. Dependencies***

The dependencies required to run scEVE, and to generate our results, are reported in the directory `setup/config/`.
Both Python3 and R dependencies are included. For each dependency, the version used is also reported.

***2. Datasets***

Run `setup.sh` to download and to set-up the datasets used for our results:
```bash
chmod +x setup.sh
./setup.sh
```
This script generates three datasets, stored in the `data/` directory: `Baron_HumPan.csv`, `Li_HumCRC.csv` and `Darmanis_HumGBM.csv`.

***3. Results***

After completing the two steps above, run `paper.R` interactively (e.g. with Rstudio) to reproduce the results of the JOBIM2024 paper.
* **Fig.2** is generated interactively.
* **Fig.3** is stored in the file `figures/C5.pdf`.
* **Fig.4** uses marker genes stored in the sheet `markers` of the file `records.xlsx`.