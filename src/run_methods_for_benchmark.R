"Run this script to apply individual clustering algorithms on the scRNA-seq dataset.

	2023/10/13 @yanisaspic"

suppressPackageStartupMessages({
  library(scater)
  library(Seurat)
  library(SeuratWrappers)
  library(SHARP)
  library(densitycut)
  library(monocle3)
  library(glue)
  library(ggplot2)
  library(cidr)
  library(scLCA)
  library(SIMLR)
  library(scCCESS)
  library(rjson)
  library(lubridate)
  library(SCpubr)
})
test <- "OK"
source("./src/utils/utils.R")


SEED <- as.integer(commandArgs(trailingOnly = TRUE)[1])
set.seed(SEED)


DATA <- read.csv("./results/tmp/ALL_DATA.csv", row.names = 1, header = TRUE)
params <- fromJSON(file = "./src/vars/params.json")
classification <- read.csv("./results/tmp/classification.csv", row.names = 1, header = TRUE)
N_HVG <- params[["N_HVG"]]


HVG <- get_hvg(DATA, n = N_HVG)
count_data <- DATA[HVG, ]
so <- CreateSeuratObject(count_data)
times <- list()



#############################################
third_party_results <- list() #
#############################################
# DIMENSIONALITY REDUCTION [SEURAT]         #
#   -> get a UMAP                           #
# required for Monocle3 method              #
# npcs=50 | n.neighbors=30                  #
#   -> error if ncols <= npcs | n.neighbors #
#############################################
so <- ScaleData(so, features = rownames(count_data))
so <- RunUMAP(so, features = rownames(count_data), seed.use = SEED)

#######################################
# CLUSTERING [SEURAT]                 #
#   -> community detection            #
# default hyperparameters (Louvain)   #
# #####################################
start <- Sys.time()
so <- RunPCA(so, features = rownames(count_data), seed.use = SEED)
so <- FindNeighbors(so, features = rownames(count_data))
so <- FindClusters(so, random.seed = SEED)
third_party_results$Seurat <- unname(Idents(so))
times$Seurat <- round(Sys.time() - start, 2)

# # Yu et al. (2022): Seurat overestimates the number of clusters
min_clusters <- 2 # min_clusters >1 to run SHARP
max_clusters <- nlevels(so$seurat_clusters)

#############################
# CLUSTERING [MONOCLE3]     #
#   -> community detection  #
# uses the UMAP of Seurat   #
#############################
start <- Sys.time()
cds <- as.cell_data_set(so)
cds <- cluster_cells(cds, random_seed = SEED)
third_party_results$monocle3 <- unname(cds@clusters@listData$UMAP$clusters)
times$monocle3 <- round(Sys.time() - start, 2)

#####################################################
# CLUSTERING [SHARP]                                #
#   -> inter/intra- cluster similarity              #
# hyperparameters for raw count data                #
#   -> SEED+1 to prevent errors specific to SEED==0 #
#####################################################
start <- Sys.time()
third_party_results$SHARP <- factor(
  SHARP(
    scExp = count_data,
    exp.type = "count",
    minN.cluster = min_clusters, maxN.cluster = max_clusters,
    rN.seed = SEED + 1, n.cores = 4,
  )$pred_clusters
)
times$SHARP <- round(Sys.time() - start, 2)



###############################################################################
# DATA TRANSFORMATION [SCATER]                                                #
#   -> get log2_tpm from raw count for:                                       #
#       .scCCESS methods (log2_tpm, see Geddes et al. 2019)                   #
###############################################################################
log2_tpm_data <- log2(calculateTPM(count_data) + 1) #
###############################################################################



#############################
# CLUSTERING [densitycut]   #
#   -> stability metric     #
# default hyperparameters   #
#############################
start <- Sys.time()
set.seed(SEED)
third_party_results$densityCut <- factor(DensityCut(t(log2_tpm_data), show.plot = FALSE)$cluster)
times$densityCut <- round(Sys.time() - start, 2)
#############################################################################
third_party_clusters <- as.data.frame(do.call(cbind, third_party_results)) #
#############################################################################


# #########################################
# other_results <- list() #
# #########################################
# # CLUSTERING [SCCESS-KMEANS]            #
# #   -> stability metric                 #
# # hyperparameters from Yu et al. (2022) #
# #########################################
# start <- Sys.time()
# estimation <- estimate_k(t(log2_tpm_data),
#   krange = min_clusters:max_clusters,
#   epochs = 20,
#   learning_rate = 0.005,
#   ensemble_sizes = 20, # benchmark
#   cluster_func = function(x, centers) {
#     kmeans(x, centers)
#   },
#   seed = SEED
# )
# other_results$scCCESS.Kmeans <- factor(
#   ensemble_cluster(t(log2_tpm_data),
#     epochs = 20,
#     learning_rate = 0.005,
#     ensemble_sizes = 20, # benchmark
#     cluster_func = function(x) {
#       kmeans(x, estimation$ngroups)
#     },
#     seed = SEED
#   )
# )
# times$scCCESS.Kmeans <- round(Sys.time() - start, 2)
# 
# ###################################################
# # CLUSTERING [SCCESS-SIMLR]                       #
# #   -> stability metric                           #
# # hyperparameters from Yu et al. (2022)           #
# #   -> min cells to prevent errors not found yet  #
# ###################################################
# start <- Sys.time()
# estimation <- estimate_k(t(log2_tpm_data),
#   krange = min_clusters:max_clusters,
#   epochs = 20,
#   learning_rate = 0.005,
#   ensemble_sizes = 20, # benchmark
#   cluster_func = function(x, centers) {
#     SIMLR_Large_Scale(t(x), c = centers, kk = 15)
#   },
#   seed = SEED
# )
# other_results$scCCESS.SIMLR <- factor(
#   ensemble_cluster(t(log2_tpm_data),
#     epochs = 20,
#     learning_rate = 0.005,
#     ensemble_sizes = 20, # benchmark
#     cluster_func = function(x) {
#       SIMLR_Large_Scale(t(x), c = estimation$ngroups, kk = 15)
#     },
#     seed = SEED
#   )
# )
# times$scCCESS.SIMLR <- round(Sys.time() - start, 2)
# 
# #########################################
# # CLUSTERING [CIDR]                     #
# #   -> inter/intra- cluster similarity  #
# #########################################
# start <- Sys.time()
# set.seed(SEED)
# scd <- scDataConstructor(as.matrix(count_data))
# scd <- determineDropoutCandidates(scd)
# scd <- wThreshold(scd)
# scd <- scDissim(scd)
# scd <- scPCA(scd, plotPC = FALSE)
# scd <- nPC(scd)
# other_results$CIDR <- factor(scCluster(scd)@clusters)
# times$CIDR <- round(Sys.time() - start, 2)
# 
# #############################################
# # CLUSTERING [scLCA]                        #
# #   -> inter/intra- cluster similarity      #
# # hyperparameters for raw count data        #
# #############################################
# start <- Sys.time()
# set.seed(SEED)
# other_results$scLCA <- factor(unname(myscLCA(DATA, clust.max = max_clusters)[[1]]))
# times$scLCA <- round(Sys.time() - start, 2)
# #################################################################
# other_clusters <- as.data.frame(do.call(cbind, other_results)) #
# #################################################################



#############################################################
# DRAW THE CLUSTERING RESULTS OF:                           #
#   -> the third party methods called in the ensemble       #
#   -> other successful scRNA-seq methods in the benchmark  #
#############################################################
pdf(file = glue("./results/clusters.pdf"))
list.clusters <- list(third_party_clusters, other_clusters)
for (i in 1:length(list.clusters)) {
  proj <- ggplot()
  clusters <- list.clusters[[i]]
  methods <- colnames(clusters)
  for (meth in methods) {
    so@active.ident <- factor(clusters[, meth])
    pi <- do_DimPlot(so) + ggtitle(meth)
    proj <- proj + pi
  }
  proj$patches$plots[[1]] <- NULL
  print(proj)
}
dev.off()



###############################################
# GET THE METRICS OF THE INDEPENDENT METHODS  #
###############################################
independent_clusters <- cbind(third_party_clusters, other_clusters)
rownames(independent_clusters) <- colnames(count_data)
metrics_wrt_method <- get_metrics_wrt_col(
  expression = count_data,
  clusters = independent_clusters,
  y = classification$ground
)
times_in_sec <- lapply(times, FUN = as.numeric, units = "secs")
metrics_wrt_method["time.s", ] <- times_in_sec

write.csv(metrics_wrt_method, "./results/tmp/metrics.csv")
