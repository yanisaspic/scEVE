"Run this script to conduct a scROB analysis.

	2024/01/10 @yanisaspic"

suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})
source("./src/utils/misc.R")
source("./src/utils/independent_cluster.R")

SEED <- 0
METHODS <- c("monocle3", "Seurat", "SHARP", "densityCut")

DATASET.INIT <- read.csv("./results/tmp/ALL_DATA.csv", row.names=1, header=TRUE)
PARAMS <- fromJSON("./src/vars/params.json")
RECORDS <- init_records(DATASET.INIT, PARAMS)

### get the global 2D-projection ###
####################################
SEUROBJ.INIT <- CreateSeuratObject(DATASET.INIT)
SEUROBJ.INIT <- FindVariableFeatures(SEUROBJ.INIT, 
                                     nfeatures=PARAMS[["n_HVGs"]])
SEUROBJ.INIT <- NormalizeData(SEUROBJ.INIT)
SEUROBJ.INIT <- ScaleData(SEUROBJ.INIT, features=VariableFeatures(SEUROBJ.INIT))
SEUROBJ.INIT <- RunUMAP(SEUROBJ.INIT, 
                        features=VariableFeatures(SEUROBJ.INIT),
                        seed.use=SEED)

### main loop ###
#################
STEP <- NA
POPULATION <- get_undug_population(RECORDS)

while (!is.na(POPULATION)) {
  RECORDS[["meta"]][POPULATION, "to_dig"] <- 0
  
  source("./src/trim_data.R")
  if(STEP=="independent_cluster"){}
  if(STEP=="seeds"){source("./src/get_seeds.R")}
  if(STEP=="genes"){source("./src/get_genes.R")}
  merge_pdfs(POPULATION)
  
  if (length(SEEDS)>0) {
    RECORDS <- update_records(RECORDS, SEEDS, POPULATION)
    write.xlsx(RECORDS, "./results/records.xlsx", rowNames=TRUE)
    SEEDS <- list()
  }
  
  POPULATION <- get_undug_population(RECORDS)
}

f <- apply(X=RECORDS$markers, MARGIN=1, FUN=sum)
RECORDS$markers <- RECORDS$markers[as.logical(f),]
write.xlsx(RECORDS, "./results/records.xlsx", rowNames=TRUE)
