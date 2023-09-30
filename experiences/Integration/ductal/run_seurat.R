setwd("D:/0AALRJ/scRNA-seq/VAE/scDisco/experiences/human_ductal")
rm(list=ls())
library(SeuratDisk)
library(Seurat)
library(stringr)
library(aricode)
memory.limit(size = 100000)

## 1. Running Complete data by seurat -------------------------------------------------
data_dir = '../../datasets/human_ductal/'
# # prepare data
data <- read.csv(paste0(data_dir, "ductal_raw.csv"), header = T, row.names = 1)
meta = read.csv(paste0(data_dir, "meta.csv"), header = T, row.names = 1)
seurat_object <- CreateSeuratObject(counts = t(data),
                                    project = "CreateSeuratObject",
                                    meta.data =meta)
gc()
simulate.list <- SplitObject(seurat_object, split.by = "batch")

# normalize and identify variable features for each dataset independently
simulate.list <- lapply(X = simulate.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)})
gc()
# run Seurat
ptm <- proc.time()
simulate.anchors <- FindIntegrationAnchors(object.list = simulate.list, dims = 1:20)
gc()
simulate.combined <- IntegrateData(anchorset = simulate.anchors, k.weight=10)
gc()
simulate.combined <- ScaleData(simulate.combined, verbose = FALSE)
simulate.combined <- RunPCA(simulate.combined, verbose = FALSE)
gc()
xpca = simulate.combined@reductions[["pca"]]@cell.embeddings
gc()
time = proc.time()-ptm
print(time)

simulate.combined <- FindNeighbors(simulate.combined)
gc()
simulate.combined <- FindClusters(simulate.combined)
gc()
clust = simulate.combined@meta.data[["seurat_clusters"]]
write.csv(clust, paste0(data_dir, "ductal_seurat_clust.csv"))

results = data.frame(xpca,
                     simulate.combined@meta.data[["celltype"]],
                     simulate.combined@meta.data[["batch"]],
                     simulate.combined@meta.data[["disease"]]
                     )
write.csv(results, paste0(data_dir, "ductal_seurat.csv"))

print(memory.profile())
print(memory.size(max=TRUE))

