setwd("D:/0AALRJ/scRNA-seq/VAE/scDisco/experiences/human_lung")
rm(list=ls())
library(SeuratDisk)
library(Seurat)
library(stringr)
library(aricode)

### 1. Running 10 sampled data ----------------------------------------------------------------
data_dir = '../../datasets/human_lung/sampling/'
seed = list(0, 1, 4, 6, 7, 8, 10, 11, 100, 111)
all_clust <- list()
for (i in seed) {
        data <- read.csv(paste0(data_dir, str_c(i), '_', "lung_raw.csv"), header = T, row.names = 1)
        meta = read.csv(paste0(data_dir, str_c(i), '_', "meta.csv"), header = T, row.names = 1)
        seurat_object <- CreateSeuratObject(counts = t(data),
                                            project = "CreateSeuratObject",
                                            meta.data =meta)
        simulate.list <- SplitObject(seurat_object, split.by = "batch")
        # normalize and identify variable features for each dataset independently
        simulate.list <- lapply(X = simulate.list, FUN = function(x) {
          x <- NormalizeData(x)
          x <- FindVariableFeatures(x)})

        simulate.anchors <- FindIntegrationAnchors(object.list = simulate.list, dims = 1:20)
        simulate.combined <- IntegrateData(anchorset = simulate.anchors, k.weight = 10)

        simulate.combined <- ScaleData(simulate.combined, verbose = FALSE)
        simulate.combined <- RunPCA(simulate.combined, verbose = FALSE)
        xpca = simulate.combined@reductions[["pca"]]@cell.embeddings

        simulate.combined <- FindNeighbors(simulate.combined)
        simulate.combined <- FindClusters(simulate.combined)

        clust = simulate.combined@meta.data[["seurat_clusters"]]
        # all_clust[[str_c(i)]] <- c(clust)
        write.csv(clust, paste0(data_dir, str_c(i), '_', "lung_seurat_clust.csv"))

        results = data.frame(xpca,
                             simulate.combined@meta.data[["celltype"]],
                             simulate.combined@meta.data[["batch"]],
                             simulate.combined@meta.data[["condition"]]
                             )
        write.csv(results, paste0(data_dir, str_c(i), '_', "lung_seurat.csv"))

}

## 2. Running Complete data by seurat -------------------------------------------------
data_dir = '../../datasets/human_lung/'
# # prepare data
meta = read.csv(paste0(data_dir, "meta.csv"))
row.names(meta)=meta[['X']]

Convert(paste0(data_dir, "adata_lung.h5ad"), dest="h5seurat",
        assay = "RNA",
        overwrite=F)
seurat_object <- LoadH5Seurat(paste0(data_dir , "adata_lung.h5seurat"),
                              meta.data = FALSE, misc = FALSE)
seurat_object <- UpdateSeuratObject(seurat_object)
colnames(seurat_object@assays[["RNA"]]@data) <- meta[['X']]
seurat_object <- AddMetaData(object = seurat_object, metadata = meta)
seurat_object@meta.data[["batch"]] <- factor(seurat_object@meta.data[["batch"]])

simulate.list <- SplitObject(seurat_object, split.by = "batch")

# normalize and identify variable features for each dataset independently
simulate.list <- lapply(X = simulate.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)})

# run Seurat
ptm <- proc.time()
simulate.anchors <- FindIntegrationAnchors(object.list = simulate.list, dims = 1:20)
simulate.combined <- IntegrateData(anchorset = simulate.anchors, k.weight=10)
simulate.combined <- ScaleData(simulate.combined, verbose = FALSE)
simulate.combined <- RunPCA(simulate.combined, verbose = FALSE)
xpca = simulate.combined@reductions[["pca"]]@cell.embeddings
time = proc.time()-ptm
print(time)

simulate.combined <- FindNeighbors(simulate.combined)
simulate.combined <- FindClusters(simulate.combined)
clust = simulate.combined@meta.data[["seurat_clusters"]]
write.csv(clust, paste0(data_dir, "lung_seurat_clust.csv"))

results = data.frame(xpca,
                     simulate.combined@meta.data[["celltype"]],
                     simulate.combined@meta.data[["batch"]],
                     simulate.combined@meta.data[["age"]]
                     )
write.csv(results, paste0(data_dir, "lung_seurat.csv"))

print(memory.profile())
print(memory.size(max=TRUE))

