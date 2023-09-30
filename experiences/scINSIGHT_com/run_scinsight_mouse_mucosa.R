setwd("D:/0AALRJ/scRNA-seq/VAE/scDisco/experiences/mouse_mucosa")
library(scINSIGHT)
library(SeuratDisk)
library(Seurat)
library(stringr)
library(aricode)

memory.limit(size=80000)

# ### 1. Running Complete data by scINSIGHT -------------------------------------------------
data_dir = '../../datasets/mouse_mucosa/'
# # prepare data
meta = read.csv(paste0(data_dir, "meta.csv"))
row.names(meta)=meta[['X']]

Convert(paste0(data_dir, "mouse_mucosa.h5ad"), dest="h5seurat",
        assay = "RNA",
        overwrite=F)
seurat_object <- LoadH5Seurat(paste0(data_dir , "mouse_mucosa.h5seurat"),
                              meta.data = FALSE, misc = FALSE)
seurat_object <- UpdateSeuratObject(seurat_object)
colnames(seurat_object@assays[["RNA"]]@data) <- meta[['X']]
seurat_object <- AddMetaData(object = seurat_object, metadata = meta)
seurat_object@meta.data[["batch"]] <- factor(seurat_object@meta.data[["batch"]])

data.list <- SplitObject(seurat_object, split.by = "batch")

batch_conditions <- lapply(data.list, function(seurat_obj) {
  seurat_obj@meta.data[["stimulus"]]
  
})
batch_unique_conditions <- lapply(batch_conditions, unique)
conditions <- unlist(batch_unique_conditions)
conditions <- as.character(conditions)


data.list <- lapply(data.list, function(x) {
  x <- subset(x, subset = nFeature_RNA > 200)
  x <- NormalizeData(x)
  all.genes <- rownames(x)
  x <- ScaleData(x, features = all.genes)
  return(x)
})

mapping.list <- lapply(data.list, function(seurat_obj) {
  data.frame(seurat_obj@meta.data[["batch"]],
             seurat_obj@meta.data[["celltype"]],
             seurat_obj@meta.data[["stimulus"]])
})
mapping = do.call(rbind, mapping.list)
colnames(mapping)<- c("batch", 'celltype', 'stimulus')


matrix_list <- list()
for (i in 1:length(data.list)) {
  matrix <- data.list[[i]]@assays[["RNA"]]@scale.data
  matrix_list[[i]] <- matrix
}
# run scINSIGHT
ptm <- proc.time()
scinsight_object <- create_scINSIGHT(matrix_list, conditions)
scinsight_result <- run_scINSIGHT(scinsight_object,
                                  K = seq(5, 15, 2),
                                  K_j = 2,
                                  LDA = c(0.001, 0.01, 0.1, 1, 10),
                                  thre.niter = 500,
                                  thre.delta = 0.01,
                                  num.cores = 1,
                                  B = 3,
                                  out.dir = NULL,
                                  method = "increase")
time = proc.time()-ptm
print(time)


clust = scinsight_result@clusters
clusts = unlist(clust)
write.csv(clusts, paste0(data_dir, "mouse_mucosa_scINSIGHT_clust.csv"))


combined_matrix <- do.call(rbind, scinsight_result@norm.W_2)
results = data.frame(combined_matrix, mapping)
write.csv(results, paste0(data_dir, "mouse_mucosa_scINSIGHT.csv"))


print(memory.profile())
print(memory.size(max=TRUE))



