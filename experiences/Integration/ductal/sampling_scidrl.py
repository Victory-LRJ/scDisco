import scanpy as sc
import torch
import pandas as pd

# # 1. Loading and preprocessing data
data_dir = 'C://0scDisco/datasets/ductal/'
device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
# ## 1.1 Read datasets from source
# adata = sc.read_h5ad(data_dir + 'human_ductal.h5ad')
# # adata.X = adata.X.todense()
# sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
# sc.pp.log1p(adata)
# sc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor='seurat_v3', inplace=True)
# adata.write_h5ad(data_dir + 'human_ductal_scidrl.h5ad')

adata = sc.read_h5ad(data_dir + 'human_ductal_scidrl.h5ad')
adata.X = adata.X.todense()

data = pd.DataFrame(adata.X, index=adata.obs_names, dtype='float64')
data.to_csv(data_dir + 'ductal_raw_scidrl.csv')

df = adata.obs
df.index = adata.obs_names
df.to_csv(data_dir + 'meta_scidrl.csv')
