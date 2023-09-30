import scanpy as sc
import anndata as ad
import torch
import numpy as np
import pandas as pd
from util.datasets import data_sample

# # 1. Loading and preprocessing data
data_dir = '../../datasets/human_ductal/'
device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
# ## 1.1 Read datasets from source
adata = sc.read_h5ad(data_dir + 'human_ductal.h5ad')
adata.X = adata.X.todense()
aspect_key = ['disease']

data = pd.DataFrame(adata.X, index=adata.obs_names, dtype='float64')
data.to_csv(data_dir + 'ductal_raw.csv')

df = adata.obs
df.index = adata.obs_names
df.to_csv(data_dir + 'meta.csv')

# # 1.2 Sampling
adata.obs['condition'] = adata.obs['disease']
seed = [0, 1, 4, 6, 7, 8, 10, 11, 100, 111]
for i in seed:
    X, y, b, c = data_sample(adata, adata.obs['celltype'], adata.obs['batch'], adata.obs['condition'], i)
    adata = ad.AnnData(X)
    adata.obsm['X_latent'] = adata.X
    adata.obs['celltype'] = np.array(y)
    adata.obs['batch'] = np.array(b)
    adata.obs['condition'] = np.array(c)

    adata.obs['celltype'] = adata.obs['celltype'].astype('category')
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    adata.obs['condition'] = adata.obs['condition'].astype('category')

    adata.raw = adata

    data = pd.DataFrame(adata.X, index=adata.obs_names, dtype='float64')
    data.to_csv(data_dir + 'sampling/' + str(i) + '_' + 'ductal_raw.csv')

    df = adata.obs
    df.index = adata.obs_names
    df.to_csv(data_dir + 'sampling/' + str(i) + '_' + 'meta.csv')

    np.savez(data_dir + 'sampling/' + str(i) + '_' + 'ductal_raw.npz', X_latent=data, meta=df)
