import os
import desc
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from time import time
from util.utils import set_seed
from memory_profiler import profile

os.environ['PYTHONHASHSEED'] = '0'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

def getdims(x=(10000,200)):
    """
    This function will give the suggested nodes for each encoder layer
    return the dims for network
    """
    assert len(x) == 2
    n_sample = x[0]
    if n_sample > 20000:# may be need complex network
        dims = [x[-1], 128, 32]
        tol = 0.001
    elif n_sample > 10000:#10000
        dims = [x[-1], 64, 32]
        tol = 0.001
    elif n_sample > 5000: #5000
        dims = [x[-1], 32, 16] #16
        tol = 0.001
    elif n_sample > 2000:
        dims = [x[-1], 128]
        tol = 0.005
    elif n_sample > 500:
        dims = [x[-1], 64]
        tol = 0.005
    else:
        dims = [x[-1], 16]
        tol = 0.01
    return dims, tol


data_dir = '../../datasets/human_lung/'
# 1. Running 10 sampled data ---------------------------------------------------------
seed = [0, 1, 4, 6, 7, 8, 10, 11, 100, 111]
for i in seed:
    set_seed(2023)
    d = np.load(data_dir + 'sampling/' + str(i) + '_' + 'lung_raw.npz', allow_pickle=True)

    adata = ad.AnnData(d['X_latent'])
    adata.obs['celltype'] = np.array(d['celltype'])
    adata.obs['batch'] = np.array(d['batch'])
    adata.obs['condition'] = np.array(d['condition'])

    adata.obs['celltype'] = adata.obs['celltype'].astype('category')
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    adata.obs['condition'] = adata.obs['condition'].astype('category')

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.obs['n_counts'] = np.matrix(adata.X.sum(axis=1)).A1
    # adata = adata[adata.obs['n_genes'] < 2500, :].copy()
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000, subset=True, inplace=True)
    adata = adata[:, adata.var['highly_variable']].copy()
    adata = desc.scale_bygroup(adata, groupby="batch")
    adata_out = adata.copy()

    dims, tol = getdims(adata.shape)
    adata_out = desc.train(adata_out,
                           dims=dims,
                           tol=tol,
                           n_neighbors=10,
                           batch_size=256,
                           louvain_resolution=[0.8],  #### 1.0
                           # not necessarily a list, you can only set one value, like, louvain_resolution=1.0
                           save_dir=str(data_dir + 'sampling/' + str(i) + '_desc'),
                           do_umap=True,
                           learning_rate=200,  # the parameter of tsne
                           use_GPU=False,
                           num_Cores=1,  # for reproducible, only use 1 cpu
                           num_Cores_tsne=4,
                           save_encoder_weights=False,
                           use_ae_weights=False,
                           )  # if do_uamp is False, it will don't compute umap coordiate

    desc_umap = adata_out.obsm['X_umap0.8']

    corrd = pd.DataFrame(desc_umap)
    adata_corrd = ad.AnnData(corrd, obs=adata.obs, dtype='float64')
    adata_corrd.obsm['X_latent'] = adata_out.obsm['X_umap0.8']
    adata_corrd.obs['clustering'] = adata_out.obs['desc_0.8']
    adata_corrd.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_corrd.obs['batch'] = np.array(adata.obs['batch'])
    adata_corrd.obs['condition'] = np.array(adata.obs['condition'])

    adata_corrd.write_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'lung_desc.h5ad')

# # 2. Running Complete data by combat ------------------------------------------------------------
@profile
def my_func(adata):
    set_seed(2023)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.obs['n_counts'] = np.matrix(adata.X.sum(axis=1)).A1
    # adata = adata[adata.obs['n_genes'] < 2500, :].copy()
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000, subset=True, inplace=True)
    adata = adata[:, adata.var['highly_variable']].copy()
    adata = desc.scale_bygroup(adata, groupby="batch")
    adata_out = adata.copy()

    start = time()
    dims, tol = getdims(adata.shape)
    adata_out = desc.train(adata_out,
                           dims=dims,  ### dims=[adata.shape[1], 128, 32], 5000-10000
                           tol=tol,  ### tol=0.001, 5000-10000
                           n_neighbors=10,
                           batch_size=256,
                           louvain_resolution=[0.8],  #### 1.0
                           # not necessarily a list, you can only set one value, like, louvain_resolution=1.0
                           save_dir=str(data_dir + '_desc'),
                           do_umap=True,
                           learning_rate=200,  # the parameter of tsne
                           use_GPU=False,
                           num_Cores=1,  # for reproducible, only use 1 cpu
                           num_Cores_tsne=4,
                           save_encoder_weights=False,
                           use_ae_weights=False,
                           )  # if do_uamp is False, it will don't compute umap coordiate
    end = time()
    print('elapsed{:.2f} seconds'.format(end - start))
    desc_umap = adata_out.obsm['X_umap0.8']
    return desc_umap, adata_out

if __name__ == '__main__':
    set_seed(2023)
    adata = sc.read_h5ad(data_dir + 'adata_lung.h5ad')

    desc_umap, adata_out = my_func(adata)

    corrd = pd.DataFrame(desc_umap)
    adata_corrd = ad.AnnData(corrd, obs=adata.obs, dtype='float64')
    adata_corrd.obsm['X_latent'] = adata_out.obsm['X_umap0.8']
    adata_corrd.obs['clustering'] = adata_out.obs['desc_0.8']
    adata_corrd.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_corrd.obs['batch'] = np.array(adata.obs['batch'])
    adata_corrd.obs['condition'] = np.array(adata.obs['age'])

    adata_corrd.write_h5ad(data_dir + "lung_desc.h5ad")

