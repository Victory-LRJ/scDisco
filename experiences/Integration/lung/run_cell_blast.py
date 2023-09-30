import os
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import Cell_BLAST as cb
from util.utils import set_seed
from memory_profiler import profile
from time import time

os.environ['PYTHONHASHSEED'] = '0'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
# ## ### conda activate cb -----------------------------------------------------------------
cb.config.N_JOBS = 4
cb.config.RANDOM_SEED = 0

data_dir = '../../datasets/human_lung/'
# # 1. Running 10 sampled data ---------------------------------------------------------
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

    axes = cb.data.find_variable_genes(adata, grouping="batch")
    adata.var["variable_genes"].sum()
    model = cb.directi.fit_DIRECTi(adata, genes=adata.var.query("variable_genes").index.to_numpy(),
                                   batch_effect="batch", latent_dim=10, cat_dim=20)
    adata.obsm['X_latent'] = model.inference(adata)

    corrd = pd.DataFrame(adata.obsm['X_latent'])
    adata_corrd = ad.AnnData(corrd, obs=adata.obs, dtype='float64')
    adata_corrd.obsm['X_latent'] = adata.obsm['X_latent']
    adata_corrd.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_corrd.obs['batch'] = np.array(adata.obs['batch'])
    adata_corrd.obs['condition'] = np.array(adata.obs['condition'])
    adata_corrd.write_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'lung_cell_blast.h5ad')

# # 2. Running Complete data by cell blast -------------------------------------------------------
@profile
def my_func(adata):
    start = time()
    axes = cb.data.find_variable_genes(adata, grouping="batch")
    adata.var["variable_genes"].sum()
    model = cb.directi.fit_DIRECTi(adata, genes=adata.var.query("variable_genes").index.to_numpy(),
                                   batch_effect="batch", latent_dim=10, cat_dim=20)
    adata.obsm['X_latent'] = model.inference(adata)
    end = time()
    print('elapsed{:.2f} seconds'.format(end - start))
    return adata

if __name__ == '__main__':
    set_seed(2023)
    data = sc.read_h5ad(data_dir + 'adata_lung.h5ad')

    adata = my_func(data)
    corrd = pd.DataFrame(adata.obsm['X_latent'])
    adata_corrd = ad.AnnData(corrd, obs=adata.obs, dtype='float64')
    adata_corrd.obsm['X_latent'] = adata.obsm['X_latent']
    adata_corrd.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_corrd.obs['batch'] = np.array(adata.obs['batch'])
    adata_corrd.obs['condition'] = np.array(adata.obs['age'])

    adata_corrd.write_h5ad(data_dir + "lung_cell_blast.h5ad")
