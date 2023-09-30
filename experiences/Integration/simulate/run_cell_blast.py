import os
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import Cell_BLAST as cb
from time import time
from memory_profiler import profile

os.environ['PYTHONHASHSEED'] = '0'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
# ## ### conda activate cb -----------------------------------------------------------------
cb.config.N_JOBS = 4
cb.config.RANDOM_SEED = 0

data_dir = '../../datasets/simulate/'
# 1. Running 10 sampled data ---------------------------------------------------------
seed = [0, 1, 4, 6, 7, 8, 10, 11, 100, 111]
for i in seed:
    adata = sc.read_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'simulate_raw.h5ad')

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

    adata_corrd.write_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'simulate_cell_blast.h5ad')

# # 2. Running Complete data
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
    data_dir_ = '../../datasets/simulate/'
    data = sc.read_h5ad(data_dir + 'sim_count.h5ad')

    adata = my_func(data)
    corrd = pd.DataFrame(adata.obsm['X_latent'])
    adata_corrd = ad.AnnData(corrd, obs=adata.obs, dtype='float64')
    adata_corrd.obsm['X_latent'] = adata.obsm['X_latent']
    adata_corrd.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_corrd.obs['batch'] = np.array(adata.obs['batch'])
    adata_corrd.obs['condition'] = np.array(adata.obs['condition'])

    adata_corrd.write_h5ad(data_dir_ + "simulate_cell_blast.h5ad")
