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

data_dir = '../../datasets/human_ductal/'
# # 1. Running Complete data by cell blast -------------------------------------------------------
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
    data = sc.read_h5ad(data_dir + 'human_ductal.h5ad')

    adata = my_func(data)
    corrd = pd.DataFrame(adata.obsm['X_latent'])
    adata_corrd = ad.AnnData(corrd, obs=adata.obs, dtype='float64')
    adata_corrd.obsm['X_latent'] = adata.obsm['X_latent']
    adata_corrd.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_corrd.obs['batch'] = np.array(adata.obs['batch'])
    adata_corrd.obs['condition'] = np.array(adata.obs['disease'])

    adata_corrd.write_h5ad(data_dir + "ductal_cell_blast.h5ad")
