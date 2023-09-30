import harmonypy as hm
import pandas as pd
import scanpy as sc
import numpy as np
from time import time
import anndata as ad
from util.utils import set_seed
from memory_profiler import profile

data_dir = '../../datasets/human_ductal/'

# # 1. Running Complete data by harmony -------------------------------------------------------
@profile
def my_func(adata):
    set_seed(2023)
    start = time()
    ho = hm.run_harmony(adata.X, adata.obs, ['batch'])
    end = time()
    print('elapsed{:.2f} seconds'.format(end - start))
    return ho

if __name__ == '__main__':
    set_seed(2023)
    adata = sc.read_h5ad(data_dir + 'human_ductal.h5ad')
    adata.X = adata.X.todense()

    ho = my_func(adata)

    corrd = pd.DataFrame(ho.Z_corr.T)
    adata_corrd = ad.AnnData(corrd, obs=adata.obs, var=adata.var, dtype='float64')
    adata_corrd.obsm['X_latent'] = adata_corrd.X
    adata_corrd.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_corrd.obs['batch'] = np.array(adata.obs['batch'])
    adata_corrd.obs['condition'] = np.array(adata.obs['disease'])

    adata_corrd.write_h5ad(data_dir + "ductal_harmony.h5ad")
