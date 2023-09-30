import harmonypy as hm
import pandas as pd
import scanpy as sc
import numpy as np
from time import time
import anndata as ad
from util.utils import set_seed
from memory_profiler import profile

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

    ho = hm.run_harmony(adata.X, adata.obs, ['batch'])

    corrd = pd.DataFrame(ho.Z_corr.T)
    adata_corrd = ad.AnnData(corrd, obs=adata.obs, var=adata.var, dtype='float64')
    adata_corrd.obsm['X_latent'] = adata_corrd.X
    adata_corrd.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_corrd.obs['batch'] = np.array(adata.obs['batch'])
    adata_corrd.obs['condition'] = np.array(adata.obs['condition'])

    adata_corrd.write_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'lung_harmony.h5ad')
#
# # # 2. Running Complete data by harmony -------------------------------------------------------
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
    adata = sc.read_h5ad(data_dir + 'adata_lung.h5ad')

    ho = my_func(adata)

    corrd = pd.DataFrame(ho.Z_corr.T)
    adata_corrd = ad.AnnData(corrd, obs=adata.obs, var=adata.var, dtype='float64')
    adata_corrd.obsm['X_latent'] = adata_corrd.X
    adata_corrd.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_corrd.obs['batch'] = np.array(adata.obs['batch'])
    adata_corrd.obs['condition'] = np.array(adata.obs['age'])

    adata_corrd.write_h5ad(data_dir + "lung_harmony.h5ad")
