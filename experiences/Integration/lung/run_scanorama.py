import numpy as np
import pandas as pd
import scanpy as sc
import scanorama
from time import time
import anndata as ad
from util.utils import set_seed
from memory_profiler import profile

# # 1. Running 10 sampled data ---------------------------------------------------------
data_dir = '../../datasets/human_lung/'
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

    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_genes=600)
    cell_norms = np.linalg.norm(adata.X, axis=1, ord=2)
    adata.X /= cell_norms[:, None]

    adata_corr = []
    genes_ = []
    all_batch = list(set(adata.obs['batch']))
    for b in all_batch:
        adata_corr.append(adata.X[adata.obs['batch'] == b, :])
        genes_.append(adata.var_names)

    integrated, corrected, genes = scanorama.correct(adata_corr, genes_, return_dimred=True)

    scanorama_res = np.concatenate(integrated)
    inted = pd.DataFrame(scanorama_res)
    adata_inted = ad.AnnData(inted, obs=adata.obs, dtype='float64')
    adata_inted.obsm['X_latent'] = adata_inted.X
    adata_inted.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_inted.obs['batch'] = np.array(adata.obs['batch'])
    adata_inted.obs['condition'] = np.array(adata.obs['condition'])

    adata_inted.write_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'lung_scanorama.h5ad')

# # 2. Running Complete data by scanorama -------------------------------------------------------
@profile
def my_func(adata):
    set_seed(2023)
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.filter_cells(adata, min_genes=600)
    cell_norms = np.linalg.norm(adata.X, axis=1, ord=2)
    adata.X /= cell_norms[:, None]
    
    adata_corr = []
    genes_ = []
    all_batch = list(set(adata.obs['batch']))
    for b in all_batch:
        adata_corr.append(adata.X[adata.obs['batch'] == b, :])
        genes_.append(adata.var_names)

    start = time()
    integrated, corrected, genes = scanorama.correct(adata_corr, genes_, return_dimred=True)
    end = time()
    print('elapsed{:.2f} seconds'.format(end - start))
    return integrated

if __name__ == '__main__':
    set_seed(2023)
    adata = sc.read_h5ad(data_dir + 'adata_lung.h5ad')

    integrated = my_func(adata)

    scanorama_res = np.concatenate(integrated)
    inted = pd.DataFrame(scanorama_res)
    adata_inted = ad.AnnData(inted, obs=adata.obs, dtype='float64')
    adata_inted.obsm['X_latent'] = adata_inted.X
    adata_inted.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_inted.obs['batch'] = np.array(adata.obs['batch'])
    adata_inted.obs['condition'] = np.array(adata.obs['age'])

    adata_inted.write_h5ad(data_dir + "lung_scanorama.h5ad")
