import os              
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi
from util.utils import set_seed
from time import time
from memory_profiler import profile

os.environ['PYTHONHASHSEED'] = '0'
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

# ## conda activate scvi-env

# # 1. Running 10 sampled data ---------------------------------------------------------
data_dir = '../../datasets/simulate/'
seed = [0, 1, 3, 5, 10, 11, 13, 15, 100, 111]
for i in seed:
    set_seed(2023)
    adata = sc.read_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'simulate_raw.h5ad')

    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep full dimension safe

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=2000,
        layer="counts",
        batch_key="batch",
        subset=True,
    )
    net_adata = adata[:, adata.var['highly_variable']].copy()
    scvi.model.SCVI.setup_anndata(net_adata, layer="counts", batch_key="batch")
    vae = scvi.model.SCVI(net_adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()

    umap_embeddings = adata.obsm["X_scVI"]
    inted = pd.DataFrame(umap_embeddings)
    adata_inted = ad.AnnData(inted, obs=adata.obs, dtype='float64')
    adata_inted.obsm['X_latent'] = adata.obsm['X_scVI']
    adata_inted.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_inted.obs['batch'] = np.array(adata.obs['batch'])
    adata_inted.obs['condition'] = np.array(adata.obs['condition'])

    adata_inted.write_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'simulate_scVI.h5ad')

# # 2. Running Complete data by scVI -------------------------------------------------------
@profile
def my_func(adata):
    start = time()
    scvi.model.SCVI.setup_anndata(net_adata, layer="counts", batch_key="batch")
    vae = scvi.model.SCVI(net_adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    end = time()
    print('elapsed{:.2f} seconds'.format(end - start))

    adata.obsm["X_scVI"] = vae.get_latent_representation()

    return adata

if __name__ == '__main__':
    set_seed(2023)
    adata = sc.read_h5ad(data_dir + 'sim_count.h5ad')
    aspect_key = ['condition']

    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep full dimension safe

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=2000,
        layer="counts",
        batch_key="batch",
        subset=True,
    )
    net_adata = adata[:, adata.var['highly_variable']].copy()

    adata = my_func(net_adata)

    umap_embeddings = adata.obsm["X_scVI"]
    inted = pd.DataFrame(umap_embeddings)
    adata_inted = ad.AnnData(inted, obs=adata.obs, dtype='float64')
    adata_inted.obsm['X_latent'] = adata.obsm['X_scVI']
    adata_inted.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_inted.obs['batch'] = np.array(adata.obs['batch'])
    adata_inted.obs['condition'] = np.array(adata.obs['condition'])

    adata_inted.write_h5ad(data_dir + "simulate_scVI.h5ad")
