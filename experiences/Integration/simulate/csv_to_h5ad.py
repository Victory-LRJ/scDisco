import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc

# ### #### read data
data_dir = '../../datasets/simulate/'
sim = sc.read_h5ad(data_dir + 'sim_count.h5ad')
seed = [0, 1, 4, 6, 7, 8, 10, 11, 100, 111]
for i in seed:
    data = pd.read_csv(data_dir + 'sampling/' + str(i) + '_' + 'simulate_seurat.csv')
    meta = pd.read_csv(data_dir + 'sampling/' + str(i) + '_' + "simulate_seurat_clust.csv")

    corrd = pd.DataFrame(data.iloc[:, 1:50], dtype=float)
    adata_corrd = ad.AnnData(corrd, dtype='float64')
    adata_corrd.obs_names = [str(i) for i in data.iloc[:, 0]]
    adata_corrd.obsm['X_latent'] = adata_corrd.X
    adata_corrd.obs['celltype'] = np.array(data['simulate.combined.meta.data...celltype...'])
    adata_corrd.obs['batch'] = np.array(data['simulate.combined.meta.data...batch...'])
    adata_corrd.obs['clustering'] = np.array(meta.iloc[:, 1])
    adata_corrd.write_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'simulate_seurat.h5ad')

data = pd.read_csv(data_dir + 'simulate_seurat.csv')
meta = pd.read_csv(data_dir + "simulate_seurat_clust.csv")

corrd = pd.DataFrame(data.iloc[:, 1:50], dtype=float)
adata_corrd = ad.AnnData(corrd, dtype='float64')
adata_corrd.obs_names = [str(i) for i in data.iloc[:, 0]]
adata_corrd.obsm['X_latent'] = adata_corrd.X
adata_corrd.obs['celltype'] = np.array(data['simulate.combined.meta.data...celltype...'])
adata_corrd.obs['batch'] = np.array(data['simulate.combined.meta.data...batch...'])
adata_corrd.obs['clustering'] = np.array(meta.iloc[:, 1])
adata_corrd.obs['condition'] = np.array(sim.obs['condition'])
adata_corrd.write_h5ad(data_dir + 'simulate_seurat.h5ad')

data = pd.read_csv(data_dir + 'simulate_scINSIGHT.csv')
meta = pd.read_csv(data_dir + "simulate_scINSIGHT_clust.csv")

corrd = pd.DataFrame(data.iloc[:, 1:9], dtype=float)
adata_corrd = ad.AnnData(corrd, dtype='float64')
adata_corrd.obs_names = [str(i) for i in data.iloc[:, 0]]
adata_corrd.obsm['X_latent'] = adata_corrd.X
adata_corrd.obs['celltype'] = np.array(data['celltype'])
adata_corrd.obs['batch'] = np.array(data['batch'])
adata_corrd.obs['condition'] = np.array(data['condition'])
adata_corrd.obs['clustering'] = np.array(meta.iloc[:, 1])
adata_corrd.write_h5ad(data_dir + 'simulate_scINSIGHT.h5ad')