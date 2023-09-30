import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc

# ### #### read data
data_dir = '../../datasets/human_ductal/'
sim = sc.read_h5ad(data_dir + 'human_ductal.h5ad')


data = pd.read_csv(data_dir + 'ductal_seurat.csv')
meta = pd.read_csv(data_dir + "ductal_seurat_clust.csv")

corrd = pd.DataFrame(data.iloc[:, 1:50], dtype=float)
adata_corrd = ad.AnnData(corrd, dtype='float64')
adata_corrd.obs_names = [str(i) for i in data.iloc[:, 0]]
adata_corrd.obsm['X_latent'] = adata_corrd.X
adata_corrd.obs['celltype'] = np.array(data['simulate.combined.meta.data...celltype...'])
adata_corrd.obs['batch'] = np.array(data['simulate.combined.meta.data...batch...'])
adata_corrd.obs['clustering'] = np.array(meta.iloc[:, 1])
adata_corrd.obs['condition'] = np.array(sim.obs['disease'])
adata_corrd.write_h5ad(data_dir + 'ductal_seurat.h5ad')
