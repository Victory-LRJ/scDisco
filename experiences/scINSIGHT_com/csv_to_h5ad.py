import pandas as pd
import anndata as ad
import numpy as np

data_dir = '../../datasets/'
#%% 1.  simulated data embedding by scINSIGHT from csv to h5ad
data_1 = pd.read_csv('simulate/simulate_scINSIGHT.csv')
meta_1 = pd.read_csv("simulate/simulate_scINSIGHT_clust.csv")

corrd = pd.DataFrame(data_1.iloc[:, 1:9], dtype=float)
adata_corrd = ad.AnnData(corrd, dtype='float64')
adata_corrd.obs_names = [str(i) for i in data_1.iloc[:, 0]]
adata_corrd.obsm['X_latent'] = adata_corrd.X
adata_corrd.obs['celltype'] = np.array(data_1['celltype'])
adata_corrd.obs['batch'] = np.array(data_1['batch'])
adata_corrd.obs['condition'] = np.array(data_1['condition'])
adata_corrd.obs['clustering'] = np.array(meta_1.iloc[:, 1])
adata_corrd.write_h5ad(data_dir + 'simulate_scINSIGHT.h5ad')

#%% 2.  human pancreas2 embedding by scINSIGHT from csv to h5ad
data_2 = pd.read_csv('pancreas2/pancreas2_scINSIGHT.csv')
meta_2 = pd.read_csv("pancreas2/pancreas2_scINSIGHT_clust.csv")

corrd = pd.DataFrame(data_2.iloc[:, 1:13], dtype=float)
adata_corrd = ad.AnnData(corrd, dtype='float64')
adata_corrd.obs_names = [str(i) for i in data_2.iloc[:, 0]]
adata_corrd.obsm['X_latent'] = adata_corrd.X
adata_corrd.obs['celltype'] = np.array(data_2['celltype'])
adata_corrd.obs['batch'] = np.array(data_2['batch'])
adata_corrd.obs['condition'] = np.array(data_2['disease'])
adata_corrd.obs['clustering'] = np.array(meta_2.iloc[:, 1])
adata_corrd.write_h5ad(data_dir + 'pancreas2_scINSIGHT.h5ad')

#%% 3.  human lung embedding by scINSIGHT from csv to h5ad
data_3 = pd.read_csv('lung/lung_scINSIGHT.csv')
meta_3 = pd.read_csv("lung/lung_scINSIGHT_clust.csv")

corrd = pd.DataFrame(data_3.iloc[:, 1:11], dtype=float)
adata_corrd = ad.AnnData(corrd, dtype='float64')
adata_corrd.obs_names = [str(i) for i in data_3.iloc[:, 0]]
adata_corrd.obsm['X_latent'] = adata_corrd.X
adata_corrd.obs['celltype'] = np.array(data_3['celltype'])
adata_corrd.obs['batch'] = np.array(data_3['batch'])
adata_corrd.obs['condition'] = np.array(data_3['age'])
adata_corrd.obs['clustering'] = np.array(meta_3.iloc[:, 1])
adata_corrd.write_h5ad(data_dir + 'lung_scINSIGHT.h5ad')

#%% 4.  human pancreas embedding by scINSIGHT from csv to h5ad
data_4 = pd.read_csv('pancreas/pancreas_scINSIGHT.csv')
meta_4 = pd.read_csv("pancreas/pancreas_scINSIGHT_clust.csv")

corrd = pd.DataFrame(data_4.iloc[:, 1:9], dtype=float)
adata_corrd = ad.AnnData(corrd, dtype='float64')
adata_corrd.obs_names = [str(i) for i in data_4.iloc[:, 0]]
adata_corrd.obsm['X_latent'] = adata_corrd.X
adata_corrd.obs['celltype'] = np.array(data_4['celltype'])
adata_corrd.obs['batch'] = np.array(data_4['batch'])
adata_corrd.obs['condition'] = np.array(data_4['type2_mellitus'])
adata_corrd.obs['clustering'] = np.array(meta_4.iloc[:, 1])
adata_corrd.write_h5ad(data_dir + 'pancreas_scINSIGHT.h5ad')

#%% 5.  mouse mucosa embedding by scINSIGHT from csv to h5ad
data_5 = pd.read_csv('mouse_mucosa/mouse_mucosa_scINSIGHT.csv')
meta_5 = pd.read_csv("mouse_mucosa/mouse_mucosa_scINSIGHT_clust.csv")

corrd = pd.DataFrame(data_5.iloc[:, 1:9], dtype=float)
adata_corrd = ad.AnnData(corrd, dtype='float64')
adata_corrd.obs_names = [str(i) for i in data_5.iloc[:, 0]]
adata_corrd.obsm['X_latent'] = adata_corrd.X
adata_corrd.obs['celltype'] = np.array(data_5['celltype'])
adata_corrd.obs['batch'] = np.array(data_5['batch'])
adata_corrd.obs['condition'] = np.array(data_5['stimulus'])
adata_corrd.obs['clustering'] = np.array(meta_5.iloc[:, 1])
adata_corrd.write_h5ad(data_dir + 'mouse_mucosa_scINSIGHT.h5ad')