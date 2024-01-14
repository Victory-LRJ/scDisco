import torch
from pylab import *
import scanpy as sc
import numpy as np
from sklearn.preprocessing import LabelEncoder
from util.utils import set_seed, kmeans


data_dir = '../../datasets/simulate/'
device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
aspect_keys = ['disease']

set_seed(2023)
Raw = sc.read_h5ad(data_dir + 'sim_count.h5ad')
Raw.obsm['X_latent'] = Raw.X
Harmony = sc.read_h5ad(data_dir + 'simulate_harmony.h5ad')
DESC = sc.read_h5ad(data_dir + 'simulate_desc.h5ad')
scVI = sc.read_h5ad(data_dir + 'simulate_scVI.h5ad')
Seurat = sc.read_h5ad(data_dir + 'simulate_seurat.h5ad')
Scanorama = sc.read_h5ad(data_dir + 'simulate_scanorama.h5ad')
Cell_BLAST = sc.read_h5ad(data_dir + 'simulate_cell_blast.h5ad')
SCIDRL = sc.read_h5ad(data_dir + 'simulate_scidrl.h5ad')
scDisInFact = sc.read_h5ad(data_dir + 'simulate_scDisinfact.h5ad')
scINSIGHT = sc.read_h5ad(data_dir + 'simulate_scINSIGHT.h5ad')
scDisco = sc.read_h5ad(data_dir + 'simulate_scDisco.h5ad')
celltypes = np.unique(Raw.obs['celltype'])

methods = ["harmony", 'scanorama', 'scVI', 'cell_blast', 'scidrl', 'scDisinfact',  "scDisco"]

for method in methods:
    adata = sc.read_h5ad(data_dir + 'simulate_' + method + ".h5ad")
    cell_type = Raw.obs['celltype'].values
    Label_enc = LabelEncoder()
    cell_type = Label_enc.fit_transform(cell_type)
    n_clusters = len(np.unique(np.array(adata.obs['celltype'])))
    adata = kmeans(adata, n_clusters, use_rep='X_latent')
    y_pred_k = np.array(adata.obs['kmeans'])
    adata.obs['clustering'] = y_pred_k
    adata.write_h5ad('C://0scDisco/datasets/simulate/' + 'simulate_' + method + '.h5ad')



