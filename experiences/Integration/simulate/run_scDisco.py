import scanpy as sc
from model.run_scDisco import run_scDisco
import torch
from util.utils import set_seed
from time import time
import numpy as np
import anndata as ad


device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")

data_dir = '../../datasets/simulate/'
model_dir = data_dir + 'simulate_model.pth'
# 1. Running 10 sampled data
seed = [0, 1, 4, 6, 7, 8, 10, 11, 100, 111]
for i in seed:
    set_seed(2023)
    d = np.load(data_dir + 'sampling/' + str(i) + '_' + 'simulate_raw.npz', allow_pickle=True)

    data = ad.AnnData(d['X_latent'])
    data.obs['celltype'] = np.array(d['celltype'])
    data.obs['batch'] = np.array(d['batch'])
    data.obs['condition'] = np.array(d['condition'])

    data.obs['celltype'] = data.obs['celltype'].astype('category')
    data.obs['batch'] = data.obs['batch'].astype('category')
    data.obs['condition'] = data.obs['condition'].astype('category')

    data.raw = data

    adata, record = run_scDisco(data, batch_key="batch", aspect_key=['condition'],
                                        cl_type='celltype', n_epochs=300, model_dir=model_dir
                                        )
    adata.obsm['X_latent'] = adata.obsm['bio_feat']
    adata.write_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'simulate_scDisco.h5ad')

# # 2. Running Complete data
data_dir_ = '../../datasets/simulate/'
data = sc.read_h5ad(data_dir_ + 'sim_count.h5ad')
data.raw = data
aspect_key = ['condition']
set_seed(2023)
start = time()
adata, record = run_scDisco(data, batch_key="batch", aspect_key=['condition'],
                              cl_type='celltype', n_epochs=300, model_dir=model_dir
                              )
end = time()
print('elapsed{:.2f} seconds'.format(end - start))


adata.obsm['X_latent'] = adata.obsm['bio_feat']
adata.obsm['X_cond'] = adata.obsm['cond_feat_condition']
adata.write_h5ad(data_dir_ + 'simulate_scDisco.h5ad')
