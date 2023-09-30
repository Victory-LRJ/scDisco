import scanpy as sc
from model.run_scDisco import run_scDisco
import torch
from pylab import *
from util.data_utils import preprocess
from util.utils import set_seed
from time import time
from memory_profiler import profile
import anndata as ad


data_dir = '../../../datasets/lung/'
model_dir = data_dir + 'lung_model.pth'
device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
aspect_keys = ['disease']

# # 1. Running 10 sampled data
seed = [0, 1, 4, 6, 7, 8, 10, 11, 100, 111]
for i in seed:
    set_seed(2023)
    d = np.load(data_dir + 'sampling/' + str(i) + '_' + 'lung_raw.npz', allow_pickle=True)

    data = ad.AnnData(d['X_latent'])
    data.obs['celltype'] = np.array(d['celltype'])
    data.obs['batch'] = np.array(d['batch'])
    data.obs['condition'] = np.array(d['condition'])

    data.obs['celltype'] = data.obs['celltype'].astype('category')
    data.obs['batch'] = data.obs['batch'].astype('category')
    data.obs['condition'] = data.obs['condition'].astype('category')

    adata = preprocess(data, size_factor=True)
    adata.raw = adata
    adata, record = run_scDisco(adata, batch_key="batch", aspect_key=['condition'],
                                        cl_type='celltype', n_epochs=50,
                                        model_dir=model_dir
                                        )
    adata.obsm['X_latent'] = adata.obsm['bio_feat']

    adata.write_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'lung_scDisco.h5ad')

# # 2. Running Complete data
@profile
def my_func(adata):
    set_seed(2023)
    start = time()
    adata, record = run_scDisco(adata, batch_key="batch", aspect_key=['condition'],
                                  cl_type='celltype', n_epochs=50,
                                  model_dir=model_dir
                                  )
    end = time()
    print('elapsed{:.2f} seconds'.format(end - start))
    adata.obsm['X_latent'] = adata.obsm['bio_feat']   

    return adata

if __name__ == '__main__':
    set_seed(2023)
    data = sc.read_h5ad(data_dir + 'human_lung.h5ad')
    data.obs['condition'] = data.obs['age']
    adata = preprocess(data, size_factor=True)
    adata.raw = adata

    adata = my_func(adata)

    adata.write_h5ad(data_dir + "lung_scDisco.h5ad")
