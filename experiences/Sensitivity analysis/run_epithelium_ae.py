import scanpy as sc
from run_scDisco_ae import run_scDisco
import torch
from pylab import *
from util.data_utils import preprocess
from util.utils import set_seed
from time import time
from memory_profiler import profile
#%% 1. read data
device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
data_dir = '../../datasets/epithelium/'
model_dir = data_dir + 'model_dic/' + 'epithelium_model_ae.pth'
aspect_keys = ['disease']
#%% 2. running human epithelium by scDisco-AE
@profile
def my_func(adata):
    set_seed(2023)
    start = time()
    adata, record = run_scDisco(adata, batch_key="batch", aspect_key=['disease'],
                                  cl_type='celltype', n_epochs=50,
                                  model_dir=model_dir
                                  )
    end = time()
    print('elapsed{:.2f} seconds'.format(end - start))
    adata.obsm['X_latent'] = adata.obsm['bio_feat']
    return adata

if __name__ == '__main__':
    set_seed(2023)
    data = sc.read_h5ad(data_dir + 'human_epithelium.h5ad')
    data.obs['condition'] = data.obs['disease']
    adata = preprocess(data, size_factor=True)
    adata.raw = adata
    adata = my_func(adata)
    adata.write_h5ad(data_dir + "epithelium_scDisco_ae.h5ad")
