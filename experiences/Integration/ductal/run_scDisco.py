import scanpy as sc
from model.run_scDisco import run_scDisco
import torch
from pylab import *
from util.data_utils import preprocess
from util.utils import set_seed
from time import time
from memory_profiler import profile

data_dir = '../../../datasets/ductal/'
model_dir = data_dir + 'model_dic/human_ductal_model.pth'
device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
aspect_keys = ['disease']

# 1. Running Complete data
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
    data = sc.read_h5ad(data_dir + 'human_ductal.h5ad')
    data.obs['condition'] = data.obs['disease']
    adata = preprocess(data, size_factor=True)
    adata.raw = adata

    adata = my_func(adata)

    adata.write_h5ad(data_dir + "ductal_scDisco.h5ad")
