# scDisco v1.0

The code for "scDisco: Integration of scRNA-seq Data by Disentangled Representation Learning with Condition Domain Adaptation".
The sources of all preprocessed data used in this work are available at https://drive.google.com/drive/folders/1OCN6UmUsM98CpsecpbmQZsXmS0HKcB4k?usp=drive_link.

The figure of scDisco is shown below.

[![](file:///model.png)](https://github.com/Victory-LRJ/scDisco/blob/main/model.png?raw=true)

#### demo

scDisco is implemented in Python 3.8, using Scanpy version 1.9.3 for data preprocessing. The code execution follows the following steps:

```python
import scanpy as sc
from model.run_scDisco import run_scDisco
import torch
from pylab import *
from util.data_utils import preprocess
from util.utils import set_seed
import anndata as ad

data_dir = '../../../datasets/lung/'
model_dir = data_dir + 'lung_model.pth'
set_seed(2023)
data = sc.read_h5ad(data_dir + 'human_lung.h5ad')
data.obs['condition'] = data.obs['age']
adata = preprocess(data, size_factor=True)
adata, record = run_scDisco(adata, batch_key="batch", aspect_key=['condition'], 
                                   cl_type='celltype', n_epochs=50, 
                                   model_dir=model_dir)
adata.obsm['X_latent'] = adata.obsm['bio_feat']
adata.write_h5ad(data_dir + "lung_scDisco.h5ad")
```
