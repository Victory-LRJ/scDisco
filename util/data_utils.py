import scanpy as sc
import numpy as np

def preprocess(adata, filter_min_counts=True, size_factor=True, logtrans_input=True, normalize_input=True, high_var_select=True):

    if filter_min_counts:
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=1)

    if size_factor:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        adata.obs['size_factor'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    else:
        adata.obs['size_factor'] = 1.0

    if logtrans_input:
        sc.pp.log1p(adata)
        if len(adata.var) >= 3000 and high_var_select:
            sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat')
            adata = adata[:, adata.var['highly_variable']]

    if normalize_input:
        sc.pp.scale(adata)

    if size_factor or normalize_input or logtrans_input:
        adata.raw = adata.copy()
    else:
        adata.raw = adata

    return adata

def getdims(x):
    n_sample = x.shape[0]
    if n_sample > 30000:
        dims = 30
    elif n_sample > 20000:
        dims = 28
    elif n_sample > 5000:
        dims = 22
    else:
        dims = 20
    return dims

