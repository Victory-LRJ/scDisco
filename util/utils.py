import typing
import numpy as np
import random
import scanpy as sc
import torch
from sklearn.cluster import KMeans
from sklearn import metrics, neighbors
from sklearn.metrics import silhouette_score

def set_seed(seed):
    np.random.seed(seed)
    random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

def calculate_metric(pred, label):
    nmi = np.round(metrics.normalized_mutual_info_score(label, pred), 5)
    ari = np.round(metrics.adjusted_rand_score(label, pred), 5)
    return nmi, ari

def kmeans(adata, n_clusters, use_rep=None):
    k_means = KMeans(n_clusters, n_init=20)
    y_pred = k_means.fit_predict(adata.obsm[use_rep])
    adata.obs['kmeans'] = y_pred
    adata.obs['kmeans'] = adata.obs['kmeans'].astype(str).astype('category')
    return adata

def louvain(adata, resolution=None, use_rep=None):
    sc.pp.neighbors(adata, use_rep=use_rep)
    sc.tl.louvain(adata, resolution=resolution)
    return adata

def silhouette_coeff_ASW(adata, percent_extract=None):
    np.random.seed(0)
    asw_f1score = []
    asw_bn = []
    asw_ctn = []
    asw_ = []
    iters = []
    for i in range(20):
        iters.append('iteration_' + str(i + 1))
        rand_cidx = np.random.choice(adata.obs_names, size=int(len(adata.obs_names) * percent_extract), replace=False)
        adata_ext = adata[rand_cidx, :]
        asw_batch = silhouette_score(adata_ext.obsm['X_latent'], adata_ext.obs['batch'])
        asw_celltype = silhouette_score(adata_ext.obsm['X_latent'], adata_ext.obs['celltype'])
        min_val = -1
        max_val = 1
        asw_batch_norm = abs(asw_batch)
        asw_celltype_norm = (asw_celltype - min_val) / (max_val - min_val)

        f1scoreASW = (2 * (1 - asw_batch_norm) * (asw_celltype_norm)) / (1 - asw_batch_norm + asw_celltype_norm)

        asw = 0.8 * (1 - asw_batch_norm) + 0.2 * asw_celltype_norm

        asw_f1score.append(f1scoreASW)
        asw_bn.append(1 - asw_batch_norm)
        asw_ctn.append(asw_celltype_norm)
        asw_.append(asw)
    BASW = np.array(asw_bn).mean()
    CASW = np.array(asw_ctn).mean()
    F1_Score = np.array(asw_f1score).mean()
    ASW = np.array(asw_).mean()
    return F1_Score, BASW, CASW, ASW


def silhouette_coeff_ASW_single(adata):
    np.random.seed(2023)
    asw_batch = silhouette_score(adata.obsm['X_latent'], adata.obs['batch'])
    asw_celltype = silhouette_score(adata.obsm['X_latent'], adata.obs['celltype'])
    min_val = -1
    max_val = 1
    asw_batch_norm = abs(asw_batch)
    asw_celltype_norm = (asw_celltype - min_val) / (max_val - min_val)

    fscoreASW = (2 * (1 - asw_batch_norm) * (asw_celltype_norm)) / (1 - asw_batch_norm + asw_celltype_norm)
    asw = 0.8 * (1 - asw_batch_norm) + 0.2 * asw_celltype_norm

    BASW = 1 - asw_batch_norm
    CASW = asw_celltype_norm
    F_Score = fscoreASW
    ASW = asw
    return F_Score, BASW, CASW, ASW
