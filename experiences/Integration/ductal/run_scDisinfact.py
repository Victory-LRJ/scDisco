import scanpy as sc
import torch
import numpy as np
import scDisInFact
from scDisInFact import scdisinfact
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import anndata as ad
from util.utils import set_seed
from time import time
from memory_profiler import profile

device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")

data_dir = '../../datasets/human_ductal/'
# # 1. Running Complete data
@profile
def my_func(data):
    set_seed(2023)
    counts = data.X
    meta_cells = data.obs[['batch', 'disease']]
    Label_enc = LabelEncoder()
    condition = Label_enc.fit_transform(data.obs['disease'])
    meta_cells['disease'] = np.array(condition)
    Label_enc = LabelEncoder()
    batch = Label_enc.fit_transform(data.obs['batch'])
    meta_cells['batch'] = np.array(batch)

    start = time()
    data_dict = scDisInFact.create_scdisinfact_dataset(counts, meta_cells, condition_key=['disease'],
                                                       batch_key="batch")
    Ks = [8, 4]
    # training device
    model = scdisinfact(data_dict=data_dict, Ks=Ks, device=device)
    losses = model.train_model(nepochs=100)
    _ = model.eval()

    z_cs = []
    z_ds = []
    zs = []
    for dataset in data_dict["datasets"]:
        with torch.no_grad():
            # pass through the encoders
            dict_inf = model.inference(counts=dataset.counts_norm.to(model.device),
                                       batch_ids=dataset.batch_id[:, None].to(model.device), print_stat=True)
            # pass through the decoder
            dict_gen = model.generative(z_c=dict_inf["mu_c"], z_d=dict_inf["mu_d"],
                                        batch_ids=dataset.batch_id[:, None].to(model.device))
            z_c = dict_inf["mu_c"]
            z_d = dict_inf["mu_d"]
            mu = dict_gen["mu"]
            z_ds.append([x.cpu().detach().numpy() for x in z_d])
            z_cs.append(z_c.cpu().detach().numpy())
    # shared-bio factor, concatenate across all training matrices
    z_cs = np.concatenate(z_cs, axis=0)
    z_ds = np.concatenate([x[0] for x in z_ds], axis=0)
    end = time()
    print('elapsed{:.2f} seconds'.format(end - start))
    return z_cs, z_ds

if __name__ == '__main__':
    set_seed(2023)
    data = sc.read_h5ad(data_dir + 'human_ductal.h5ad')
    z_cs, z_ds = my_func(data)

    inted = pd.DataFrame(z_cs)
    adata_inted = ad.AnnData(inted, obs=data.obs, dtype='float64')
    adata_inted.obsm['X_latent'] = adata_inted.X
    adata_inted.obsm['X_cond'] = z_ds
    adata_inted.obs['batch'] = np.array(data.obs['batch'])
    adata_inted.obs['celltype'] = np.array(data.obs['celltype'])
    adata_inted.obs['condition'] = np.array(data.obs['disease'])

    adata_inted.write_h5ad(data_dir + 'ductal_scDisinfact.h5ad')



