import torch
from pylab import *
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from util.utils import set_seed, calculate_metric
from collections import Counter

data_dir = 'C://0scDisco/datasets/simulate/'
device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")

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

methods = ["seurat", "harmony", 'scanorama', 'desc', 'scVI', 'cell_blast', 'scidrl',
           'scDisinfact', 'scINSIGHT', "scDisco"]

Raw = sc.read_h5ad(data_dir + 'sim_count.h5ad')
celltypes = np.unique(Raw.obs['celltype'])
cell_type = Raw.obs['celltype'].values
Label_enc = LabelEncoder()
cell_type = Label_enc.fit_transform(cell_type)
Raw.obs['clustering'] = np.array(cell_type)

def count_same_labels(cell_type, y_pred_k):
    num_same_labels = sum(1 for x, y in zip(cell_type, y_pred_k) if x == y if np.isclose(x, y))
    return num_same_labels


for celltype in celltypes:
    raw = Raw[Raw.obs['celltype'] == celltype, :]
    cell_type = np.array(raw.obs['clustering'])
    label_encoder = LabelEncoder()
    cell_type_encoded = label_encoder.fit_transform(cell_type)
    all_same_labels_count = []
    all_ari = []
    for method in methods:
        data = sc.read_h5ad(data_dir + 'simulate_' + method + '.h5ad')
        adata = data[data.obs['celltype'] == celltype, :]


        y_pred_k = np.array(adata.obs['clustering'])
        label_counts = Counter(y_pred_k)
        label_mapping = {label: idx for idx, (label, _) in
                         enumerate(sorted(label_counts.items(), key=lambda x: x[1], reverse=True))}
        y_pred = np.array([label_mapping[label] for label in y_pred_k])


        nmi_k, ari_k = calculate_metric(cell_type_encoded, y_pred)
        all_ari.append(ari_k)
        same_labels_count = count_same_labels(cell_type_encoded, y_pred)/len(y_pred)
        all_same_labels_count.append(same_labels_count)

    A = pd.DataFrame({"same_labels_count": np.array(all_same_labels_count)})
    A.to_csv(data_dir + 'same_labels_count_' + celltype + ".csv")
    B = pd.DataFrame({"ARI": np.array(all_ari)})
    B.to_csv(data_dir + 'ari_' + celltype + ".csv")
