import torch
from torch.utils.data import Dataset
from sklearn.preprocessing import OneHotEncoder
import numpy as np
import pandas as pd

class Batch_Dataset(Dataset):
    def __init__(self, raw_mat, exp_mat, aspects, batches, size_factor, condition_types, batch_labels):
        super(Batch_Dataset).__init__()
        # self.raw_mat = raw_mat.toarray()
        self.raw_mat = raw_mat.astype(np.float32)
        self.exp_mat = exp_mat.astype(np.float32)

        self.aspects = {}
        for (aspect, condition) in aspects.items():
            OneHot_enc = OneHotEncoder()
            onehot_code = OneHot_enc.fit_transform(condition.reshape(-1, 1))
            self.aspects[aspect] = onehot_code.toarray()

        self.batch_labels = torch.tensor(batch_labels)


        OneHot_enc = OneHotEncoder()
        onehot_code = OneHot_enc.fit_transform(batches.reshape(-1, 1))
        batches = onehot_code.toarray()
        self.batches = batches.astype(np.float32)

        self.size_factor = size_factor.astype(np.float32)

        self.condition_types = {}
        for (aspect, data) in condition_types.items():
            self.condition_types[aspect] = torch.tensor(data)
        self.n_conditions = [len(set(condition)) for (aspect, condition) in aspects.items()]


    def __len__(self):
        return len(self.size_factor)

    def __getitem__(self, idx):
        aspects = {aspect: cond[idx, :].astype(np.float32) for (aspect, cond) in self.aspects.items()}
        condition_types = {aspect: data[idx] for (aspect, data) in self.condition_types.items()}

        return idx, self.raw_mat[idx, :], self.exp_mat[idx, :], aspects, self.batches[idx, :], self.size_factor[idx], condition_types, self.batch_labels[idx]

def data_sample(adata, celltypes, b_labels, c_labels, seed):
    x_sample = []
    y_sample = []
    b_sample = []
    c_sample = []
    for celltype in list(np.unique(celltypes)):
        adata_celltype = adata[adata.obs["celltype"] == celltype, :]
        index = adata_celltype.obs_names
        data = pd.DataFrame(adata_celltype.X, index=index)
        data = data.sample(frac=0.95, replace=False, weights=None, random_state=seed, axis=0)
        index_ = data.index
        x_sample.append(data.values)
        y_sample.extend(celltypes[index_])
        b_sample.extend(b_labels[index_])
        c_sample.extend(c_labels[index_])
    y_sample = pd.Series(y_sample, dtype='category')
    b_sample = pd.Series(b_sample, dtype='category')
    c_sample = pd.Series(c_sample, dtype='category')

    return np.concatenate(x_sample, axis=0), y_sample, b_sample,  c_sample