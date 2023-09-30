import pandas as pd
import torch
from sklearn.preprocessing import OneHotEncoder
from pylab import *

device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
# ## Gradients #######################################################
seed = [2, 0, 1, 3, 2013]
def gradients(adata,
              model=None, model_dir=None,
              aspect=None, c_labels=None,
              ):

    normalize_deviation = True
    used_condition_types = np.unique(adata.obs["condition"])

    batches = adata.obs['batch'].values.to_numpy()
    OneHot_enc = OneHotEncoder()
    onehot_code = OneHot_enc.fit_transform(batches.reshape(-1, 1))
    batches_onehot = onehot_code.toarray()
    batch_onehot = batches_onehot.astype(np.float32)
    batch_onehot = pd.DataFrame(batch_onehot)
    batch_onehot['condition'] = np.array(adata.obs['condition'])

    c_types = pd.DataFrame(c_labels['condition'])
    c_types['condition'] = np.array(adata.obs['condition'])

    gene_grad = []
    gene_grad_colums_0 = pd.DataFrame()
    gene_grad_colums_1 = pd.DataFrame()
    for condition_type in used_condition_types:
        print(f"Dealing with {condition_type}...")
        used_ref = adata[adata.obs["condition"] != condition_type, :]
        used_query = adata[adata.obs["condition"] == condition_type, :]
        b = batch_onehot[batch_onehot['condition'] == condition_type]
        b = b.drop(b.columns[[-1]], axis=1)

        c = c_types[c_types['condition'] == condition_type]
        c = c.drop(c.columns[[-1]], axis=1)
        c = torch.tensor(np.array(c)).squeeze()

        ref_latent = np.ones(shape=(len(used_query), 8))
        for k in range(5):
            np.random.seed(seed[k])
            index = np.random.choice(len(used_ref), size=len(used_query), replace=True)
            s = 'cond_feat' + '_' + aspect
            ref_latent_ = np.array(adata.obsm[s][index, :])
            ref_latent = ref_latent + ref_latent_
        ref_latent = ref_latent/5

        data = pd.DataFrame(adata.obsm[s])
        data['condition'] = np.array(adata.obs['condition'])
        query_latent = np.array(data.loc[data['condition'] == condition_type])
        query_latent = np.delete(query_latent, -1, axis=1)

        deviation = query_latent - ref_latent
        deviation = deviation.astype('float')
        if normalize_deviation:
            deviation /= np.linalg.norm(deviation, axis=1, keepdims=True)  # 按行标准化?还是按列 之前是1

        model.load_state_dict(torch.load(model_dir))

        _gene_grad = model.fetch_grad(torch.from_numpy(used_query.X.astype('float32')).to(device),
                                      torch.from_numpy(np.array(b).astype('float32')).to(device),
                                      torch.from_numpy(deviation).to(device),
                                      aspect='condition', c_labels=c,
                                      )
        _gene_grad = np.average(_gene_grad.cpu(), axis=0)
        gene_grad.append(_gene_grad)
        gene_grad_colums_ = pd.DataFrame({'index' + '_' + condition_type: list(adata.var.T), condition_type: np.array(_gene_grad).T}, index=list(adata.var.T))
        gene_grad_colums_['index' + '_' + condition_type] = pd.Series(gene_grad_colums_['index' + '_' + condition_type])
        gene_grad_colums_0 = pd.concat([gene_grad_colums_0, gene_grad_colums_], axis=1)

        gene_grad_colums__ = pd.DataFrame(
            {'index': list(adata.var.T), 'gene': np.array(_gene_grad).T}, index=list(adata.var.T))
        gene_grad_colums__['index'] = condition_type + '_' + gene_grad_colums__['index']
        gene_grad_colums_1 = pd.concat([gene_grad_colums_1, gene_grad_colums__], axis=0)
    gene_grad_colums = pd.DataFrame(np.array(gene_grad).T, index=list(adata.var.T), columns=list(used_condition_types))
    gene_grad_index = pd.DataFrame(np.array(gene_grad), index=list(used_condition_types), columns=list(adata.var.T))
    return gene_grad_colums, gene_grad_index, gene_grad_colums_0, gene_grad_colums_1

# ## Genes sort #######################################################
def mutigenes_sort(used_condition_types=None,
                   gene_grad_colums_1=None
                   ):
    index_double_f = []
    # l = len(used_condition_types)
    # A1, A2, A3 = [], [], []
    A1, A2 = [], []
    sorted_gene_grad = gene_grad_colums_1.sort_values(by='gene', ascending=False)
    for i in range(len(sorted_gene_grad['gene'])):
        # if used_condition_types[0] in sorted_gene_grad['index'][i] and sorted_gene_grad.index[i] not in (set(A2)|set(A3)) and len(A1) < 15:
        if used_condition_types[0] in sorted_gene_grad['index'][i] and sorted_gene_grad.index[i] not in (set(A2)) and len(A1) < 15:
            A = sorted_gene_grad.index[i]
            A1.append(A)
        # if used_condition_types[1] in sorted_gene_grad['index'][i] and sorted_gene_grad.index[i] not in (set(A1)|set(A3)) and len(A2) < 15:
        if used_condition_types[1] in sorted_gene_grad['index'][i] and sorted_gene_grad.index[i] not in (set(A1)) and len(A2) < 15:
            A = sorted_gene_grad.index[i]
            A2.append(A)
        if len(A1) == 15 and len(A2) == 15:
            index_double_f.append(A1)
            index_double_f.append(A2)
            break
    return index_double_f


