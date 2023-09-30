import scanpy as sc
from model.scDisco import scDisco
from condition_genes import *
from sklearn.preprocessing import LabelEncoder
from util.utils import set_seed
import matplotlib
from util.data_utils import getdims

plt.rc('font', family='Arial')
matplotlib.use('Agg')
device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
# %% 1. loading integrated data
data_dir = '../../datasets/ductal/'
model_dir = data_dir + 'model_dic/human_ductal_model.pth'
data = sc.read_h5ad(data_dir + 'human_ductal.h5ad')
adata = sc.read_h5ad(data_dir + 'ductal_scDisco.h5ad')
used_condition_types = np.unique(adata.obs["condition"])

# %% 2. plot condition density by embedding
sc.pp.neighbors(adata, use_rep='cond_feat_condition')
sc.tl.umap(adata)
sc.tl.embedding_density(adata, groupby='condition', key_added='condition_density')
sc.pl.embedding_density(adata, key='condition_density')
plt.savefig(data_dir + 'figures/condition/figure12A.png', dpi=600, format='png', bbox_inches='tight')

# %% 3. loading model dict
set_seed(2023)
aspect_key = ['condition']
aspects = {aspect: adata.obs[aspect].values.to_numpy() for aspect in aspect_key}
n_conditions = {aspect: len(np.unique(conditions)) for (aspect, conditions) in aspects.items()}

model = scDisco(input_dim=adata.X.shape[1],
                       bio_z_dim=getdims(adata.X),
                       cond_dim=8,
                       EncLayer=[512, 256],
                       DecLayer=[256, 512],
                       n_batches=len(np.unique(adata.obs['batch'])),
                       n_conditions=n_conditions,
                       aspect_key=aspect_key,
                       ).to(device)

condition_types = {}
condition_labels = {}
for (aspect, conditions) in aspects.items():
    Label_enc = LabelEncoder()
    condition_types[aspect] = Label_enc.fit_transform(conditions)
    condition_labels[aspect] = torch.tensor(condition_types[aspect])

# %% 4. important genes select
gene_grad_colums, gene_grad_index, gene_grad_colums_0, gene_grad_colums_1 = gradients(adata,
                                                                                      model=model,
                                                                                      model_dir=model_dir,
                                                                                      aspect='condition',
                                                                                      c_labels=condition_labels)

index_double_f = mutigenes_sort(used_condition_types=used_condition_types,
                                gene_grad_colums_1=gene_grad_colums_1)
print(index_double_f)
marker_genes_dict = {
    'N': index_double_f[0],
    'T': index_double_f[1],
}
sc.pl.matrixplot(adata, marker_genes_dict, 'condition', dendrogram=True,
                 cmap='Blues', standard_scale='var', colorbar_title='column scaled\nexpression')
plt.savefig(data_dir + 'figures/condition/figure12B.svg', dpi=300, format='svg', bbox_inches='tight')



