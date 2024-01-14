import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

col1 = ["#E64B35CC", "#0072B5CC", "#00A087CC", "#3C5488CC", "#F39B7FCC", "#F7DC05FF", "#FD7446E5",
       "#8491B4CC", "#7E6148CC", "#B09C85CC", "#E18727CC", "#FFDC91E5", "#6A6599E5", "#9467BDB2",
       "#0000FFFF", "#FF0000FF", "#00FF00FF", "#000033FF", "#FF00B6FF", "#005300FF", "#FFD300FF",
       "#009FFFFF", "#9A4D42FF", "#00FFBEFF", "#783FC1FF", "#1F9698FF", "#FFACFDFF", "#B1CC71FF", "#F1085CFF",
       "#FE8F42FF", "#DD00FFFF", "#201A01FF", "#720055FF", "#766C95FF", "#02AD24FF", "#C8FF00FF", "#886C00FF",
       "#FFB79FFF", "#858567FF", "#A10300FF", "#14F9FFFF", "#00479EFF", "#DC5E93FF", "#93D4FFFF", "#004CFFFF"]
col2 = ["#FF7F00", "#7570B3", "#B2DF8A", "#E31A1C"]

def plot_umap_all(adata, name, method, ax, show):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    sc.pp.neighbors(adata, use_rep='X_latent')
    sc.tl.umap(adata)
    if name == 'celltype':
        sc.pl.umap(adata, color=['celltype'], ax=ax, show=show, legend_loc='None')
        ax.set_title(method, fontsize=24, family='Arial')
    elif name == 'batch':
        sc.pl.umap(adata, color=['batch'], ax=ax, show=show, legend_loc='None', palette=col1)
    elif name == 'condition':
        sc.pl.umap(adata, color=['condition'], ax=ax, show=show, legend_loc='None', palette=col2)

    x_axis = ax.get_xaxis()
    y_axis = ax.get_yaxis()
    x_axis.set_visible(False)
    y_axis.set_visible(False)

# %%
data_dir = '../../../datasets/ductal/'
fig = plt.figure(figsize=(16, 10))
datasets = ['by cell type', 'by batch', 'by condition']
sub_figs = fig.subfigures(3, 1)
axs = []
for i, sub_fig in enumerate(sub_figs):
    axs.append(sub_fig.subplots(1, 4))
    sub_fig.supylabel(datasets[i], x=0.07, fontsize=24)
axs = np.array(axs)

# %%
Raw = sc.read_h5ad(data_dir + 'human_ductal.h5ad')
Raw.obsm['X_latent'] = Raw.X
Raw.obs['condition'] = Raw.obs['disease']
Harmony = sc.read_h5ad(data_dir + 'ductal_harmony.h5ad')
Scanorama = sc.read_h5ad(data_dir + 'ductal_scanorama.h5ad')
DESC = sc.read_h5ad(data_dir + 'ductal_desc.h5ad')
scVI = sc.read_h5ad(data_dir + 'ductal_scVI.h5ad')
cell_BLAST = sc.read_h5ad(data_dir + 'ductal_cell_blast.h5ad')
SCIDRL = sc.read_h5ad(data_dir + 'ductal_scidrl.h5ad')
scDisInFact = sc.read_h5ad(data_dir + 'ductal_scDisinfact.h5ad')
scDisco = sc.read_h5ad(data_dir + 'ductal_scDisco.h5ad')
# %%
plot_umap_all(Scanorama, 'celltype', 'Scanorama', axs[0][0], show=False)
plot_umap_all(Scanorama, 'batch', 'Scanorama', axs[1][0], show=False)
plot_umap_all(Scanorama, 'condition', 'Scanorama', axs[2][0], show=False)
plot_umap_all(Harmony, 'celltype', 'Harmony', axs[0][1], show=False)
plot_umap_all(Harmony, 'batch', 'Harmony', axs[1][1], show=False)
plot_umap_all(Harmony, 'condition', 'Harmony', axs[2][1], show=False)
plot_umap_all(scDisInFact, 'celltype', 'scDisInFact', axs[0][2], show=False)
plot_umap_all(scDisInFact, 'batch', 'scDisInFact', axs[1][2], show=False)
plot_umap_all(scDisInFact, 'condition', 'scDisInFact', axs[2][2], show=False)
plot_umap_all(SCIDRL, 'celltype', 'SCIDRL', axs[0][3], show=False)
plot_umap_all(SCIDRL, 'batch', 'SCIDRL', axs[1][3], show=False)
plot_umap_all(SCIDRL, 'condition', 'SCIDRL', axs[2][3], show=False)
plt.savefig(data_dir + 'figures/umap/figure8S.png', dpi=600, format='png', bbox_inches='tight')
