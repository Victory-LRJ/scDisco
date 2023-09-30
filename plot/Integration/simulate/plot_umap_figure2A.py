import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# %%
col1 = ["#E64B35CC", "#0072B5CC", "#00A087CC", "#3C5488CC", "#F39B7FCC", "#F7DC05FF", "#FD7446E5",
       "#8491B4CC", "#7E6148CC", "#B09C85CC", "#E18727CC", "#FFDC91E5", "#6A6599E5", "#9467BDB2",
       "#FFFFFFFF", "#0000FFFF", "#FF0000FF", "#00FF00FF", "#000033FF", "#FF00B6FF", "#005300FF", "#FFD300FF",
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
data_dir = '../../../datasets/simulate/'
fig = plt.figure(figsize=(20, 10))
datasets = ['by cell type', 'by batch', 'by condition']
sub_figs = fig.subfigures(3, 1)
axs = []
for i, sub_fig in enumerate(sub_figs):
    axs.append(sub_fig.subplots(1, 5))
    sub_fig.supylabel(datasets[i], x=0.08, fontsize=24, family='Arial')
axs = np.array(axs)

# %%
Raw = sc.read_h5ad(data_dir + 'sim_count.h5ad')
Raw.obsm['X_latent'] = Raw.X
Harmony = sc.read_h5ad(data_dir + 'simulate_harmony.h5ad')
DESC = sc.read_h5ad(data_dir + 'simulate_desc.h5ad')
scVI = sc.read_h5ad(data_dir + 'simulate_scVI.h5ad')
Seurat = sc.read_h5ad(data_dir + 'simulate_seurat.h5ad')
Scanorama = sc.read_h5ad(data_dir + 'simulate_scanorama.h5ad')
cell_BLAST = sc.read_h5ad(data_dir + 'simulate_cell_blast.h5ad')
SCIDRL = sc.read_h5ad(data_dir + 'simulate_scidrl.h5ad')
scDisInFact = sc.read_h5ad(data_dir + 'simulate_scDisinfact.h5ad')
scDisco = sc.read_h5ad(data_dir + 'simulate_scDisco.h5ad')
# %%
plot_umap_all(Raw, 'celltype', 'Raw', axs[0][0], show=False)
plot_umap_all(Raw, 'batch', 'Raw', axs[1][0], show=False)
plot_umap_all(Raw, 'condition', 'Raw', axs[2][0], show=False)
plot_umap_all(scDisco, 'celltype', 'scDisco', axs[0][1], show=False)
plot_umap_all(scDisco, 'batch', 'scDisco', axs[1][1], show=False)
plot_umap_all(scDisco, 'condition', 'scDisco', axs[2][1], show=False)
plot_umap_all(Harmony, 'celltype', 'Harmony', axs[0][2], show=False)
plot_umap_all(Harmony, 'batch', 'Harmony', axs[1][2], show=False)
plot_umap_all(Harmony, 'condition', 'Harmony', axs[2][2], show=False)
plot_umap_all(scVI, 'celltype', 'scVI', axs[0][3], show=False)
plot_umap_all(scVI, 'batch', 'scVI', axs[1][3], show=False)
plot_umap_all(scVI, 'condition', 'scVI', axs[2][3], show=False)
plot_umap_all(Seurat, 'celltype', 'Seurat', axs[0][4], show=False)
plot_umap_all(Seurat, 'batch', 'Seurat', axs[1][4], show=False)
plot_umap_all(Seurat, 'condition', 'Seurat', axs[2][4], show=False)
plt.savefig(data_dir + 'figures/umap/figure2A.png', dpi=600, format='png')

sc.pl.umap(Raw, color=['celltype'])
plt.savefig(data_dir + 'figures/umap/plot_celltype.png', dpi=600, format='png', bbox_inches='tight')
sc.pl.umap(Raw, color=['batch'], palette=col1)
plt.savefig(data_dir + 'figures/umap/plot_batch.png', dpi=600, format='png', bbox_inches='tight')
sc.pl.umap(Raw, color=['condition'], palette=col2)
plt.savefig(data_dir + 'figures/umap/plot_condition.png', dpi=600, format='png', bbox_inches='tight')