import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

col1 = ["#E64B35CC", "#0072B5CC", "#00A087CC", "#3C5488CC", "#F39B7FCC", "#F7DC05FF", "#FD7446E5",
       "#8491B4CC", "#7E6148CC", "#B09C85CC", "#E18727CC", "#FFDC91E5", "#6A6599E5", "#9467BDB2",
       "#FFFFFFFF", "#0000FFFF", "#FF0000FF", "#00FF00FF", "#000033FF", "#FF00B6FF", "#005300FF", "#FFD300FF",
       "#009FFFFF", "#9A4D42FF", "#00FFBEFF", "#783FC1FF", "#1F9698FF", "#FFACFDFF", "#B1CC71FF", "#F1085CFF",
       "#FE8F42FF", "#DD00FFFF", "#201A01FF", "#720055FF", "#766C95FF", "#02AD24FF", "#C8FF00FF", "#886C00FF",
       "#FFB79FFF", "#858567FF", "#A10300FF", "#14F9FFFF", "#00479EFF", "#DC5E93FF", "#93D4FFFF", "#004CFFFF"]
col2 = ["#FF7F00", "#7570B3", "#B2DF8A", "#E31A1C"]

def plot_umap_all(adata, name, ax, show, data=None):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    sc.pp.neighbors(adata, use_rep='X_latent')
    sc.tl.umap(adata)
    if name == 'celltype':
        sc.pl.umap(adata, color=['celltype'], ax=ax, show=show, legend_loc='None')
        ax.set_title(data, fontsize=16, family='Arial')
    elif name == 'batch':
        sc.pl.umap(adata, color=['batch'], ax=ax, show=show, legend_loc='None', palette=col1)

    x_axis = ax.get_xaxis()
    y_axis = ax.get_yaxis()
    # 隐藏x轴和y轴的标签
    x_axis.set_visible(False)
    y_axis.set_visible(False)
# %%
fig = plt.figure(figsize=(20, 6.5))
labels = ['scINSIGHT', 'scDisco']

sub_figs = fig.subfigures(2, 1)
axs = []
for i, sub_fig in enumerate(sub_figs):
    axs.append(sub_fig.subplots(1, 5))
    sub_fig.supylabel(labels[i], x=0.09, fontsize=16)
axs = np.array(axs)

# %%
simulate = sc.read_h5ad('../../datasets/simulate/simulate_scINSIGHT.h5ad')
pancreas2 = sc.read_h5ad('../../datasets/pancreas2/pancreas2_scINSIGHT.h5ad')
lung = sc.read_h5ad('../../datasets/lung/lung_scINSIGHT.h5ad')
pancreas = sc.read_h5ad('../../datasets/pancreas/pancreas_scINSIGHT.h5ad')
mouse_mucosa = sc.read_h5ad('../../datasets/mouse_mucosa/mouse_mucosa_scINSIGHT.h5ad')

simulate_X = sc.read_h5ad('../../datasets/simulate/simulate_scDisco.h5ad')
pancreas2_X = sc.read_h5ad('../../datasets/pancreas2/pancreas2_scDisco.h5ad')
lung_X = sc.read_h5ad('../../datasets/lung/lung_scDisco.h5ad')
pancreas_X = sc.read_h5ad('../../datasets/pancreas/pancreas_scDisco.h5ad')
mouse_mucosa_X = sc.read_h5ad('../../datasets/mouse_mucosa/mouse_mucosa_scDisco.h5ad')

# %%
plot_umap_all(simulate, 'celltype', axs[0][0], show=False, data='simulated data')
plot_umap_all(simulate_X, 'celltype', axs[1][0], show=False)
plot_umap_all(pancreas2, 'celltype', axs[0][1], show=False, data='human pancreas2')
plot_umap_all(pancreas2_X, 'celltype', axs[1][1], show=False)
plot_umap_all(lung, 'celltype', axs[0][2], show=False, data='human lung')
plot_umap_all(lung_X, 'celltype', axs[1][2], show=False)
plot_umap_all(pancreas, 'celltype', axs[0][3], show=False, data='human pancreas')
plot_umap_all(pancreas_X, 'celltype', axs[1][3], show=False)
plot_umap_all(mouse_mucosa, 'celltype', axs[0][4], show=False, data='mouse mucosa')
plot_umap_all(mouse_mucosa_X, 'celltype', axs[1][4], show=False)


plt.savefig('../../datasets/scINSIGHT_to_scDisco/figures/celltype_umap_scINSIGHT_to_scDisco.png',
            dpi=600, format='png', bbox_inches='tight')

# %%
plot_umap_all(simulate, 'batch', axs[0][0], show=False, data='simulated data')
plot_umap_all(simulate_X, 'batch', axs[1][0], show=False)
plot_umap_all(pancreas2, 'batch', axs[0][1], show=False, data='human pancreas2')
plot_umap_all(pancreas2_X, 'batch', axs[1][1], show=False)
plot_umap_all(lung, 'batch', axs[0][2], show=False, data='human lung')
plot_umap_all(lung_X, 'batch', axs[1][2], show=False)
plot_umap_all(pancreas, 'batch', axs[0][3], show=False, data='human pancreas')
plot_umap_all(pancreas_X, 'batch', axs[1][3], show=False)
plot_umap_all(mouse_mucosa, 'batch', axs[0][4], show=False, data='mouse mucosa')
plot_umap_all(mouse_mucosa_X, 'batch', axs[1][4], show=False)


plt.savefig('../../datasets/scINSIGHT_to_scDisco/figures/batch_umap_scINSIGHT_to_scDisco.png',
            dpi=600, format='png', bbox_inches='tight')