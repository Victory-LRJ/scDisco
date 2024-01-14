import scanpy as sc
import numpy as np
import anndata as ad
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
def plot_pca_all(adata, name, method, ax, show):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    adata_pca = ad.AnnData(adata.obsm['X_latent'])
    adata_pca.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_pca.obs['batch'] = np.array(adata.obs['batch'])
    adata_pca.obs['condition'] = np.array(adata.obs['condition'])
    if method not in ['DESC']:
        sc.tl.pca(adata_pca, n_comps=2, svd_solver='arpack')
    elif method in ['DESC']:
        adata_pca.obsm['X_pca'] = adata.obsm['X_latent']
    if name == 'celltype':
        sc.pl.pca(adata_pca, color=['celltype'], ax=ax, show=show, legend_loc='None')
        ax.set_title(method, fontsize=24, family='Arial')
    elif name == 'batch':
        sc.pl.pca(adata_pca, color=['batch'], ax=ax, show=show, legend_loc='None', palette=col1)
    elif name == 'condition':
        sc.pl.pca(adata_pca, color=['condition'], ax=ax, show=show, legend_loc='None', palette=col2)

    x_axis = ax.get_xaxis()
    y_axis = ax.get_yaxis()
    x_axis.set_visible(False)
    y_axis.set_visible(False)

# %%
data_dir = 'C://0scDisco/datasets/simulate/'
fig = plt.figure(figsize=(22, 10))
datasets = ['by cell type', 'by batch', 'by condition']
sub_figs = fig.subfigures(3, 1)
axs = []
for i, sub_fig in enumerate(sub_figs):
    axs.append(sub_fig.subplots(1, 6))
    sub_fig.supylabel(datasets[i], x=0.08, fontsize=24)
axs = np.array(axs)

# %%
Raw = sc.read_h5ad(data_dir + 'sim_count.h5ad')
Raw.obsm['X_latent'] = Raw.X
Seurat = sc.read_h5ad(data_dir + 'simulate_seurat.h5ad')
Harmony = sc.read_h5ad(data_dir + 'simulate_harmony.h5ad')
Scanorama = sc.read_h5ad(data_dir + 'simulate_scanorama.h5ad')
DESC = sc.read_h5ad(data_dir + 'simulate_desc.h5ad')
scVI = sc.read_h5ad(data_dir + 'simulate_scVI.h5ad')
cell_BLAST = sc.read_h5ad(data_dir + 'simulate_cell_blast.h5ad')
SCIDRL = sc.read_h5ad(data_dir + 'simulate_scidrl.h5ad')
scDisInFact = sc.read_h5ad(data_dir + 'simulate_scDisinfact.h5ad')
scINSIGHT = sc.read_h5ad(data_dir + 'simulate_scINSIGHT.h5ad')
scDisco = sc.read_h5ad(data_dir + 'simulate_scDisco.h5ad')
# %%
plot_pca_all(cell_BLAST, 'celltype', 'Cell BLAST', axs[0][0], show=False)
plot_pca_all(cell_BLAST, 'batch', 'Cell BLAST', axs[1][0], show=False)
plot_pca_all(cell_BLAST, 'condition', 'Cell BLAST', axs[2][0], show=False)
plot_pca_all(SCIDRL, 'celltype', 'SCIDRL', axs[0][1], show=False)
plot_pca_all(SCIDRL, 'batch', 'SCIDRL', axs[1][1], show=False)
plot_pca_all(SCIDRL, 'condition', 'SCIDRL', axs[2][1], show=False)
plot_pca_all(DESC, 'celltype', 'DESC', axs[0][2], show=False)
plot_pca_all(DESC, 'batch', 'DESC', axs[1][2], show=False)
plot_pca_all(DESC, 'condition', 'DESC', axs[2][2], show=False)
plot_pca_all(Scanorama, 'celltype', 'Scanorama', axs[0][3], show=False)
plot_pca_all(Scanorama, 'batch', 'Scanorama', axs[1][3], show=False)
plot_pca_all(Scanorama, 'condition', 'Scanorama', axs[2][3], show=False)
plot_pca_all(scDisInFact, 'celltype', 'scDisInFact', axs[0][4], show=False)
plot_pca_all(scDisInFact, 'batch', 'scDisInFact', axs[1][4], show=False)
plot_pca_all(scDisInFact, 'condition', 'scDisInFact', axs[2][4], show=False)
plt.savefig(data_dir + 'figures/pca2.png', dpi=600, format='png')