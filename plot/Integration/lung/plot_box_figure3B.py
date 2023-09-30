import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')

# %% read data
data_dir = '../../../datasets/lung/'
methods = ["harmony", 'scanorama', 'scDisinfact', 'scidrl', "seurat", 'scVI', 'cell_blast', "scDisco", 'desc']
m = ["NMI", "ARI", "F1_Score", "SAS", 'BASW', 'CASW']

NMI, ARI, F_Score, SAS, BASW, CASW, ASW = [], [], [], [], [], [], []
for method in methods:
    data = pd.read_csv(data_dir + 'metircs_' + method + ".csv")
    df1 = pd.DataFrame(data, columns=['NMI'], dtype=float)
    NMI.append(np.array(df1))
    df2 = pd.DataFrame(data, columns=['ARI'], dtype=float)
    ARI.append(np.array(df2))
    df3 = pd.DataFrame(data, columns=['F1_Score'], dtype=float)
    F_Score.append(np.array(df3))
    df4 = pd.DataFrame(data, columns=['SAS'], dtype=float)
    SAS.append(np.array(df4))
    df5 = pd.DataFrame(data, columns=['BASW'], dtype=float)
    BASW.append(np.array(df5))
    df6 = pd.DataFrame(data, columns=['CASW'], dtype=float)
    CASW.append(np.array(df6))
    df7 = pd.DataFrame(data, columns=['ASW'], dtype=float)
    ASW.append(np.array(df7))
NMI = np.array(NMI).squeeze().T
ARI = np.array(ARI).squeeze().T
F1_Score = np.array(F_Score).squeeze().T
SAS = np.array(SAS).squeeze().T
BASW = np.array(BASW).squeeze().T
CASW = np.array(CASW).squeeze().T
ASW = np.array(ASW).squeeze().T

colors = ['pink', 'lightgray', 'paleturquoise', 'darkseagreen', 'lightsalmon', 'powderblue', 'thistle', 'gold', 'mistyrose']
labels = ["Harmony", 'Scanorama', 'scDisInFact', 'SCIDRL', "Seurat", 'scVI', 'Cell BLAST', "scDisco", 'DESC']
# %%  box plot ARI
fig, axes = plt.subplots(figsize=(6, 3))
bplot = axes.boxplot(ARI, showmeans=False, widths=0.5, notch=False, vert=True, patch_artist=True, showfliers=False,
                              medianprops={'linestyle': '-', 'color': 'black'}, whiskerprops={'linestyle': '--'},
                              labels=labels)
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
axes.set_xticklabels(labels, fontsize=15, rotation=45, family='Arial', ha='right')
axes.set_ylabel('ARI', fontsize=24, family='Arial', labelpad=20)
axes.set_ylim(0, 1)
fig.tight_layout()
plt.savefig(data_dir + 'figures/evaluation/figure4B1.svg', dpi=300, format='svg', bbox_inches='tight')

# %% box plot F1 Score
fig, axes = plt.subplots(figsize=(6, 3))
bplot = axes.boxplot(F1_Score, showmeans=False, widths=0.5, notch=False, vert=True, patch_artist=True, showfliers=False,
                              medianprops={'linestyle': '-', 'color': 'black'}, whiskerprops={'linestyle': '--'},
                              labels=labels)
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
axes.set_xticklabels(labels, fontsize=15, rotation=45, family='Arial', ha='right')
axes.set_ylabel('F1 Score', fontsize=24, family='Arial', labelpad=20)
axes.set_ylim(0, 1)
fig.tight_layout()
plt.savefig(data_dir + 'figures/evaluation/figure4B2.svg', dpi=300, format='svg', bbox_inches='tight')

# %% box plot NMI
fig, axes = plt.subplots(figsize=(6, 3))
bplot = axes.boxplot(NMI, showmeans=False, widths=0.5, notch=False, vert=True, patch_artist=True, showfliers=False,
                              medianprops={'linestyle': '-', 'color': 'black'}, whiskerprops={'linestyle': '--'},
                              labels=labels)
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
axes.set_xticklabels(labels, fontsize=15, rotation=45, family='Arial', ha='right')
axes.set_ylabel('NMI', fontsize=24, family='Arial', labelpad=20)
axes.set_ylim(0, 1)
fig.tight_layout()
plt.savefig(data_dir + 'figures/evaluation/lung_nmi.svg', dpi=300, format='svg', bbox_inches='tight')

# %% box plot SAS
fig, axes = plt.subplots(figsize=(6, 3))
bplot = axes.boxplot(SAS, showmeans=False, widths=0.5, notch=False, vert=True, patch_artist=True, showfliers=False,
                              medianprops={'linestyle': '-', 'color': 'black'}, whiskerprops={'linestyle': '--'},
                              labels=labels)
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
axes.set_xticklabels(labels, fontsize=15, rotation=45, family='Arial', ha='right')
axes.set_ylabel('SAS', fontsize=24, family='Arial', labelpad=20)
axes.set_ylim(0, 1)
fig.tight_layout()
plt.savefig(data_dir + 'figures/evaluation/lung_SAS.svg', dpi=300, format='svg', bbox_inches='tight')

# %% box plot BASW
fig, axes = plt.subplots(figsize=(6, 3))
bplot = axes.boxplot(BASW, showmeans=False, widths=0.5, notch=False, vert=True, patch_artist=True, showfliers=False,
                              medianprops={'linestyle': '-', 'color': 'black'}, whiskerprops={'linestyle': '--'},
                              labels=labels)
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
axes.set_xticklabels(labels, fontsize=15, rotation=45, family='Arial', ha='right')
axes.set_ylabel('BASW', fontsize=24, family='Arial', labelpad=20)
axes.set_ylim(0, 1)
fig.tight_layout()
plt.savefig(data_dir + 'figures/evaluation/lung_basw.svg', dpi=300, format='svg', bbox_inches='tight')

# %% box plot CASW
fig, axes = plt.subplots(figsize=(6, 3))
bplot = axes.boxplot(CASW, showmeans=False, widths=0.5, notch=False, vert=True, patch_artist=True, showfliers=False,
                              medianprops={'linestyle': '-', 'color': 'black'}, whiskerprops={'linestyle': '--'},
                              labels=labels)
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
axes.set_xticklabels(labels, fontsize=15, rotation=45, family='Arial', ha='right')
axes.set_ylabel('CASW', fontsize=24, family='Arial', labelpad=20)
axes.set_ylim(0, 1)
fig.tight_layout()
plt.savefig(data_dir + 'figures/evaluation/lung_casw.svg', dpi=300, format='svg', bbox_inches='tight')

# %%  box plot ASW
fig, axes = plt.subplots(figsize=(6, 3))
bplot = axes.boxplot(ASW, showmeans=False, widths=0.5, notch=False, vert=True, patch_artist=True, showfliers=False,
                              medianprops={'linestyle': '-', 'color': 'black'}, whiskerprops={'linestyle': '--'},
                              labels=labels)
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
axes.set_xticklabels(labels, fontsize=15, rotation=45, family='Arial', ha='right')
axes.set_ylabel('ASW', fontsize=24, family='Arial', labelpad=20)
axes.set_ylim(0, 1)
fig.tight_layout()
plt.savefig(data_dir + 'figures/evaluation/lung_asw.svg', dpi=300, format='svg', bbox_inches='tight')
