import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# %%
datasets_name = ["pancreas2", "lung", 'pancreas', 'mouse_mucosa', 'epithelium', 'ductal']
datasets = ['human pancreas2', 'human lung', 'human pancreas',
            'mouse mucosa', 'human epithelium', 'human ductal']

def plot_comp(data_dir, i, ax):
    plt.subplot(ax)
    data = pd.read_csv(data_dir + 'metircs_' + datasets_name[i] + '_ae' + ".csv")
    df1 = pd.DataFrame(data, columns=['ARI'], dtype=float)
    df2 = pd.DataFrame(data, columns=['F1_Score'], dtype=float)

    ARI = np.array(df1).squeeze()
    F_Score = np.array(df2).squeeze()
    ARI = ARI.tolist()
    F_Score = F_Score.tolist()

    VAE = []
    AE = []
    VAE.append(ARI[0])
    AE.append(ARI[1])
    VAE.append(F_Score[0])
    AE.append(F_Score[1])

    vae_1 = VAE
    ae_1 = AE

    plt.ylim(0, 1)
    plt.xticks(fontsize=7)
    plt.tick_params(labelsize=24)

    bar_width = 0.1
    index = np.array([0, 0.3])

    plt.bar(index - 0.5 * bar_width, vae_1, bar_width, alpha=1, color='#FB8930', label='VAE')
    plt.bar(index + 0.5 * bar_width, ae_1, bar_width, alpha=1, color='#5C89A0', label='AE')
    x_labels = ['ARI', 'F1 Score']
    plt.xticks(index + bar_width/50, x_labels, fontsize=24, family='Arial')  # index+bar_width/2 to centre the label
    plt.ylim(0, 1)
    plt.tight_layout()

# %%
fig = plt.figure(figsize=(15, 5))
plot_comp('../../datasets/pancreas2/', 0, 131)
plot_comp('../../datasets/lung/', 1, 132)
plot_comp('../../datasets/pancreas/', 2, 133)
plt.savefig('figures/figure13C.svg', dpi=300, format='svg', bbox_inches='tight')
plot_comp('../../datasets/mouse_mucosa/', 3, 131)
plot_comp('../../datasets/epithelium/', 4, 132)
plot_comp('../../datasets/ductal/', 5, 133)
plt.savefig('figures/figure13D.svg', dpi=300, format='svg', bbox_inches='tight')






