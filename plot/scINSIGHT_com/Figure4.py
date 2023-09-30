# %%
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# %%
metrics_name = ['ARI', 'F1_Score']

def plot_comp(data_dir, i, ax, y):
    plt.subplot(ax)
    data = pd.read_csv(data_dir + 'metircs_' + 'scINSIGHT' + ".csv")
    df1 = pd.DataFrame(data, columns=[metrics_name[i]], dtype=float)

    data2 = pd.read_csv(data_dir + 'metircs_' + 'scDisco' + ".csv")
    df2 = pd.DataFrame(data2, columns=[metrics_name[i]], dtype=float)


    scINSIGHT = np.array(df1).squeeze()
    scDisco = np.array(df2).squeeze()
    scINSIGHT = scINSIGHT.tolist()
    scDisco = scDisco.tolist()

    plt.xticks(fontsize=15)
    plt.tick_params(labelsize=15)

    bar_width = 0.35
    index = np.arange(5)

    plt.bar(index - bar_width / 2, scINSIGHT, bar_width, alpha=0.8, color='#FB8930', label='scINSIGHT')
    plt.bar(index + bar_width / 2, scDisco, bar_width, alpha=1, color='#5C89A0', label='scDisco')

    x_labels = ['simulated data', 'human pancreas2', 'human lung', 'human pancreas', 'mouse mucosa']
    plt.xticks(index, x_labels, fontsize=20, rotation=45, family='Arial', ha='right')
    plt.ylim(0, 1)
    plt.ylabel(y, fontsize=24, family='Arial',  labelpad=20)
    plt.legend(loc='best')


# %%
fig = plt.figure(figsize=(10, 5))
data_dir = '../../datasets/scINSIGHT_to_scDisco/'
plot_comp(data_dir, 0, 121, 'ARI')
plot_comp(data_dir, 1, 122, 'F1 Score')
plt.tight_layout()
plt.savefig(data_dir + 'figures/figure12.svg', dpi=300, format='svg', bbox_inches='tight')






