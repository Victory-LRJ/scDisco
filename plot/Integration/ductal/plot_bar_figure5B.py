# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# %%  read data
data_dir = '../../../datasets/ductal/'
m = ["NMI", "ARI", "F1_Score", "SAS", 'BASW', 'CASW']
method_name = 'ductal_'

data = pd.read_csv(data_dir + 'metircs_' + method_name + ".csv")
df1 = pd.DataFrame(data, columns=['NMI'], dtype=float)
df2 = pd.DataFrame(data, columns=['ARI'], dtype=float)
df3 = pd.DataFrame(data, columns=['F1_Score'], dtype=float)
df4 = pd.DataFrame(data, columns=['SAS'], dtype=float)
df5 = pd.DataFrame(data, columns=['BASW'], dtype=float)
df6 = pd.DataFrame(data, columns=['CASW'], dtype=float)
df7 = pd.DataFrame(data, columns=['ASW'], dtype=float)

NMI = np.array(df1).squeeze().T
ARI = np.array(df2).squeeze().T
F1_Score = np.array(df3).squeeze().T
SAS = np.array(df4).squeeze().T
BASW = np.array(df5).squeeze().T
CASW = np.array(df6).squeeze().T
ASW = np.array(df7).squeeze().T

# %%  plot bar
def plot_bar(x, y, label, color, metrics):
    plt.bar(x, y, label=label, color=color, width=0.6)
    plt.xticks(rotation=45, fontsize=48, family='Arial', ha='right')
    plt.yticks(fontsize=36, family='Arial', ha='right')
    plt.ylabel(metrics, fontsize=60, family='Arial', labelpad=20)
    plt.rcParams.update({'font.size': 36})


color = ['lawngreen', 'peru', 'cornflowerblue', 'limegreen', 'rosybrown', 'gold', 'deepskyblue', 'orange']
x = ['SCIDRL', 'scDisInFact', "Harmony", "Scanorama", 'Cell BLAST', 'DESC', 'scVI', "scDisco"]

y_NMI = NMI
y_ARI = ARI
y_F1_Score = F1_Score
y_SAS = SAS
y_BASW = BASW
y_CASW = CASW
y_ASW = ASW


fig1 = plt.figure(figsize=[21, 12])
plot_bar(x, y_ARI, x, color, 'ARI')
plt.savefig(data_dir + 'figures/evaluation/figure8B1.svg', dpi=300, format='svg', bbox_inches='tight')

fig2 = plt.figure(figsize=[21, 12])
plot_bar(x, y_F1_Score, x, color, 'F1 Score')
plt.savefig(data_dir + 'figures/evaluation/figure8B2.svg', dpi=300, format='svg', bbox_inches='tight')

fig3 = plt.figure(figsize=[21, 12])
plot_bar(x, y_NMI, x, color, 'NMI')
plt.savefig(data_dir + 'figures/evaluation/ductal_nmi.svg', dpi=300, format='svg', bbox_inches='tight')

fig4 = plt.figure(figsize=[21, 12])
plot_bar(x, y_BASW, x, color, 'BASW')
plt.savefig(data_dir + 'figures/evaluation/ductal_basw.svg', dpi=300, format='svg', bbox_inches='tight')

fig5 = plt.figure(figsize=[21, 12])
plot_bar(x, y_CASW, x, color, 'CASW')
plt.savefig(data_dir + 'figures/evaluation/ductal_casw.svg', dpi=300, format='svg', bbox_inches='tight')

fig6 = plt.figure(figsize=[21, 12])
plot_bar(x, y_ASW, x, color, 'ASW')
plt.savefig(data_dir + 'figures/evaluation/ductal_asw.svg', dpi=300, format='svg', bbox_inches='tight')



