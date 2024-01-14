import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')
os.environ['CUDA_VISIBLE_DEVICES'] = '1'
plt.rcParams['font.sans-serif'] = ['Arial']


df_t = pd.DataFrame()
df_m = pd.DataFrame()

datasets = ['human pancreas2', 'human lung', 'human pancreas',
            'mouse mucosa', 'human epithelium', 'human ductal']

labels = ["Seurat", "Harmony", 'Scanorama', 'DESC', 'scVI', 'cell Blast', 'SCIDRL','scDisInFact',
          'scINSIGHT', "scDisco"]
methods = ["simulate_seurat", "simulate_harmony", 'simulate_scanorama', 'simulate_desc',
           'simulate_scVI', 'simulate_cell_blast', 'simulate_scidrl', 'simulate_scDisinfact',
           "simulate_scINSIGHT", "simulate_scDisco"]
methods_name = ["Seurat", "Harmony", 'Scanorama', 'DESC', 'scVI', 'Cell BLAST', 'SCIDRL',
                'scDisInFact', 'scINSIGHT', "scDisco"]
colors = ['#D9412B', '#F67948', '#FABB6E', '#F9E07F', '#8f8ce7', '#92C5DE', '#5C90C2', '#3951A2', '#A52A2A', '#E79796']
markers = ['X', 'P', 'd', '^', '1', 'v', '*', 's', '+', 'o']

Seurat_t = [73.83, 38.09, 160.11, 179.01, 2245.13, None]
Harmony_t = [5.95, 5.82, 18.91, 24.53, 55.93, 144.01]
Scanorma_t = [98.10, 30.28, 207.56, 1319.99, 3461.55, 14231.54]
DESC_t = [64.10, 75.96, 199.71, 403.95, 482.22, 1126.53]
scVI_t = [276.89, 405.31, 1019.40, 1374.40, 2483.61, 2415.75]
cell_BLAST_t = [214.48, 99.18, 384.48, 203.65, 981.81, 2627.47]
SCIDRL_t = [1116.68, 5401.86, 2766.30, 5424.46, 11787.49, 19984.86]
scDisInFact_t = [1530.45, 1863.18, 4599.78, 12652.22, 7326.38, 13019.08]
scINSIGHT_t = [3557.21, 3587.29, 9452.02, 9215.27, None, None]
scDisco_t = [25.10, 37.53, 109.70, 145.63, 337.01, 841.05]

Seurat_m = [7078.06, 7133.38, 13071.06, 22274.75, 52420.44, None]
Harmony_m = [649.3, 621.1, 1056.9, 1709.7, 2914.5, 4375.1]
Scanorma_m = [1593.1, 1333.7, 4236.8, 3714.3, 7627.7, 9550.0]
DESC_m = [1221.90, 1181.7, 2185.0, 4009.6, 7682.6, 2646.9]
scVI_m = [1420.90, 1380.2, 2754.7, 5983.9, 10386.7, 7377.5]
cell_BLAST_m = [824.20, 667.4, 1145.6, 2316.4, 3355.0, 1891.6]
SCIDRL_m = [2802.5, 2809.3, 4708.6, 8779.4, 5851.5, 11259.3]
scDisInFact_m = [1365.4, 2948.9, 4014.0, 5479.2, 8891.7, 12250.0]
scINSIGHT_m = [4689.75, 4841.12, 14355.88, 26809.38, None, None]
scDisco_m = [622.2, 618.9, 1001.7, 1728.4, 3322.8, 2122.1]

x_index = np.arange(6)

data1 = list(zip(Seurat_t, Harmony_t, Scanorma_t, DESC_t, scVI_t, cell_BLAST_t, SCIDRL_t,
                 scDisInFact_t, scINSIGHT_t, scDisco_t))
t_df = pd.DataFrame(data1, columns=methods_name, index=datasets)
n, m = t_df.shape

plt.figure(figsize=(10, 6))

plt.rcParams["font.family"] = "Arial"
plt.rcParams.update({'font.size': 28})
plt.rc('legend', fontsize=15)

for i in range(m):
    y = t_df.iloc[:, i]
    x = x_index
    plt.plot(x, y, marker=markers[i], color=colors[i], markersize=8, label=methods_name[i])

plt.xticks(x, datasets, rotation=30, family='Arial', fontsize=20)
plt.yticks(family='Arial', size=15)
plt.legend(loc='upper left', framealpha=0.7)
plt.ylabel("time / s", family='Arial', size=20)
plt.tight_layout()
plt.savefig('figures/time.svg', dpi=300, format='svg', bbox_inches='tight')

data2 = list(zip(Seurat_m, Harmony_m, Scanorma_m, DESC_m, scVI_m, cell_BLAST_m, SCIDRL_m, scDisInFact_m,
                 scINSIGHT_m, scDisco_m))
m_df = pd.DataFrame(data2, columns=methods_name, index=datasets)
n, m = m_df.shape

plt.figure(figsize=(10, 6))

plt.rcParams["font.family"] = "Arial"
plt.rcParams.update({'font.size': 28})
plt.rc('legend', fontsize=15)


for i in range(m):
    y = m_df.iloc[:, i]
    x = x_index
    plt.plot(x, y, marker=markers[i], color=colors[i], markersize=8, label=methods_name[i])

plt.xticks(x, datasets, rotation=30, family='Arial', fontsize=20)
plt.yticks(family='Arial', size=15)
plt.legend(loc='upper left', framealpha=0.7)
plt.ylabel("memory usage / MB", family='Arial', size=20)
plt.tight_layout()
plt.savefig('figures/memory.svg', dpi=300, format='svg', bbox_inches='tight')



