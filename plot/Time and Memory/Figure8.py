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
Harmony_t = [279.69, 370.53, 1576.69, 2041.26, 4994.23, 20208.74]
Scanorma_t = [84.70, 87.61, 435.85, 1202.40, 2028.87, 14140.63]
DESC_t = [64.10, 75.96, 199.71, 403.95, 482.22, 1126.53]
scVI_t = [276.89, 405.31, 1019.40, 1374.40, 2483.61, 2415.75]
cell_BLAST_t = [214.48, 99.18, 384.48, 203.65, 981.81, 2627.47]
SCIDRL_t = [962.44, 947.76, 2730.56, 5215.26, 11232.05, 20769.08]
scDisInFact_t = [1530.45, 1863.18, 4599.78, 12652.22, 7326.38, 13019.08]
scINSIGHT_t = [50538.42, 24965.26, 87148.31, 169311.07, None, None]
scDisco_t = [25.10, 37.53, 109.70, 145.63, 337.01, 841.05]

Seurat_m = [7078.06, 7133.38, 13071.06, 22274.75, 52420.44, None]
Harmony_m = [1437.5, 1352.1, 1963.5, 3838.0, 7376.0, 9300.0]
Scanorma_m = [825.40, 798.4, 1683.3, 3755.7, 5225.4,  8253.3]
DESC_m = [1221.90, 1181.7, 2185.0, 4009.6, 7682.6, 2646.9]
scVI_m = [1420.90, 1380.2, 2754.7, 5983.9, 10386.7, 7377.5]
cell_BLAST_m = [824.20, 667.4, 1145.6, 2316.4, 3355.0, 1891.6]
SCIDRL_m = [2358.9, 2531.6, 4245.2, 7691.1, 6998.4, 13042.3]
scDisInFact_m = [1365.4, 2948.9, 4014.0, 5479.2, 8891.7, 12250.0]
scINSIGHT_m = [5056.38, 5251.75, 17201.31, 27644.88, None, None]
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



