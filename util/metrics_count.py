import numpy as np
import scanpy as sc
import pandas as pd
from util.utils import calculate_metric, silhouette_coeff_ASW, kmeans
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

def count(data_dir, method_name, methods):
    seed = [0, 1, 4, 6, 7, 8, 10, 11, 100, 111]
    for method in methods:
        all_nmi, all_ari = [], []
        all_f1_score = []
        all_basw, all_casw, all_asw = [], [], []
        for i in seed:
            adata = sc.read_h5ad(data_dir + str(i) + '_' + method_name + method + ".h5ad")

            cell_type = adata.obs['celltype'].values
            Label_enc = LabelEncoder()
            cell_type = Label_enc.fit_transform(cell_type)
            batches = np.array(adata.obs['batch'])
            OneHot_enc = OneHotEncoder()
            onehot_code = OneHot_enc.fit_transform(batches.reshape(-1, 1))
            b = onehot_code.toarray()
            b = b.astype(np.float32)

            # calculate acc, nmi, ari
            if method in ['seurat', 'desc', 'scINSIGHT']:
                y_pred_k = np.array(adata.obs['clustering'])
                nmi_k, ari_k = calculate_metric(cell_type, y_pred_k)
            else:
                n_clusters = len(np.unique(np.array(adata.obs['celltype'])))
                adata = kmeans(adata, n_clusters, use_rep='X_latent')
                y_pred_k = np.array(adata.obs['kmeans'])
                nmi_k, ari_k = calculate_metric(cell_type, y_pred_k)
                print(ari_k)

            # calculate basw, casw, f_score
            F1_Score, BASW, CASW, ASW = silhouette_coeff_ASW(adata, percent_extract=0.9)

            all_nmi.append(nmi_k)
            all_ari.append(ari_k)
            all_f1_score.append(F1_Score)
            all_basw.append(BASW)
            all_casw.append(CASW)
            all_asw.append(ASW)


        A = pd.DataFrame({'NMI': np.array(all_nmi), "ARI": np.array(all_ari),
                          'F1_Score': np.array(all_f1_score),
                          'BASW': np.array(all_basw), 'CASW': np.array(all_casw), 'ASW': np.array(all_asw)
                          })
        A.to_csv(data_dir + 'metircs_' + method + ".csv")




