import numpy as np
import scanpy as sc
import pandas as pd
from util.utils import set_seed, calculate_metric, silhouette_coeff_ASW, kmeans
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

def count(data_dir, data_dir_, datasets, method_names):
    for method_name in method_names:
        all_nmi, all_ari = [], []
        all_f1_score = []
        all_basw, all_casw, all_asw = [], [], []
        set_seed(2023)
        for dataset in datasets:
            adata = sc.read_h5ad(data_dir + dataset + '/' + dataset + '_' + method_name + ".h5ad")

            cell_type = adata.obs['celltype'].values
            Label_enc = LabelEncoder()
            cell_type = Label_enc.fit_transform(cell_type)
            batches = np.array(adata.obs['batch'])
            OneHot_enc = OneHotEncoder()
            onehot_code = OneHot_enc.fit_transform(batches.reshape(-1, 1))
            b = onehot_code.toarray()
            b = b.astype(np.float32)

            # calculate nmi, ari
            if method_name == 'scINSIGHT':
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
        A.to_csv(data_dir_ + 'metircs_' + method_name + ".csv")




