from scidrl.main import *
import argparse
import scanpy as sc
import anndata as ad
from time import time
from memory_profiler import profile

os.environ["CUDA_VISIBLE_DEVICES"] = "1"
my_seed = 1234
os.environ['PYTHONHASHSEED']=str(my_seed)

# # 1. Running 10 sampled data ---------------------------------------------------------
data_dir = '../../datasets/human_lung/'
seed = [0, 1, 4, 6, 7, 8, 10, 11, 100, 111]
for i in seed:
    random.seed(my_seed)
    np.random.seed(my_seed)
    tf.random.set_seed(my_seed)
    session_conf = tf.compat.v1.ConfigProto()
    session_conf.gpu_options.allow_growth = True
    sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)
    tf.compat.v1.keras.backend.set_session(sess)

    d = np.load(data_dir + 'sampling/' + str(i) + '_' + 'lung_raw.npz', allow_pickle=True)

    adata = ad.AnnData(d['X_latent'])
    adata.obs['celltype'] = np.array(d['celltype'])
    adata.obs['batch'] = np.array(d['batch'])
    adata.obs['condition'] = np.array(d['condition'])

    adata.obs['celltype'] = adata.obs['celltype'].astype('category')
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    adata.obs['condition'] = adata.obs['condition'].astype('category')

    data_file = data_dir + 'sampling/' + str(i) + '_' + 'lung_raw.csv'
    meta_file = data_dir + 'sampling/' + str(i) + '_' + 'meta.csv'

    parser = argparse.ArgumentParser()
    parser.add_argument('--zdim', type=int, default=16, help='Dim of embedding.')
    parser.add_argument('--znoise_dim', type=int, default=2, help='Dim of noise embedding.')
    parser.add_argument('--epochs', type=int, default=200, help='Number of epochs to train.')  ### 500
    parser.add_argument('--batch_size', type=int, default=100, help='Size of batches to train.')
    parser.add_argument('--lr', type=float, default=1e-3, help='Initial learning rate.')
    parser.add_argument('--gamma', type=float, default=1, help='Weight of classifier loss.')
    parser.add_argument('--fg_lambda', type=float, default=6,
                        help='Weight of GRL.')  #####!!! if the number of datasets is 2, fg_lambda=1;else fg_lambda=the number of datasets#######
    parser.add_argument('--acts', type=str, default='sigmoid', help='Activity function of classifier and discriminator')
    parser.add_argument('--minmaxscale', type=bool, default=True, help='minmax scaling of data')
    params, unknown = parser.parse_known_args()

    if not os.path.isdir(data_dir + "model/"):
        os.makedirs(data_dir + "model/")
    model_file = data_dir + "model/model_" + str(params.fg_lambda) + '-' + str(params.batch_size) + '-' + str(
        params.epochs) + ".h5"

    scidrl = SCIDRL_train(params, data_file, meta_file)
    loss = scidrl.train()
    embed, correct = scidrl.infer()

    inted = pd.DataFrame(embed)
    adata_inted = ad.AnnData(inted, obs=adata.obs, dtype='float64')
    adata_inted.obsm['X_latent'] = adata_inted.X
    adata_inted.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_inted.obs['batch'] = np.array(adata.obs['batch'])
    adata_inted.obs['condition'] = np.array(adata.obs['condition'])

    adata_inted.write_h5ad(data_dir + 'sampling/' + str(i) + '_' + 'lung_scidrl.h5ad')

# # 2. Running Complete data by scidrl -------------------------------------------------------
@profile
def my_func(params, data_file, meta_file):
    start = time()
    scidrl = SCIDRL_train(params, data_file, meta_file)
    loss = scidrl.train()
    embed, correct = scidrl.infer()
    end = time()
    print('elapsed{:.2f} seconds'.format(end - start))
    return embed

if __name__ == '__main__':
    random.seed(my_seed)
    np.random.seed(my_seed)
    tf.random.set_seed(my_seed)
    session_conf = tf.compat.v1.ConfigProto()
    session_conf.gpu_options.allow_growth = True
    sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)
    tf.compat.v1.keras.backend.set_session(sess)

    data_file = data_dir+'lung_raw.csv'
    meta_file = data_dir+'meta.csv'
    adata = sc.read_h5ad(data_dir + 'adata_lung.h5ad')

    parser = argparse.ArgumentParser()
    parser.add_argument('--zdim', type=int, default=16, help='Dim of embedding.')
    parser.add_argument('--znoise_dim', type=int, default=2, help='Dim of noise embedding.')
    parser.add_argument('--epochs', type=int, default=200, help='Number of epochs to train.') ### 500
    parser.add_argument('--batch_size', type=int, default=100, help='Size of batches to train.')
    parser.add_argument('--lr', type=float, default=1e-3, help='Initial learning rate.')
    parser.add_argument('--gamma', type=float, default=1, help='Weight of classifier loss.')
    parser.add_argument('--fg_lambda', type=float, default=6, help='Weight of GRL.')#####!!! if the number of datasets is 2, fg_lambda=1;else fg_lambda=the number of datasets#######
    parser.add_argument('--acts', type=str, default='sigmoid', help='Activity function of classifier and discriminator')
    parser.add_argument('--minmaxscale', type=bool, default=True, help='minmax scaling of data')
    params,unknown=parser.parse_known_args()

    if not os.path.isdir(data_dir+"model/"):
        os.makedirs(data_dir+"model/")
    model_file = data_dir+"model/model_"+str(params.fg_lambda)+'-'+str(params.batch_size)+'-'+str(params.epochs)+".h5"

    embed = my_func(params, data_file, meta_file)

    inted = pd.DataFrame(embed)
    adata_corrd = ad.AnnData(inted, obs=adata.obs, dtype='float64')
    adata_corrd.obsm['X_latent'] = adata_corrd.X
    adata_corrd.obs['celltype'] = np.array(adata.obs['celltype'])
    adata_corrd.obs['batch'] = np.array(adata.obs['batch'])
    adata_corrd.obs['condition'] = np.array(adata.obs['age'])

    adata_corrd.write_h5ad(data_dir + "lung_scidrl.h5ad")



