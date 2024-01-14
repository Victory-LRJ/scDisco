from scidrl.main import *
import argparse
import scanpy as sc
import anndata as ad
from time import time
from memory_profiler import profile

os.environ["CUDA_VISIBLE_DEVICES"] = "1"
my_seed = 1234
os.environ['PYTHONHASHSEED']=str(my_seed)

data_dir = '../../datasets/human_ductal/'

# # 1. Running Complete data by scidrl -------------------------------------------------------
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

    data_file = data_dir+'ductal_raw_scidrl.csv'
    meta_file = data_dir+'meta_scidrl.csv'
    adata = sc.read_h5ad(data_dir + 'human_ductal_scidrl.h5ad')

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
    adata_corrd.obs['condition'] = np.array(adata.obs['disease'])

    adata_corrd.write_h5ad(data_dir + "ductal_scidrl.h5ad")



