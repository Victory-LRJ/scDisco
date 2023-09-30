from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from model.scDisco import scDisco
from model.network import *
from torch import nn, optim
from torch.optim.lr_scheduler import StepLR
from torch.utils.data import DataLoader
from util.datasets import Batch_Dataset
from util.losses import ELOBkldloss
from pylab import *
from util.utils import calculate_metric, kmeans
from util.data_utils import getdims


def run_scDisco(adata,
                  n_clusters: int = None,
                  batch_key: str = "batch", aspect_key: list = None, cl_type: str = None,
                  EncLayer: list = [512, 256], DecLayer: list = [256, 512],
                  cond_dim: int = 8, n_epochs: int = None, batch_size: int = 256,
                  lr_c: float = 0.0001, lr_m: int = 0.0001, lda: float = 0.001, lda_c: float = 0.001,
                  model_dir: str = None,
                  ):
    device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")
    ####################   Prepare data for training   ######################################################################
    raw_mat = adata.raw.X
    exp_mat = adata.X
    bio_z_dim = getdims(exp_mat)

    Label_enc = LabelEncoder()
    batch_labels = Label_enc.fit_transform(adata.obs[batch_key])

    batches = adata.obs[batch_key].values.to_numpy()
    n_batches = len(np.unique(batches))

    OneHot_enc = OneHotEncoder()
    onehot_code = OneHot_enc.fit_transform(batches.reshape(-1, 1))
    batches_onehot = onehot_code.toarray()
    batch_onehot = batches_onehot.astype(np.float32)

    aspects = {aspect: adata.obs[aspect].values.to_numpy() for aspect in aspect_key}
    n_aspects = len(aspect_key)
    n_conditions = {aspect: len(np.unique(conditions)) for (aspect, conditions) in aspects.items()}

    condition_types = {}
    condition_labels = {}
    for (aspect, conditions) in aspects.items():
        Label_enc = LabelEncoder()
        condition_types[aspect] = Label_enc.fit_transform(conditions)
        condition_labels[aspect] = torch.tensor(condition_types[aspect])

    print(' '.join(["Number of batches is {},".format(n_batches)]
                   + ["number of conditions for aspect {} is {}".format(aspect, n_cond)
                      for (aspect, n_cond) in n_conditions.items()]))

    cell_type = adata.obs[cl_type].values if cl_type else None
    Label_enc = LabelEncoder()
    cell_type = Label_enc.fit_transform(cell_type)

    if n_clusters is None:
        if cell_type is not None:
            n_clusters = len(np.unique(cell_type))
        else:
            raise Exception('Please input number of clusters or set cell type information in adata.obs[cl_type]')

    size_factor = 'size_factor'
    if size_factor not in adata.obs:
        size_factor = np.ones((exp_mat.shape[0], 1), dtype=np.float32)
    else:
        size_factor = adata.obs['size_factor'].values

    train_dataset = Batch_Dataset(raw_mat, exp_mat, aspects, batches, size_factor, condition_types, batch_labels)

    train_loader = DataLoader(train_dataset,
                              batch_size=batch_size,
                              shuffle=True,
                              drop_last=False)
    n_iters = len(train_loader)

    batch_labels = torch.tensor(batch_labels)
    #######################   Prepare models   #############################################################################
    input_dim = adata.X.shape[1]
    cond_cls = {}
    for (aspect, n_cond) in n_conditions.items():
        cond_cls[aspect] = Classifier(z_dim=cond_dim, n_batch=n_cond).to(device) ##原来用的Classifier

    model = scDisco(input_dim=input_dim,
                       bio_z_dim=bio_z_dim,
                       cond_dim=cond_dim,
                       EncLayer=EncLayer,
                       DecLayer=DecLayer,
                       n_batches=n_batches,
                       n_conditions=n_conditions,
                       aspect_key=aspect_key,
                       ).to(device)


    #################################   Prepare optimzers & loss   #################################
    optimizer_m = optim.Adam(model.parameters(), lr=lr_m)
    scheduler_m = StepLR(optimizer_m, step_size=1, gamma=0.99)

    optimizer_cc = {aspect: optim.Adam(cond_cls[aspect].parameters(), lr=lr_c) for aspect in aspect_key}
    scheduler_cc = {aspect: StepLR(optimizer_cc[aspect], step_size=1, gamma=0.99) for aspect in aspect_key}

    MSE = torch.nn.MSELoss().to(device)
    DisB = nn.CrossEntropyLoss().to(device)
    KLD = ELOBkldloss().to(device)

    ###################################   some hyper-parameters   ###################################
    update_interval = 5

    ########################################   Training   ##########################################
    all_ari_k, all_ari_l = [], []
    all_nmi_k, all_nmi_l = [], []
    all_loss_ = []
    best_all_ari_k = 0
    best_all_ari_l = 0
    best_epoch_k = 0
    best_epoch_l = 0

    for epoch in range(1, n_epochs+1):

        iter_mse_loss = 0.
        iter_kld_loss = 0.
        iter_cls_loss = 0.
        iter_all_loss = 0.

        pred_c_all = {}

        for iters, (idx, raw, exp, aspects, batch, sf, c_labels, b_labels) in enumerate(train_loader):
            idx, raw, exp, batch, sf= idx.to(device), raw.to(device), exp.to(device), batch.to(device), sf.to(device)
            aspects = {aspect: cond.to(device) for (aspect, cond) in aspects.items()}
            b_labels = b_labels.to(device)
            c_labels = {aspect: c_label.to(device) for (aspect, c_label) in c_labels.items()}

            # Step 2: Train encoder, decoder and batch-classifier
            z_mu_bio, z_logvar_bio, z_bio, z_mu_cond, z_logvar_cond, z_cond, x_x, z_cond3 = model(exp, batch, c_labels)

            # MSE losses
            mse_loss = MSE(raw, x_x)

            # KLD Losses
            kld_loss_1 = KLD(z_mu_bio, z_logvar_bio)
            kld_loss_2 = 0
            for aspect in aspect_key:
                kld_loss_2 += KLD(z_mu_cond[aspect], z_logvar_cond[aspect])
            kld_loss_2 = kld_loss_2 / n_aspects
            kld_loss = kld_loss_1 + kld_loss_2

            # condition-classifier
            cls_loss = 0.
            for (aspect, cls) in cond_cls.items():
                cond_logit = cls(z_cond3[aspect])
                cls_loss = cls_loss + DisB(cond_logit, aspects[aspect])
                pred_c = torch.max(cond_logit, dim=1)[1]
                pred_c_all.setdefault(aspect, []).extend(pred_c.detach().cpu().numpy())

            # Step 3: Train all loss
            # All losses
            all_loss = mse_loss + lda * kld_loss + lda_c * cls_loss / n_aspects

            optimizer_m.zero_grad()
            for aspect in aspect_key:
                optimizer_cc[aspect].zero_grad()

            all_loss.backward()
            optimizer_m.step()
            for aspect in aspect_key:
                optimizer_cc[aspect].step()


            iter_mse_loss += mse_loss
            iter_kld_loss += kld_loss
            iter_cls_loss += cls_loss
            iter_all_loss += all_loss

        epoch_mse_loss = iter_mse_loss/n_iters
        epoch_kld_loss = iter_kld_loss/n_iters
        epoch_cls_loss = iter_cls_loss/n_iters
        epoch_all_loss = iter_all_loss/n_iters

        print('Epoch [{}/{}],  mse_loss:{:.4f}, kld_loss:{:.4f}, cls_loss:{:.4f}, all_loss:{:.4f}'.format(
               epoch, n_epochs, epoch_mse_loss, epoch_kld_loss, epoch_cls_loss, epoch_all_loss))

        all_loss_0 = all_loss
        all_loss_.append(all_loss_0.detach().cpu().numpy())

        if epoch % 5 == 0:
            scheduler_m.step()
            {aspect: scheduler_cc[aspect].step() for aspect in aspect_key}

    # # # ############################   Return results   ###################################
    record = None
    if cl_type:
         record = {
             'loss': all_loss_
         }

    adata.obsm['bio_feat'], cond_feat = model.EncodeAll(exp_mat, batch_onehot, condition_labels, device=device)
    for (aspect, feat) in cond_feat.items():
        adata.obsm['cond_feat' + '_' + aspect] = feat
    torch.save(model.state_dict(), model_dir)

    return adata, record
