import math
import torch
import torch.nn as nn
from model.network import buildNetwork
from util.data_utils import *

class scDisco(nn.Module):
    def __init__(self, input_dim,
                 bio_z_dim,
                 cond_dim,
                 EncLayer,
                 DecLayer,
                 n_batches,
                 n_conditions,
                 aspect_key,
                 activation="relu"):
        super(scDisco, self).__init__()

        device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")

        self.bio_z_dim = bio_z_dim
        self.cond_zim = cond_dim
        self.n_batches = n_batches
        self.aspect_key = aspect_key
        self.enc_dim = input_dim + n_batches
        self.dec_dim = bio_z_dim + n_batches + cond_dim * len(aspect_key)
        self.n_conditions = n_conditions
        self.activation = activation

        ########################  Encoder & Decoder #########################
        self.Encoder = buildNetwork([self.enc_dim] + EncLayer, activation=activation).to(device)
        self.Decoder = buildNetwork([self.dec_dim] + DecLayer, activation=activation).to(device)

        self.enc_bio = nn.Linear(EncLayer[-1], bio_z_dim).to(device)
        self.enc_batch = nn.Linear(bio_z_dim, n_batches).to(device)

        self.bn_list = nn.ModuleList()
        for i in range(n_batches):
            self.bn_list.append(
                nn.BatchNorm1d(bio_z_dim, eps=1e-5, momentum=0.1, affine=True, track_running_stats=True)).to(device)

        self.cn_list = {}
        for aspect in self.aspect_key:
            self.cn_list[aspect] = nn.ModuleList()
            for i in range(n_conditions[aspect]):
                self.cn_list[aspect].append(
                    nn.BatchNorm1d(cond_dim, eps=1e-5, momentum=0.1, affine=True, track_running_stats=True)).to(device)

        self.enc_cond = {}
        for aspect in self.aspect_key:
            self.enc_cond[aspect] = nn.Linear(EncLayer[-1], cond_dim).to(device)

        self.dec_x = nn.Sequential(nn.Linear(DecLayer[-1], input_dim), nn.ReLU()).to(device)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def Encoder_all(self, x, b, c_labels):
        x_ = torch.cat([x, b], 1)
        h = self.Encoder(x_)

        z_bio = self.enc_bio(h)

        z_cond = {}
        for (aspect, enc) in self.enc_cond.items():
            z_cond[aspect] = enc(h)

        z_cond2 = {}
        for aspect in self.aspect_key:
            z_cond2[aspect] = z_cond[aspect]

        z2_ = {}
        for aspect in self.aspect_key:
            z2_[aspect] = []
            for i in range(self.n_conditions[aspect]):
                locals()['index_' + str(i) + aspect] = torch.argwhere(c_labels[aspect] == i).squeeze()
                locals()['z2_' + str(i) + aspect] = z_cond[aspect][locals()['index_' + str(i) + aspect]]
                if locals()['z2_' + str(i) + aspect].dim() == 2:
                    z2_[aspect].append(self.cn_list[aspect][i](locals()['z2_' + str(i) + aspect]))
                else:
                    z2_[aspect].append((locals()['z2_' + str(i) + aspect]).unsqueeze(-1).T)

            z_cond[aspect] = torch.cat(z2_[aspect], 0)

        z_cond3 = {}
        for aspect in self.aspect_key:
            z_cond3[aspect] = z_cond2[aspect] - z_cond[aspect]

        return z_bio, z_cond, z_cond2, z_cond3

    def forward(self, x, b, b_labels):
        z_bio, z_cond, z_cond2, z_cond3 = self.Encoder_all(x, b, b_labels)

        z_cond_ = []
        for aspect in z_cond2:
            z_cond_.append(z_cond2[aspect])
        z_cond_all = torch.cat(z_cond_, dim=1)

        z = torch.cat([z_bio, b, z_cond_all], dim=1)

        # Decode
        h_d = self.Decoder(z)
        x_x = self.dec_x(h_d)

        return z_bio, z_cond, x_x, z_cond3

    def EncodeAll(self, X, B, c_labels, device, batch_size=256):
        all_z_bio, all_z_cond = [], {}
        all_z_cond_f = {}
        c_l = {}
        for aspect in self.aspect_key:
            all_z_cond[aspect] = []
            all_z_cond_f[aspect] = []
            c_l[aspect] = []

        num_cells = X.shape[0]
        num_batch = int(math.ceil(1.0 * X.shape[0] / batch_size))
        for batch_idx in range(num_batch):
            exp = X[batch_idx * batch_size: min((batch_idx + 1) * batch_size, num_cells)]
            exp = torch.from_numpy(np.float32(exp)).to(device)

            b = B[batch_idx * batch_size: min((batch_idx + 1) * batch_size, num_cells)]
            b = torch.tensor(np.float32(b)).to(device)
            for aspect in self.aspect_key:
                c_l[aspect] = c_labels[aspect][batch_idx * batch_size: min((batch_idx + 1) * batch_size, num_cells)]
            with torch.no_grad():
                z_bio, _, z_cond_f, z_cond = self.Encoder_all(exp, b, c_l)

            all_z_bio.append(z_bio.data)
            for (aspect, cond) in z_cond_f.items():
                all_z_cond_f[aspect].append(cond.data)
            for (aspect, cond) in z_cond.items():
                all_z_cond[aspect].append(cond.data)

        all_z_bio = torch.cat(all_z_bio, dim=0)
        all_z_cond_f_ = {}
        for (aspect, cond) in all_z_cond_f.items():
            all_z_cond_f[aspect] = torch.cat(cond, dim=0)
            all_z_cond_f_[aspect] = all_z_cond_f[aspect].cpu().numpy()
        all_z_cond_ = {}
        for (aspect, cond) in all_z_cond.items():
            all_z_cond[aspect] = torch.cat(cond, dim=0)
            all_z_cond_[aspect] = all_z_cond[aspect].cpu().numpy()

        return all_z_bio.cpu().numpy(), all_z_cond_f_, all_z_cond_
