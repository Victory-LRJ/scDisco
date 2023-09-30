from torch import nn
import torch

device = torch.device("cuda:1" if torch.cuda.is_available() else "cpu")

class ZINBLoss(nn.Module):
    def __init__(self, ridge_lambda=0.0):
        self.ridge_lambda = ridge_lambda

        super(ZINBLoss, self).__init__()

    def forward(self, x, mu, disp, pi, size_factor):
        eps = 1e-10
        mu = (mu.T * size_factor).T
        # torch.lgamma是gamma变换
        t1 = torch.lgamma(disp+eps) + torch.lgamma(x+1.0) - torch.lgamma(x+disp+eps)
        t2 = (disp+x) * torch.log(1.0 + (mu/(disp+eps))) + (x * (torch.log(disp+eps) - torch.log(mu+eps)))
        nb_final = t1 + t2

        nb_case = nb_final - torch.log(1.0-pi+eps)
        zero_nb = torch.pow(disp/(disp+mu+eps), disp)
        zero_case = -torch.log(pi + ((1.0-pi)*zero_nb)+eps)
        result = torch.where(torch.le(x, 1e-8), zero_case, nb_case)

        if self.ridge_lambda > 0:
            ridge = self.ridge_lambda*torch.square(pi)
            result += ridge

        return result.mean()


class ELOBkldloss(nn.Module):
    def __init__(self):
        super(ELOBkldloss, self).__init__()

    def forward(self, mu, logvar):
        result = -((0.5 * logvar) - (torch.exp(logvar) + mu ** 2) / 2. + 0.5)
        result = torch.sum(result, dim=-1)

        return result.mean()

