import torch
import torch.nn as nn
import torch.nn.functional as F

class Mish(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x):
        return x * torch.tanh(F.softplus(x))

def buildNetwork(layers, activation="relu"):
    net = []
    for i in range(1, len(layers)):
        net.append(nn.Linear(layers[i-1], layers[i]))
        if activation == "relu":
            net.append(nn.ReLU())
        elif activation == "sigmoid":
            net.append(nn.Sigmoid())
        elif activation == "mish":
            net.append(Mish())
        elif activation == "tanh":
            net.append(nn.Tanh())
    return nn.Sequential(*net)

# Different Layer Initialization Parameters
def weights_init_normal(m):
    classname = m.__class__.__name__
    if classname.find("Linear") != -1:
        torch.nn.init.normal_(m.weight.data, 0.0, 0.02)
    elif classname.find("BatchNorm") != -1:
        torch.nn.init.normal_(m.weight.data, 1.0, 0.02)
        torch.nn.init.constant_(m.bias.data, 0.0)

# Classifier
class Classifier(nn.Module):
    def __init__(self, z_dim, n_batch):
        super(Classifier, self).__init__()
        self.z_dim = z_dim
        self.n_batch = n_batch
        self.model = nn.Sequential(
            nn.Linear(z_dim, 100),
            nn.ReLU(),
            nn.Linear(100, 100),
            nn.ReLU(),
        )

        # Classify layers
        self.cls_layer = nn.Sequential(
            nn.Linear(100, n_batch)
        )

    def forward(self, data):
        out = self.model(data)
        classify = self.cls_layer(out)
        return classify

class MeanAct(nn.Module):
    def __init__(self):
        super(MeanAct, self).__init__()

    def forward(self, x):
        return torch.clamp(torch.exp(x), min=1e-5, max=1e6)


class DispAct(nn.Module):
    def __init__(self):
        super(DispAct, self).__init__()

    def forward(self, x):
        return torch.clamp(F.softplus(x), min=1e-4, max=1e4)

