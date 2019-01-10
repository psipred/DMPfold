import torch
import torch.nn as nn
import torch.nn.functional as F
from math import sqrt

NUM_CHANNELS = 60

class Maxout2d(nn.Module):

    def __init__(self, in_channels, out_channels, pool_size):
        super(Maxout2d, self).__init__()
        self.in_channels, self.out_channels, self.pool_size = in_channels, out_channels, pool_size
        self.lin = nn.Conv2d(in_channels=in_channels, out_channels=out_channels * pool_size, kernel_size=1, bias=False)
        nn.init.xavier_uniform_(self.lin.weight, gain=sqrt(2.0))
        self.norm = nn.BatchNorm2d(out_channels * pool_size, affine=True)

    def forward(self, inputs):
        shape = list(inputs.size())
        out = self.norm(self.lin(inputs))
        m = out.view(shape[0], self.out_channels, self.pool_size, shape[2], shape[3]).max(dim=2)[0]
        return m

class Maxout1d(nn.Module):

    def __init__(self, in_channels, out_channels, pool_size):
        super(Maxout1d, self).__init__()
        self.in_channels, self.out_channels, self.pool_size = in_channels, out_channels, pool_size
        self.lin = nn.Conv1d(in_channels=in_channels, out_channels=out_channels * pool_size, kernel_size=1, bias=False)
        nn.init.xavier_uniform_(self.lin.weight, gain=sqrt(2.0))
        self.norm = nn.BatchNorm1d(out_channels * pool_size, affine=True)

    def forward(self, inputs):
        shape = list(inputs.size())
        out = self.norm(self.lin(inputs))
        m = out.view(shape[0], self.out_channels, self.pool_size, shape[2]).max(dim=2)[0]
        return m

# ConvNet Module
class ConvNet(nn.Module):
    def __init__(self,width):
        super(ConvNet, self).__init__()
        self.redim = Maxout2d(in_channels=NUM_CHANNELS, out_channels=width, pool_size=2)
        self.resblocks = nn.ModuleList()

        for i in range(8):
            layer = nn.Conv2d(in_channels=width, out_channels=width, kernel_size=5, dilation=1, padding=2, bias=False)
            nn.init.xavier_uniform_(layer.weight, gain=sqrt(2.0))
            self.resblocks.append(layer)
            self.resblocks.append(nn.BatchNorm2d(width, affine=False))

        self.outlayer1 = Maxout1d(width*2, width, 3)
        self.outlayer2 = nn.Conv1d(in_channels=width, out_channels=6, kernel_size=1)
#        self.lastnorm = nn.BatchNorm1d(width, affine=True)
        nn.init.xavier_uniform_(self.outlayer2.weight)

    def forward(self, x):
        out = self.redim(x)
        preconvnet = out.sum(dim=3)

        for i in range(len(self.resblocks)//2):
            out = self.resblocks[i*2](out)
            out = F.relu(self.resblocks[i*2+1](out))
        
        out = out.sum(dim=3)
        out = torch.cat((preconvnet, out), dim=1)
        out = F.dropout(out, p=0.5)
        out = self.outlayer1(out)
        out = torch.tanh(self.outlayer2(out))
        return out
