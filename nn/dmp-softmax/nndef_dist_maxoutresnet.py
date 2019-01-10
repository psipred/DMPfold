import torch
import torch.nn as nn
import torch.nn.functional as F
from math import sqrt

NUM_CHANNELS = 441+60

class Maxout2d(nn.Module):

    def __init__(self, in_channels, out_channels, pool_size, kernel_size=1, dilation=1):
        super(Maxout2d, self).__init__()
        self.in_channels, self.out_channels, self.pool_size = in_channels, out_channels, pool_size
        self.lin = nn.Conv2d(in_channels=in_channels, out_channels=out_channels * pool_size, kernel_size=kernel_size, dilation=dilation, padding=dilation*(kernel_size-1)//2)
        self.norm = nn.InstanceNorm2d(out_channels, affine=True)
        nn.init.xavier_uniform_(self.lin.weight, gain=sqrt(2.0))

    def forward(self, inputs):
        x = self.lin(inputs)

        N,C,H,W = x.size()

        x = x.view(N, C//self.pool_size, self.pool_size, H, W)
        x = x.max(dim=2)[0]
        x = self.norm(x)

        return x

# ResNet Module
class ResNet(nn.Module):
    def __init__(self,width):
        super(ResNet, self).__init__()
        self.redim = Maxout2d(in_channels=NUM_CHANNELS, out_channels=width, pool_size=3)
        self.resblocks = nn.ModuleList()

        for fsize,dilv in [(5,1), (5,2), (5,1), (5,4), (5,1), (5,8), (5,1), (5,16), (5,1), (5,32), (5,1), (5,64), (5,1), (5,128), (5,1), (5,1), (5,1), (5,1)]:
             if fsize > 0:
                layer = Maxout2d(in_channels=width, out_channels=width, pool_size=4, kernel_size=fsize, dilation=dilv)
                self.resblocks.append(layer)

        self.outlayer = nn.Conv2d(in_channels=width, out_channels=20, kernel_size=1)
        nn.init.xavier_uniform_(self.outlayer.weight)

    def forward(self, x):

        x = self.redim(x)

        for block in self.resblocks:
            residual = x
            x = block(x)
            x += residual

        out = self.outlayer(x)

        return out
