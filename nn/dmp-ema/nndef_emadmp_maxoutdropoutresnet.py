import torch.nn as nn
import torch.nn.functional as F
from math import sqrt

NUM_CHANNELS = 441+60

class Maxout2d(nn.Module):

    def __init__(self, in_channels, out_channels, pool_size, kernel_size=1, dilation=1):
        super(Maxout2d, self).__init__()
        self.in_channels, self.out_channels, self.pool_size = in_channels, out_channels, pool_size
        self.lin = nn.Conv2d(in_channels=in_channels, out_channels=out_channels * pool_size, kernel_size=kernel_size, dilation=dilation, padding=dilation*(kernel_size-1)//2)
        nn.init.xavier_uniform_(self.lin.weight, gain=sqrt(2.0))

    def forward(self, inputs):
        shape = list(inputs.size())
        out = self.lin(inputs)
        m = out.view(shape[0], self.out_channels, self.pool_size, shape[2], shape[3]).max(dim=2)[0]
        return m

# ResNet Module
class ResNet(nn.Module):
    def __init__(self,width):
        super(ResNet, self).__init__()
        self.redim = Maxout2d(in_channels=NUM_CHANNELS, out_channels=width, pool_size=3)
        self.firstnorm = nn.InstanceNorm2d(width, affine=True)
        self.resblocks = nn.ModuleList()

        for fsize,dilv in [(5,1), (5,2), (5,1), (5,4), (5,1), (5,8), (5,1), (5,16), (5,1), (5,32), (5,1), (5,64), (5,1), (5,128), (5,1), (5,1), (5,1), (5,1)]:
             if fsize > 0:
                layer = Maxout2d(in_channels=width, out_channels=width, pool_size=2, kernel_size=fsize, dilation=dilv)
                self.resblocks.append(layer)
                self.resblocks.append(nn.InstanceNorm2d(width, affine=True))
                self.resblocks.append(nn.Dropout2d(p=0.5))

        self.pooling = nn.AdaptiveAvgPool2d((4,4))
        self.outlayer = nn.Linear(width*4*4, 6)
        nn.init.xavier_uniform_(self.outlayer.weight)

    def forward(self, x):
        out = self.redim(x)
        out = self.firstnorm(out)
        for i in range(len(self.resblocks)//3):
            residual = out
            out = self.resblocks[i*3](out)
            out = self.resblocks[i*3+1](out)
            out = self.resblocks[i*3+2](out)
            out += residual

        out = self.pooling(out)
        out = out.view(out.size(0), -1)
        out = self.outlayer(out)

        return out
