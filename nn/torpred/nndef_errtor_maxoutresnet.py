# Neural network definition file

# By David T. Jones, Jan 2019

# Copyright (C) 2019 University College London

# License: GPLv3

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
class ResNet_Block(nn.Module):
    def __init__(self,width,fsize,dilv):
        super(ResNet_Block, self).__init__()
        self.layer1 = Maxout2d(in_channels=width, out_channels=width, pool_size=4, kernel_size=fsize, dilation=dilv)

    def forward(self, x):

        residual = x
        out = self.layer1(x)
        out += residual

        return out

# ResNet Module
class ResNet(nn.Module):
    def __init__(self,width):
        super(ResNet, self).__init__()
        self.redim = Maxout2d(in_channels=NUM_CHANNELS, out_channels=width, pool_size=3)

        layers = []

        for fsize,dilv in [(5,1), (5,2), (5,1), (5,4), (5,1), (5,8), (5,1), (5,16), (5,1), (5,32), (5,1), (5,64), (5,1), (5,128), (5,1), (5,1), (5,1), (5,1)]:
             if fsize > 0:
                layer = ResNet_Block(width, fsize, dilv)
                layers.append(layer)

        self.lstm = nn.LSTM(width, width * 2, num_layers=1, bidirectional=True)
        self.outlayer = nn.Conv1d(in_channels=width*4, out_channels=3, kernel_size=1)
        self.outlayer.weight.data.zero_()
        self.outlayer.bias.data.zero_()

        self.resnet_model = nn.Sequential(*layers)

    def forward(self, x):

        x = self.redim(x)

        x = self.resnet_model(x)

        nrowcols = x.size()[2]
        x = x[0,:,:,:].permute(2,1,0)
        lstm_out, hidden = self.lstm(x)
        x = lstm_out[[-1],:,:].permute(0,2,1)
        x = self.outlayer(x)

        #out = torch.log(1 + torch.exp(50 * x) / (1 + torch.exp(50 * (x - 1)))) / 50
        out = x

        return out
