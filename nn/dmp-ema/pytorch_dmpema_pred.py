#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import time

from math import sqrt, log

import numpy as np

import torch
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F

from nndef_emadmp_maxoutdropoutresnet import ResNet

MAP_CHANNELS = 60

# ############################## Main program ################################
# Everything else will be handled in our main program now. We could pull out
# more functions to better separate the code, but it wouldn't make it any
# easier to read.

def main():

    # Create neural network model (depending on first command line parameter)
    network = ResNet(32).cpu().eval()

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    network.load_state_dict(torch.load(scriptdir + '/FINAL_fullmap_ema_model.pt', map_location=lambda storage, loc: storage))

    mapdata = np.fromfile(sys.argv[1], dtype=np.float32)
    length = int(sqrt(mapdata.shape[0]/441))
    inputs = mapdata.reshape(1,441,length,length)

    mapdata = np.fromfile(sys.argv[2], dtype=np.float32)

    inputs = torch.from_numpy(np.concatenate((mapdata.reshape(1,MAP_CHANNELS,length,length), inputs), axis=1)).type(torch.FloatTensor)
    
    with torch.no_grad():
        output = network(inputs)
        result = F.softmax(output, dim=1)

    for p in result.data[0,:]:
        print(p.item(), end='\t')

    print()

    print(torch.sum(result.data[0,3:]).item())

if __name__=="__main__":
    main()
