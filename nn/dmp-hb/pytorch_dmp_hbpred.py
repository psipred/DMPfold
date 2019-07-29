#!/usr/bin/env python

# Predict main chain H-bond donors & acceptors

# By David T. Jones, Jan 2019

# Copyright (C) 2019 University College London

# License: GPLv3

from __future__ import print_function

import sys
import os
import time

from math import sqrt, log

import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

from nndef_dmp_resnet import ResNet

MAP_CHANNELS = 60

# ############################## Main program ################################
# Everything else will be handled in our main program now. We could pull out
# more functions to better separate the code, but it wouldn't make it any
# easier to read.

def main():

    device = torch.device("cpu")
     
     # Create neural network model (depending on first command line parameter)
    network = ResNet(64).eval().to(device)

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    network.load_state_dict(torch.load(scriptdir + '/FINAL_fullmap_hb_model.pt', map_location=lambda storage, loc: storage))

    mapdata = np.fromfile(sys.argv[1], dtype=np.float32)
    length = int(sqrt(mapdata.shape[0]/441))
    inputs = mapdata.reshape(1,441,length,length)

    mapdata = np.fromfile(sys.argv[2], dtype=np.float32)

    with torch.no_grad():
        inputs = torch.from_numpy(np.concatenate((mapdata.reshape(1,MAP_CHANNELS,length,length), inputs), axis=1)).type(torch.FloatTensor).to(device)
        result = torch.sigmoid(network(inputs).data)

    for wi in range(length):
#        print("{} {} 3.8 3.8".format(wi+1,wi+2))
        for wj in range(length):
            if abs(wi - wj) > 2:
                print("{} {} 0 3.5 ".format(wi+1,wj+1), end='')
                print(" {}".format(result[0,0,wi,wj]))

if __name__=="__main__":
    main()
