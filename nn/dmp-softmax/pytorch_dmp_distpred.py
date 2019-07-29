#!/usr/bin/env python

# Predict first stage distance probabilities

# By David T. Jones, Jan 2019 */

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

from nndef_dist_maxoutresnet import ResNet

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

    network.load_state_dict(torch.load(scriptdir + '/FINAL_fullmap_dist_model.pt', map_location=lambda storage, loc: storage))

    mapdata = np.fromfile(sys.argv[1], dtype=np.float32)
    length = int(sqrt(mapdata.shape[0]/441))
    inputs = mapdata.reshape(1,441,length,length)

    mapdata = np.fromfile(sys.argv[2], dtype=np.float32)

    inputs = torch.from_numpy(np.concatenate((mapdata.reshape(1,MAP_CHANNELS,length,length), inputs), axis=1)).type(torch.FloatTensor)
    
    network.eval()
    with torch.no_grad():
        output = network(inputs)
        result = F.softmax(0.5 * (output + output.transpose(2,3)), dim=1)

    for wi in range(0, length-1):
        for wj in range(wi+5, length):
            probs = result.data[0,:,wi,wj]
            print("{} {}".format(wi+1,wj+1), end='')
            for dbin in range(20):
                print(" {}".format(probs[dbin]), end='')
            print('')

if __name__=="__main__":
    main()
