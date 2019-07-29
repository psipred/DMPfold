#!/usr/bin/env python

# Predict distance probabilities with previous seed model

# By David T. Jones, Jan 2019 */

# Copyright (C) 2019 University College London

# License: GPLv3

from __future__ import print_function

import sys
import os
import time

from math import sqrt, log

import numpy as np
from scipy.spatial.distance import pdist, squareform

import torch
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F

from nndef_iterdist_maxoutresnet import ResNet

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

    network.load_state_dict(torch.load(scriptdir + '/FINAL_fullmap_iterdist_model.pt', map_location=lambda storage, loc: storage))

    mapdata = np.fromfile(sys.argv[1], dtype=np.float32)
    length = int(sqrt(mapdata.shape[0]/441))
    inputs = mapdata.reshape(1,441,length,length)

    mapdata = np.fromfile(sys.argv[2], dtype=np.float32)

    coords = []
    with open(sys.argv[3], 'r') as refpdbfile:
        for line in refpdbfile:
            if line[:4] == 'ATOM' and line[12:16] == ' CA ':
                # Split the line
                pdb_fields = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
                coords.append(np.array([float(pdb_fields[6]), float(pdb_fields[7]), float(pdb_fields[8])]))
            if line[:3] == 'END':
                assert length == len(coords), "Number of residues in PDB and 21c files differ - does your FASTA file match the map/21c file?"
                cv = np.asarray(coords)
                initdmat = squareform(0.1*pdist(cv, 'euclidean')).astype(np.float32).reshape(1,1,length,length)

    with torch.no_grad():
        inputs = torch.from_numpy(np.concatenate((initdmat, mapdata.reshape(1,MAP_CHANNELS,length,length), inputs), axis=1)).type(torch.FloatTensor).to(device)

        output = network(inputs)

        result = F.softmax(0.5 * (output + output.transpose(2,3)), dim=1)

        for wi in range(0, length-1):
            for wj in range(wi+5, length):
                print("{} {} ".format(wi+1,wj+1), end='')
                probs = result.data[0,:,wi,wj]
                for bn in range(20):
                    print(" {}".format(probs[bn]), end='')
                print("")

if __name__=="__main__":
    main()
