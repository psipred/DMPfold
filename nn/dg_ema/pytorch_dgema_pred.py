#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import time
import random

import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

# ############################## Main program ################################

def main():

    # Create neural network model
    network = nn.Sequential(
            nn.Linear(4, 10, bias=True),
            nn.SELU(),
            nn.Linear(10, 10, bias=True),
            nn.SELU(),
            nn.Linear(10, 10, bias=True),
        )

    # Load current model snapshot
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    network.load_state_dict(torch.load(scriptdir + '/FINAL_dgema_model.pt', map_location=lambda storage, loc: storage))

    network.eval()
    with torch.no_grad():
        values = ((float(sys.argv[1])-192.2)/139.3, (float(sys.argv[2])-1276.6)/3714.8, (float(sys.argv[3])+0.1072)/0.0419, (float(sys.argv[4])+151.0)/155.4)
        inputs = torch.tensor([values])
        output = F.softmax(network(inputs)/1.2515, dim=1)

        #        print(output)
        #        print(torch.sum(output[0,5:]).item())
        tmscore = 0.0
        for c in range(10):
            tmscore += output[0,c].item() * (0.05 + c * 0.1)
            print(output[0,c].item(), end=' ')
        print(tmscore, torch.sum(output[0,5:]).item())

if __name__=="__main__":
    main()
