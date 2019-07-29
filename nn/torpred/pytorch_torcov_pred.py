#!/usr/bin/env python

# Predict main chain torsion angles from DMP input maps

# By David T. Jones, Jan 2019

# Copyright (C) 2019 University College London

# License: GPLv3


from __future__ import print_function

import sys
import os
import time
import random

from math import sqrt,atan2

import numpy as np

import torch
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from nndef_sincostor_maxoutresnet import ResNet
from nndef_errtor_maxoutresnet import ResNet as ErrNet

# ############################## Main program ################################

def main():#
    device = torch.device("cpu")
    network = ResNet(64).eval().to(device)
    network2 = ErrNet(64).eval().to(device)

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    network.load_state_dict(torch.load(scriptdir + '/FINAL_torcovpred_model.pt', map_location=lambda storage, loc: storage))
    network2.load_state_dict(torch.load(scriptdir + '/FINAL_torcoverrpred_model.pt', map_location=lambda storage, loc: storage))

    mapdata = np.fromfile(sys.argv[1], dtype=np.float32)
    length = int(sqrt(mapdata.shape[0]/441))
    inputs = mapdata.reshape(1,441,length,length)

    mapdata = np.fromfile(sys.argv[2], dtype=np.float32)

    inputs = torch.from_numpy(np.concatenate((mapdata.reshape(1,60,length,length), inputs), axis=1)).type(torch.FloatTensor)
    
    network.eval()
    with torch.no_grad():
        result = network(inputs)
        prederr = torch.clamp(network2(inputs).data * 180.0, 1.0, 180.0) 

        print("! DMPtorsions VFORMAT (DMPtorsions V0.2)")

        for wi in range(0, length):
            sin_phi, sin_psi, sin_omega, cos_phi, cos_psi, cos_omega = result[0, :, wi]
            phi = atan2(sin_phi, cos_phi) * 180.0 / np.pi
            psi = atan2(sin_psi, cos_psi) * 180.0 / np.pi
            omega = atan2(sin_omega, cos_omega) * 180.0 / np.pi
            err_phi, err_psi, err_omega = prederr[0, :, wi]
            #print("%4d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f" % (wi+1, phi, err_phi, psi, err_psi, omega, err_omega))
            if wi > 0:
                print("assign (resid {} and name c) (resid {} and name n) (resid {} and name ca) (resid {} and name c) 1.0 {:.2f} {:.2f} 2  ! phi".format(wi, wi+1, wi+1, wi+1, phi, err_phi))
            if wi < length-1:
                print("assign (resid {} and name n) (resid {} and name ca) (resid {} and name c) (resid {} and name n) 1.0 {:.2f} {:.2f} 2  ! psi".format(wi+1, wi+1, wi+1, wi+2, psi, err_psi))
                # To our knowledge the below line is never actually called
                if abs(omega) < 10.0:
                    print("assign (resid {} and name ca) (resid {} and name c) (resid {} and name n) (resid {} and name ca) 1.0 {:.2f} {:.2f} 2  ! omega".format(wi+1, wi+1, wi+2, wi+2, omega, err_omega))
                    
if __name__=="__main__":
    main()
