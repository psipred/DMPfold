#!/usr/bin/env python

# Secondary structure prediction with long range information using Deep Residual Network
# by David T. Jones 2018

from __future__ import print_function

import sys
import os

from math import sqrt,atan2

import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

from nndef_sincostor_conv import ConvNet
from nndef_errtor_conv import ConvNet as ErrConvNet


# ############################## Main program ################################
# Everything else will be handled in our main program now. We could pull out
# more functions to better separate the code, but it wouldn't make it any
# easier to read.

def main():
    # Create neural network model
    network = ConvNet(64).cpu().eval()
    network2 = ErrConvNet(64).cpu().eval()

    scriptdir = os.path.dirname(os.path.realpath(__file__))

    network.load_state_dict(torch.load(scriptdir + '/FINAL_torcovpred_model.pt', map_location=lambda storage, loc: storage))
    network2.load_state_dict(torch.load(scriptdir + '/FINAL_torcoverrpred_model.pt', map_location=lambda storage, loc: storage))

    mapdata = np.fromfile(sys.argv[1], dtype=np.float32)

    length = int(sqrt(mapdata.shape[0]/60))

    with torch.no_grad():
        inputs = torch.from_numpy(mapdata.reshape(1,60,length,length)).type(torch.FloatTensor)

        result = network(inputs).data
        prederr = torch.clamp(network2(inputs).data * 180.0, 1.0, 180.0) 

    print("! TORPRED CNSsolve FORMAT (TORCOVPRED V0.1)")

    for wi in range(0, length):
        sin_phi, cos_phi, sin_psi, cos_psi, sin_omega, cos_omega = result[0, :, wi]
        phi = atan2(sin_phi, cos_phi) * 180.0 / np.pi
        psi = atan2(sin_psi, cos_psi) * 180.0 / np.pi
        omega = atan2(sin_omega, cos_omega) * 180.0 / np.pi
        err_phi, err_psi, err_omega = prederr[0, :, wi]
        #print("%4d %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f" % (wi+1, phi, err_phi, psi, err_psi, omega, err_omega))
        if wi > 0:
            print("assign (resid {} and name c) (resid {} and name n) (resid {} and name ca) (resid {} and name c) 1.0 {:.2f} {:.2f} 2  ! phi".format(wi, wi+1, wi+1, wi+1, phi, err_phi))
        if wi < length-1:
            print("assign (resid {} and name n) (resid {} and name ca) (resid {} and name c) (resid {} and name n) 1.0 {:.2f} {:.2f} 2  ! psi".format(wi+1, wi+1, wi+1, wi+2, psi, err_psi))

if __name__=="__main__":
    main()
