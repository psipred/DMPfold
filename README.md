# DMPfold

## NB this repo is under development and is not ready for production use

[![Build Status](https://travis-ci.org/psipred/DMPfold.svg?branch=master)](https://travis-ci.org/psipred/DMPfold)

Fast de novo protein model generation from covarying sequences using predicted distances and iterative model building.

See the [pre-print](https://arxiv.org/abs/1811.12355) for more.

## Installation

As it makes use of a lot of different software, installation can be a little fiddly.
However we have aimed to make it as straightforward as possible.
These instructions should work for a Linux system:
- Make sure you have Python 3 with [PyTorch](https://pytorch.org) 0.4 or later, NumPy and SciPy installed. DMPfold has been tested on Python 3.6 and 3.7.
- Install [HH-suite](https://github.com/soedinglab/hh-suite) and the uniclust30 database, unless you are getting your alignments from elsewhere.
- Install [FreeContact](https://rostlab.org/owiki/index.php/FreeContact).
- Install [CCMPred](https://github.com/soedinglab/CCMpred).
- Install [CNS](http://cns-online.org/v1.3).
- Install [MODELLER](https://salilab.org/modeller), which requires a license key. Only the Python package is required so this can be installed with `conda`.
- Other software is pre-compiled and included here (PSIPRED, PSICOV, various utility scripts with the code in `src`). This should run okay but may need separate compilation if issues arise.
- Change lines 10/13/14/15/18/21 in `seq2maps.csh` and lines 10/13 in `run_dmpfold.sh` to point to the installed locations of the above software.
Check the continuous integration [logs](https://travis-ci.org/psipred/DMPfold) for additional tips and a step-by-step install on Ubuntu.

## Usage

Example of running on CASP12 target T0864:

- `csh seq2maps.csh T0864.fasta` to run HHblits, PSIPRED, SOLVPRED, PSICOV, FreeContact, CCMpred and alnstats. Outputs 21c and map files.
- `sh run_dmpfold.sh T0864.fasta T0864.21c T0864.map ./T0864` to run DMPfold, where the last parameter is an output directory that will be created.

Running `sh run_dmpfold.sh T0864.fasta T0864.21c T0864.map ./T0864 5 20` instead runs 5 iterations with 20 models per iteration (default is 3 and 50).
