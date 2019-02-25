# DMPfold

## NB this repo is under development and is not ready for production use

[![Build Status](https://travis-ci.org/psipred/DMPfold.svg?branch=master)](https://travis-ci.org/psipred/DMPfold)

Fast de novo protein model generation from covarying sequences using predicted distances and iterative model building.

See the [pre-print](https://arxiv.org/abs/1811.12355) for more.

## Installation

As it makes use of a lot of different software, installation can be a little fiddly.
However we have aimed to make it as straightforward as possible.
These instructions should work for a Linux system:
- Make sure you have Python 3 with [PyTorch](https://pytorch.org) 0.4 or later, NumPy and SciPy installed. GPU setup is optional for Pytorch. DMPfold has been tested on Python 3.6 and 3.7.
- Install [HH-suite](https://github.com/soedinglab/hh-suite) and the uniclust30 database, unless you are getting your alignments from elsewhere.
- Install [FreeContact](https://rostlab.org/owiki/index.php/FreeContact).
- Install [CCMPred](https://github.com/soedinglab/CCMpred).
- Install [CNS](http://cns-online.org/v1.3).
- Install [MODELLER](https://salilab.org/modeller), which requires a license key. Only the Python package is required so this can be installed with `conda install modeller -c salilab`.
- Other software is pre-compiled and included here (PSIPRED, PSICOV, various utility scripts with the code in `src`). This should run okay but may need separate compilation if issues arise.
- You also need access to the legacy BLAST software. We will look to update this to BLAST+ soon.
- Change lines 10/13/14/15/18/21 in `seq2maps.csh`, lines 10/13 in `run_dmpfold.sh` and lines 4/7 in `bin/runpsipredandsolvwithdb` to point to the installed locations of the above software. You can also set the number of cores to use in `seq2maps.csh`.

You may need to set `ulimit -s unlimited` to get `seq2maps.csh` to work.
Check the continuous integration [logs](https://travis-ci.org/psipred/DMPfold) for additional tips and a step-by-step install on Ubuntu.

## Usage

Example of running on CASP12 target T0864:

- `csh seq2maps.csh T0864.fasta` to run HHblits, PSIPRED, SOLVPRED, PSICOV, FreeContact, CCMpred and alnstats. Outputs `.21c` and `.map` files.
- `sh run_dmpfold.sh T0864.fasta T0864.21c T0864.map ./T0864` to run DMPfold, where the last parameter is an output directory that will be created. The final model is `final_1.pdb` and other structures may be generated as `final_2.pdb` to `final_5.pdb`.

Running `sh run_dmpfold.sh T0864.fasta T0864.21c T0864.map ./T0864 5 20` instead runs 5 iterations with 20 models per iteration (default is 3 and 50).
