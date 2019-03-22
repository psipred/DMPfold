# DMPfold

## NB this repo is under development and is not ready for production use

[![Build Status](https://travis-ci.org/psipred/DMPfold.svg?branch=master)](https://travis-ci.org/psipred/DMPfold)

Fast de novo protein model generation from covarying sequences using predicted distances and iterative model building.

See the [pre-print](https://arxiv.org/abs/1811.12355) for more.

## Installation

As it makes use of a lot of different software, installation can be a little fiddly.
However we have aimed to make it as straightforward as possible.
These instructions should work for a Linux system:
- Make sure you have Python 3 with [PyTorch](https://pytorch.org) 0.4 or later, NumPy and SciPy installed. GPU setup is optional for Pytorch - it won't speed things up much because running the network isn't a time-consuming step. DMPfold has been tested on Python 3.6 and 3.7.
- Install [HH-suite](https://github.com/soedinglab/hh-suite) and the uniclust30 database, unless you are getting your alignments from elsewhere.
- Install [FreeContact](https://rostlab.org/owiki/index.php/FreeContact).
- Install [CCMpred](https://github.com/soedinglab/CCMpred).
- Install [CNS](http://cns-online.org/v1.3). Change the `nrestraints = 20000` line in `cns_solve_1.3/modules/nmr/readdata` to a larger number, e.g. `nrestraints = 30000`, to allow DMPfold to run on larger structures.
- Install [MODELLER](https://salilab.org/modeller), which requires a license key. Only the Python package is required so this can be installed with `conda install modeller -c salilab`.
- Other software is pre-compiled and included here (PSIPRED, PSICOV, various utility scripts with the code in `src`). This should run okay but may need separate compilation using the makefile if issues arise. Some other standard programs, such as csh shell, are assumed.
- You also need access to the legacy BLAST software, in particular `blastpgp` and `makemat`. We will look to update this to BLAST+ soon.
- Change lines 10/13/14/15/18/21 in `seq2maps.csh`, lines 4/7 in `bin/runpsipredandsolvwithdb` and lines 10/13 in `run_dmpfold.sh` to point to the installed locations of the above software. You can also set the number of cores to use in `seq2maps.csh`.

You may need to set `ulimit -s unlimited` to get `seq2maps.csh` to work.
Check the continuous integration [setup script](.travis.yml) and [logs](https://travis-ci.org/psipred/DMPfold) for additional tips and a step-by-step installation on Ubuntu.

## Usage

Example of running on CASP12 target T0864.
First you need to generate the `.21c` and `.map` files.
This can be done in one of two ways:
- From a single sequence: `csh seq2maps.csh T0864.fasta` to run HHblits, PSIPRED, SOLVPRED, PSICOV, FreeContact, CCMpred and alnstats.
- From an alignment: `csh aln2maps.csh T0864.aln` to run PSIPRED, SOLVPRED, PSICOV, FreeContact, CCMpred and alnstats. The file `T0864.aln` has one sequence per line with the ungapped target sequence as the first line.

Then run `sh run_dmpfold.sh T0864.fasta T0864.21c T0864.map ./T0864` to run DMPfold, where the last parameter is an output directory that will be created.
The final model is `final_1.pdb` and other structures may or may not be generated as `final_2.pdb` to `final_5.pdb` if they are significantly different.
Running `sh run_dmpfold.sh T0864.fasta T0864.21c T0864.map ./T0864 5 20` instead runs 5 iterations with 20 models per iteration (default is 3 and 50).

## Data

Models for the 1,475 [Pfam](http://pfam.xfam.org) families modelled in the paper can be downloaded [here](http://bioinf.cs.ucl.ac.uk/downloads/dmpfold/pfam_models.tgz).

The directory [pfam](pfam) in this repository contains text files with the lists from Figure 4A of the paper and target sequences for modelled families.
