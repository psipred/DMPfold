# DMPfold

[![Build Status](https://travis-ci.org/psipred/DMPfold.svg?branch=master)](https://travis-ci.org/psipred/DMPfold)

Extending genome-scale de novo protein modelling coverage using iterative deep learning-based prediction of structural constraints

See the [pre-print](https://arxiv.org/abs/1811.12355) for more.

## Installation

As it makes use of a lot of different software, installation can be a little fiddly.
However we have aimed to make it as straightforward as possible.
These instructions should work for a Linux system:
- Make sure you have Python 3 with [PyTorch](https://pytorch.org) 0.4 or later, NumPy and SciPy installed. GPU setup is optional for Pytorch - it won't speed things up much because running the network isn't a time-consuming step. DMPfold has been tested on Python 3.6 and 3.7.
- Install [HH-suite](https://github.com/soedinglab/hh-suite) and the uniclust30 database, unless you are getting your alignments from elsewhere.
- Install [FreeContact](https://rostlab.org/owiki/index.php/FreeContact).
- Install [CCMpred](https://github.com/soedinglab/CCMpred).
- Install [MODELLER](https://salilab.org/modeller), which requires a license key. Only the Python package is required so this can be installed with `conda install modeller -c salilab`.
- Install [CNS](http://cns-online.org/v1.3). We found we had to follow all of the steps in [this comment](https://ask.bioexcel.eu/t/cns-errors-before-after-recompilation/54/14) to get this working, and had to install the `flex-devel` package via our system package manager. Change the `nrestraints = 20000` line in `cns_solve_1.3/modules/nmr/readdata` to a larger number, e.g. `nrestraints = 30000`, and rebuild CNS to allow DMPfold to run on larger structures. Also set set MXFPEPS2 value to 8192 in machvar.inc. Edit MXRTP to 4000 in rtf.inc in the source directory. 

WRITE (6,’(I6,E10.3,E10.3)’) I, ONEP, ONEM
just above line 67, which looks like this:
IF (ONE .EQ. ONEP .OR. ONE .EQ. ONEM) THEN

- Download and patch the required CNS scripts by changing into the `cnsfiles` directory and running `sh installscripts.sh`.
- Install the [legacy BLAST](https://tinyurl.com/y57hq2wo) software, in particular `formatdb`, `blastpgp` and `makemat`. We may update this to BLAST+ in the future.
- Other software is pre-compiled and included here (PSIPRED, PSICOV, various utility scripts with the code in `src`). This should run okay but may need separate compilation using the makefile if issues arise. Some other standard programs, such as csh shell, are assumed.
- Change lines 10/13-15/18/21/24 in `seq2maps.csh`, lines 11/14/17/20 in `aln2maps.csh`, lines 4/7 in `bin/runpsipredandsolvwithdb` and lines 10/13 in `run_dmpfold.sh` to point to the installed locations of the above software. You can also set the number of cores to use in `seq2maps.csh` and `aln2maps.csh`. This sets the number of cores for HHblits, PSICOV, FreeContact and CCMpred - the script will run faster with this set to a value larger than 1 (e.g. 4 or 8).

Check the continuous integration [setup script](.travis.yml) and [logs](https://travis-ci.org/psipred/DMPfold) for additional tips and a step-by-step installation on Ubuntu.

## Usage

Here we give an example of running DMPfold on Pfam family PF10963.
First you need to generate the `.21c` and `.map` files.
This can be done in one of two ways:
- From a single sequence: `csh seq2maps.csh PF10963.fasta` to run HHblits, PSIPRED, SOLVPRED, PSICOV, FreeContact, CCMpred and alnstats.
- From an alignment: `csh aln2maps.csh PF10963.aln` to run PSIPRED, SOLVPRED, PSICOV, FreeContact, CCMpred and alnstats. The file `PF10963.aln` has one sequence per line with the ungapped target sequence as the first line.

Then run `sh run_dmpfold.sh PF10963.fasta PF10963.21c PF10963.map ./PF10963` to run DMPfold, where the last parameter is an output directory that will be created.
The final model is `final_1.pdb` and other structures may or may not be generated as `final_2.pdb` to `final_5.pdb` if they are significantly different.
Running `sh run_dmpfold.sh PF10963.fasta PF10963.21c PF10963.map ./PF10963 5 20` instead runs 5 iterations with 20 models per iteration (default is 3 and 50).

## Data

Models for the 1,475 [Pfam](http://pfam.xfam.org) families modelled in the paper can be downloaded [here](http://bioinf.cs.ucl.ac.uk/downloads/dmpfold/pfam_models.tgz).
Additional models for the remainder of the dark Pfam families can be downloaded [here](http://bioinf.cs.ucl.ac.uk/downloads/dmpfold/pfam_lowconf_models.tgz) (some were not modelled due to small sequence alignments).
Alignments for the Pfam families without available templates can be downloaded [here](http://bioinf.cs.ucl.ac.uk/downloads/dmpfold/pfam_alignments.tgz).
The format is one sequence per line with the ungapped target sequence as the first line.

The directory [pfam](pfam) in this repository contains text files with the lists from Figure 4A of the paper, target sequences for modelled families and data for modelled families (sequence length, effective sequence count, distogram satisfaction scores, estimated TM-score and probability TM-score >= 0.5).
