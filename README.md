# DMPfold

[![Build Status](https://travis-ci.org/psipred/DMPfold.svg?branch=master)](https://travis-ci.org/psipred/DMPfold)

Deep learning extends de novo protein modelling coverage of genomes using iteratively predicted structural constraints.

See our [paper in Nature Communications](https://www.nature.com/articles/s41467-019-11994-0) for more.
Please cite the paper if you use DMPfold.

## Installation

As it makes use of a lot of different software, installation can be a little fiddly.
However we have aimed to make it as straightforward as possible.
These instructions should work for a Linux system:
- Make sure you have Python 3 with [PyTorch](https://pytorch.org) 0.4 or later, NumPy and SciPy installed. GPU setup is optional for Pytorch - it won't speed things up much because running the network isn't a time-consuming step. DMPfold has been tested on Python 3.6 and 3.7. The command `python3` should point to the Python that you want to use.
- Install [HH-suite](https://github.com/soedinglab/hh-suite) and the uniclust30 database, unless you are getting your alignments from elsewhere.
- Install [FreeContact](https://rostlab.org/owiki/index.php/FreeContact).
- Install [CCMpred](https://github.com/soedinglab/CCMpred).
- Install [MODELLER](https://salilab.org/modeller), which requires a license key. Only the Python package is required so this can be installed with `conda install modeller -c salilab`.
- Install [CNS](http://cns-online.org/v1.3). We found we had to follow all of the steps in [this comment](https://ask.bioexcel.eu/t/cns-errors-before-after-recompilation/54/14) to get CNS working: set `MXFPEPS2` in machvar.inc to 8192, remove `-fastm` flag in the make file, set `MXRTP` in rtf.inc in the source directory to 4000 and in machvar.f add `WRITE (6,’(I6,E10.3,E10.3)’) I, ONEP, ONEM` just above line 67, which looks like `IF (ONE .EQ. ONEP .OR. ONE .EQ. ONEM) THEN`. We also had to install the `flex-devel` package via our system package manager. In addition, you should change two values in `cns_solve_1.3/modules/nmr/readdata` to larger numbers to allow DMPfold to run on larger structures. Change the `nrestraints = 20000` line to something like `nrestraints = 50000` and the `nassign 1600` line to something like `nassign 3000`.
- Download and patch the required CNS scripts by changing into the `cnsfiles` directory and running `sh installscripts.sh`.
- Install [CD-HIT](https://github.com/weizhongli/cdhit), which is usually as simple as a clone and `make`. CD-HIT is not required if you don't need to predict the TM-score of generated models.
- Install the [legacy BLAST](https://tinyurl.com/y57hq2wo) software, in particular `formatdb`, `blastpgp` and `makemat`. We may update this to BLAST+ in the future.
- Other software is pre-compiled and included here (PSIPRED, PSICOV, various utility scripts with the code in `src`). This should run okay but may need separate compilation using the makefile if issues arise. Some other standard programs, such as csh shell, are assumed.
- Change lines 10/13-15/18/21/24 in `seq2maps.csh`, lines 11/14/17/20 in `aln2maps.csh`, lines 4/7 in `bin/runpsipredandsolvwithdb`, lines 10/13 in `run_dmpfold.sh` and lines 7/10 in `predict_tmscore.sh` to point to the installed locations of the above software. You can also set the number of cores to use in `seq2maps.csh` and `aln2maps.csh`. This sets the number of cores for HHblits, PSICOV, FreeContact and CCMpred - the script will run faster with this set to a value larger than 1 (e.g. 4 or 8).

Check the continuous integration [setup script](.travis.yml) and [logs](https://travis-ci.org/psipred/DMPfold) for additional tips and a step-by-step installation on Ubuntu.

## Usage

Here we give an example of running DMPfold on Pfam family PF10963.
First you need to generate the `.21c` and `.map` files.
This can be done in one of two ways:
- From a single sequence: `csh seq2maps.csh example/PF10963.fasta` to run HHblits, PSIPRED, SOLVPRED, PSICOV, FreeContact, CCMpred and alnstats.
- From an alignment: `csh aln2maps.csh example/PF10963.aln` to run PSIPRED, SOLVPRED, PSICOV, FreeContact, CCMpred and alnstats. The file `PF10963.aln` has one sequence per line with the ungapped target sequence as the first line.

Then run `sh run_dmpfold.sh example/PF10963.fasta PF10963.21c PF10963.map ./PF10963` to run DMPfold, where the last parameter is an output directory that will be created.
Running `sh run_dmpfold.sh example/PF10963.fasta PF10963.21c PF10963.map ./PF10963 5 20` instead runs 5 iterations with 20 models per iteration (default is 3 and 50).
The final model is `final_1.pdb` and other structures may or may not be generated as `final_2.pdb` to `final_5.pdb` if they are significantly different.
Many other files are generated totalling around 100 MB - these should be deleted to save disk space if you are running DMPfold on many sequences.

To predict the TM-score of a DMPfold model using our trained predictor, run `sh predict_tmscore.sh example/PF10963.fasta PF10963.aln PF10963/final_1.pdb PF10963/rawdistpred.1`.
If this predictor estimates that a model has a TM-score of at least 0.5 then there is an 83% chance of this being the case according to cross-validation of the Pfam validation set.

See Supplementary Figure 1 in the paper for estimations on run time.
It takes around 3 hours on a single core to carry out a complete DMPfold run for a 200 residue protein, but this can occasionally be much longer due to PSICOV not converging.
8 GB memory is generally sufficient to run DMPfold but more may be required for larger proteins.

Figure 5 in the paper gives some data on how DMPfold performs with respect to sequence length.
Sequences up to around 600 residues in length can be modelled accurately, with performance degrading above this.

## Data

Models for the 1,475 [Pfam](http://pfam.xfam.org) families modelled in the paper can be downloaded [here](http://bioinf.cs.ucl.ac.uk/downloads/dmpfold/pfam_models.tgz).
Additional models for the remainder of the dark Pfam families can be downloaded [here](http://bioinf.cs.ucl.ac.uk/downloads/dmpfold/pfam_lowconf_models.tgz) (some were not modelled due to small sequence alignments).
Alignments for the Pfam families without available templates can be downloaded [here](http://bioinf.cs.ucl.ac.uk/downloads/dmpfold/pfam_alignments.tgz).
The format is one sequence per line with the ungapped target sequence as the first line.

The directory [pfam](pfam) in this repository contains text files with the lists from Figure 4A of the paper, target sequences for modelled families and data for modelled families (sequence length, effective sequence count, distogram satisfaction scores, estimated TM-score and probability TM-score >= 0.5).

The list of PDB chains used for training can be found [here](data/train_list.txt).
