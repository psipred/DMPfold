#!/bin/bash

# Predict the TM-score of a model generated using DMPfold
# Arguments are target sequence, aln file, model to score and initial distance predictions

# By Joe Greener, Jun 2019

# Copyright (C) 2019 University College London

# License: GPLv3


# Set this to point to the DMPfold directory
dmpfolddir=~/DMPfold

# Set this to point to the CD-HIT command
cdhitcmd=~/cdhit/cd-hit

targseq=$1
alnfile=$2
model=$3
distpred=$4

bindir=$dmpfolddir/bin
target=$(basename "${1%.*}")

# Get number of residues
nres=$(grep -v "^>" $targseq | tr -d "\n" | wc -c)

# Calculate effective number of sequences in alignment using cd-hit
awk -v m=$target 'BEGIN {count=1} {printf (">%s_%d\n", m, count); print; count++}' $alnfile | sed -e 's/-//g' > $target.afa
$cdhitcmd -i $target.afa -o $target.cdhit.out -T 4 -c 0.62 -n 4 > /dev/null
neff=$(grep -Ec "^>" $target.cdhit.out.clstr)
rm $target.afa $target.cdhit.out $target.cdhit.out.clstr

# Calculate satisfaction of distance constraint scores
dpscore=$($bindir/distcompareprob $model $distpred | awk '{print $2}')
dpscorenorm=$($bindir/distcompareprob $model $distpred | awk '{print $3}')

# Predict TM-score of final model
python3 $dmpfolddir/nn/dg_ema/pytorch_dgema_pred.py $nres $neff $dpscorenorm $dpscore > $target.ema.txt
tmest=$(awk '{print $11}' $target.ema.txt)
tmround=$(python3 -c "print(round($tmest, 2))")

echo "nres is $nres"
echo "neff is $neff"
echo "dpscore is $dpscore"
echo "dpscorenorm is $dpscorenorm"
echo "Probability of each TM-score bin (cols 1-10), estimated TM-score (col 11) and prob TM > 0.5 (col 12) written to $target.ema.txt"
echo "Estimated TM-score of $model is $tmround"
