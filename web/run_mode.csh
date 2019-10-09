#!/bin/csh

# Set this to point to the DMPfold directory
set dmpfolddir = /data/DMPfold
set bindir = $dmpfolddir/bin
set seqfile = $1

set seq_count = `grep '>' $seqfile | wc -l`

if($seq_count > 1) then
  echo "/bin/csh $dmpfolddir/aln2maps.csh"
  /bin/csh $dmpfolddir/aln2maps.csh
else
  echo "/bin/csh $dmpfolddir/seq2maps.csh"
  /bin/csh $dmpfolddir/seq2maps.csh
endif
