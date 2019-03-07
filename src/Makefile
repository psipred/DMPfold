######################################################################
#	Makefile for DMPfold utilities                               #
######################################################################

#### CC		Command to execute C compiler
#### CFLAGS	Flags to pass to C compiler.

INCDIR		= .
CC		= cc

CFLAGS		= -O
LIBS		= -lm

all:		aln2fasta contact2noe cov21stats fasta2tlc hbond2noe hbond2ssnoe pdbhfilt qmodcheck qmodcheck_mainens tmclust alnstats dmp_makepredmaps qmodope_mainens distcompareprob

install:
		/bin/cp aln2fasta contact2noe cov21stats fasta2tlc hbond2noe hbond2ssnoe pdbhfilt qmodcheck qmodcheck_mainens tmclust alnstats dmp_makepredmaps qmodope_mainens ../bin

clean:
		/bin/rm -f aln2fasta contact2noe cov21stats fasta2tlc hbond2noe hbond2ssnoe pdbhfilt qmodcheck qmodcheck_mainens tmclust alnstats dmp_makepredmaps qmodope_mainens *.o

psipred:	sspred_avpred.c ssdefs.h sspred_net.h
		$(CC) $(CFLAGS) sspred_avpred.c $(LIBS) -o psipred

psipass2:	sspred_hmulti.c ssdefs.h sspred_net2.h
		$(CC) $(CFLAGS) sspred_hmulti.c $(LIBS) -o psipass2

chkparse:	chkparse.c
		$(CC) $(CFLAGS) chkparse.c $(LIBS) -o chkparse

seq2mtx:	seq2mtx.c
		$(CC) $(CFLAGS) seq2mtx.c $(LIBS) -o seq2mtx

pfilt:		pfilt.c
		$(CC) $(CFLAGS) pfilt.c $(LIBS) -o pfilt

aln2fasta:	aln2fasta.c
		$(CC) $(CFLAGS) aln2fasta.c $(LIBS) -o aln2fasta

contact2noe:	contact2noe.c
		$(CC) $(CFLAGS) contact2noe.c $(LIBS) -o contact2noe

cov21stats:	cov21stats.c
		$(CC) $(CFLAGS) cov21stats.c $(LIBS) -o cov21stats

fasta2tlc:	fasta2tlc.c
		$(CC) $(CFLAGS) fasta2tlc.c $(LIBS) -o fasta2tlc

hbond2noe:	hbond2noe.c
		$(CC) $(CFLAGS) hbond2noe.c $(LIBS) -o hbond2noe

hbond2ssnoe:	hbond2ssnoe.c
		$(CC) $(CFLAGS) hbond2ssnoe.c $(LIBS) -o hbond2ssnoe

pdbhfilt:	pdbhfilt.c
		$(CC) $(CFLAGS) pdbhfilt.c $(LIBS) -o pdbhfilt

qmodcheck:	qmodcheck.c
		$(CC) $(CFLAGS) qmodcheck.c $(LIBS) -o qmodcheck

qmodcheck_mainens:	qmodcheck_mainens.c
		$(CC) $(CFLAGS) qmodcheck_mainens.c $(LIBS) -o qmodcheck_mainens

tmclust:	tmclust.c
		$(CC) $(CFLAGS) tmclust.c $(LIBS) -o tmclust

alnstats:	alnstats.c
		$(CC) $(CFLAGS) alnstats.c $(LIBS) -o alnstats

dmp_makepredmaps:	dmp_makepredmaps.c
		$(CC) $(CFLAGS) dmp_makepredmaps.c $(LIBS) -o dmp_makepredmaps

qmodope_mainens:	qmodope_mainens.c
		$(CC) $(CFLAGS) qmodope_mainens.c $(LIBS) -o qmodope_mainens

distcompareprob:	distcompareprob.c
		$(CC) $(CFLAGS) distcompareprob.c $(LIBS) -o distcompareprob
