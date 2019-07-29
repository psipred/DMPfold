/* Convert aln format file to fasta format file */

/* David Jones, July 2011 */

/* Copyright (C) 2011 University College London */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#define MAXSEQLEN 1000000

#define FALSE 0
#define TRUE 1

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    GAP, UNK
};

const char     *rescodes = "ARNDCQEGHILKMFPSTWYVBZX";

/* Dump a rude message to standard error and exit */
void
                fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

/* Convert AA letter to numeric code (0-22) */
int
                aanum(int ch)
{
    static int      aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
	22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 22);
}


int main(int argc, char **argv)
{
    int             i, j, ngap, naa, seqlen = 0, seqcount=0, nid;
    char            desc[MAXSEQLEN], refseq[MAXSEQLEN], seq[MAXSEQLEN], buf[MAXSEQLEN], *p;
    FILE           *ifp;

    if (argc < 2)
	fail("Usage: aln2fasta aln-file");

    ifp = fopen(argv[1], "r");
    if (!ifp)
	fail("Unable to open aln file!");

    nid = 1;
    for (;;)
    {
	if (!fgets(buf, MAXSEQLEN, ifp))
	    break;
	printf(">SEQ%06d\n%s", nid++, buf);
    }

    return 0;
}
