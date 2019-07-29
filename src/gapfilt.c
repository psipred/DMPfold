/* Filter gaps from alignment file - David Jones, 2011 */

/* David Jones, June 2011 */

/* Copyright (C) 2011 University College London */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define FALSE 0
#define TRUE 1

#define SQR(x) ((x)*(x))

#define MAXSEQLEN 50000

char  **aln;
int nseqs;

const char     *rescodes = "ARNDCQEGHILKMFPSTWYVXXX";

/* Dump a rude message to standard error and exit */
void
                fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

/* Convert AA letter to numeric code (0-21) */
int
                aanum(int ch)
{
    const static int aacvs[] =
    {
	999, 0, 3, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10, 12, 2,
	21, 14, 5, 1, 15, 16, 21, 19, 17, 21, 18, 6
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

int             main(int argc, char **argv)
{
    int             a, b, c, d, i, j, ii, jj, k, l, n, llen, maxnps=3, posn, qstart, seqlen, endflg = FALSE, tots, maxtots, lstart, nids, s;
    double sumi[21], sumj[21], sum, score, pab[MAXSEQLEN][MAXSEQLEN][21][21], eij[MAXSEQLEN][MAXSEQLEN][21][21], gij[MAXSEQLEN][MAXSEQLEN][21][21], pij[MAXSEQLEN][MAXSEQLEN][21], pa[MAXSEQLEN][21], pav[MAXSEQLEN][21], di, pdir, hahb, hx, hy, hxy, mi, z, oldvec[21], change, wtsum;
    float *weight, maxgapf;
    char            buf[4096], name[512], seq[MAXSEQLEN], qseq[160], *cp, ch;
    FILE *ifp;

    if (argc < 2)
	fail("Usage: gapfilt alnfile {fract}");

    if (argc > 2)
	maxgapf = atof(argv[2]);
    else
	maxgapf = 0.3;
    
    ifp = fopen(argv[1], "r");
    if (!ifp)
	fail("Unable to open alignment file!");

    for (nseqs=0;; nseqs++)
    {
	if (!fgets(seq, MAXSEQLEN, ifp))
	    break;

	a = b = 0;
	
	cp = seq;
	while (ch = *cp++)
	{
	    if (ch == '\n')
		break;
	    if (!isalpha(ch))
		a++;
	    b++;
	}

	if ((float)a/b <= maxgapf)
	    printf("%s", seq);
    }
    
    return 0;
}
