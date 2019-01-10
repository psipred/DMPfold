/* METAPSICOV - Neural Network Prediction of Contacts */

/* Copyright (C) 2000 David T. Jones - Created : January 2000 */
/* Original Neural Network code Copyright (C) 1990 David T. Jones */

/* Average Prediction Module */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>


#define MAXSEQLEN 5000

/* logistic 'squashing' function (+/- 1.0) */
#define logistic(x) (1.0F / (1.0F + (float)exp(-(x))))

#define MAX(x,y) ((x)>(y)?(x):(y))

char           *wtfnm;

float **cmutmat;

int             profile[MAXSEQLEN][20];

char seq[MAXSEQLEN];

struct entry
{
    char           *id, *seq, **contmap;
    float         **mimat, **mipmat, **potmat, **psicov, **evfold, **ccmpred, *cprob, *hprob, *eprob, *entropy, *psolv, **profile, effnseq;
    float         *aacomp, cmean, hmean, emean, smean, entropmean;
    int           length, nseq;
} target;

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    UNK
};

void
err(char *s)
{
    fprintf(stderr, "%s\n", s);
}

void
fail(char *s)
{
    err(s);
    exit(1);
}

/* Allocate matrix */
void           *allocmat(int rows, int columns, int size)
{
    int             i;
    void          **p, *rp;

    rp = malloc(rows * sizeof(void *) + sizeof(int));

    if (rp == NULL)
	fail("allocmat: malloc [] failed!");

    *((int *)rp) = rows;

    p = rp + sizeof(int);

    for (i = 0; i < rows; i++)
	if ((p[i] = calloc(columns, size)) == NULL)
	    fail("allocmat: malloc [][] failed!");

    return p;
}

/* Allocate vector */
void           *allocvec(int columns, int size)
{
    void          *p;

    p = calloc(columns, size);

    if (p == NULL)
	fail("allocvec: calloc failed!");

    return p;
}

/* Convert AA letter to numeric code (0-20) */
int
aanum(ch)
    int             ch;
{
    static const int      aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2,
	20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

/* Make 1st level prediction averaged over specified weight sets */
void
gen_output(char *mapfname, char *fixfname)
{
    int             aa, i, j, k, n, winpos, winpos2, midpos, ws, seqsep, nres;
    float val;
    double prob, x;
    FILE *ofp;

    ofp = fopen(fixfname, "wb");
    if (!ofp)
	fail("Fixed param file open error!");

    val = target.length;
    if (fwrite(&val, sizeof(float), 1, ofp) != 1)
	fail("Fixed param file write!");

    val = target.nseq;
    if (fwrite(&val, sizeof(float), 1, ofp) != 1)
	fail("Fixed param file write!");

    if (fwrite(&target.effnseq, sizeof(float), 1, ofp) != 1)
	fail("Bad fixed param file write!");

    if (fwrite(&target.cmean, sizeof(float), 1, ofp) != 1)
	fail("Bad fixed param file write!");

    if (fwrite(&target.hmean, sizeof(float), 1, ofp) != 1)
	fail("Bad fixed param file write!");
    
    if (fwrite(&target.emean, sizeof(float), 1, ofp) != 1)
	fail("Bad fixed param file write!");

    if (fwrite(&target.smean, sizeof(float), 1, ofp) != 1)
	fail("Bad fixed param file write!");

    if (fwrite(&target.entropmean, sizeof(float), 1, ofp) != 1)
	fail("Bad fixed param file write!");

    if (fwrite(target.aacomp, sizeof(float), 21, ofp) != 21)
	fail("Bad fixed param file write!");
    
    fclose(ofp);
    
    ofp = fopen(mapfname, "wb");
    
    if (!ofp)
	fail("Map file open error!");
    
    nres = target.length;
    
    for (k=0; k<21; k++)
	for (i=0; i<nres; i++)
	    for (j=0; j<nres; j++)
	    {
		if (fwrite(&target.profile[i][k], sizeof(float), 1, ofp) != 1)
		    fail("Map write error!");
	    }
    
    for (k=0; k<21; k++)
	for (i=0; i<nres; i++)
	    for (j=0; j<nres; j++)
	    {
		if (fwrite(&target.profile[j][k], sizeof(float), 1, ofp) != 1)
		    fail("Map write error!");
	    }
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.mimat[i][j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.mipmat[i][j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.potmat[i][j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.psicov[i][j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.evfold[i][j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.ccmpred[i][j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.cprob[i], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.cprob[j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.hprob[i], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.hprob[j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.eprob[i], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.eprob[j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.entropy[i], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.entropy[j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.psolv[i], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}

    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    if (fwrite(&target.psolv[j], sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}

    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	{
	    val = log1p((float)abs(j-i));
	    if (j < i)
		val = -val;
	    if (fwrite(&val, sizeof(float), 1, ofp) != 1)
		fail("Map write error!");
	}
    
    val = 1.0F;
    for (i=0; i<nres; i++)
	for (j=0; j<nres; j++)
	    if (fwrite(&val, sizeof(float), 1, ofp) != 1)
		fail("Map write error!");

    fclose(ofp);
}

/* Read PSI AA frequency data */
int             getmtx(FILE *lfil)
{
    int             aa, i, j, naa;
    char            buf[256], *p;

    if (fscanf(lfil, "%d", &naa) != 1)
      fail("Bad mtx file - no sequence length!");

    if (naa > MAXSEQLEN)
      fail("Input sequence too long!");

    if (fscanf(lfil, "%s", seq) != 1)
      fail("Bad mtx file - no sequence!");

    while (!feof(lfil))
      {
	if (!fgets(buf, 65536, lfil))
	  fail("Bad mtx file!");
	if (!strncmp(buf, "-32768 ", 7))
	  {
	    for (j=0; j<naa; j++)
	      {
		if (sscanf(buf, "%*d%d%*d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%*d%d", &profile[j][ALA],  &profile[j][CYS], &profile[j][ASP],  &profile[j][GLU],  &profile[j][PHE],  &profile[j][GLY],  &profile[j][HIS],  &profile[j][ILE],  &profile[j][LYS],  &profile[j][LEU],  &profile[j][MET],  &profile[j][ASN],  &profile[j][PRO],  &profile[j][GLN],  &profile[j][ARG],  &profile[j][SER],  &profile[j][THR],  &profile[j][VAL],  &profile[j][TRP],  &profile[j][TYR]) != 20)
		  fail("Bad mtx format!");
		aa = aanum(seq[j]);
		if (!fgets(buf, 65536, lfil))
		  break;
	      }
	  }
      }

    return naa;
}

/* Read in data */
void            read_dat(char *colname, char *pairname, char *psiname, char *evfname, char *ccmpfname, char *ssname, char *solvname)
{
    FILE           *ifp, *tfp, *cfp, *pfp, *sfp;
    char            buf[512];
    int             i, j, k, l, m, ncon, nncon, nres, npairs, ctoggle = 0;
    float val, consv[22], aventropy;

    if (!(cfp = fopen(colname, "r")))
	fail("Cannot open column stats file!");
    
    if (!fgets(buf, 512, cfp))
	fail("Bad column stats file!");
    
    sscanf(buf, "%d", &nres);
    
    target.length = nres;
    target.cprob = allocvec(nres, sizeof(float));
    target.hprob = allocvec(nres, sizeof(float));
    target.eprob = allocvec(nres, sizeof(float));
    target.entropy = allocvec(nres, sizeof(float));
    target.psolv = allocvec(nres, sizeof(float));
    target.aacomp = allocvec(21, sizeof(float));
    target.profile = allocmat(nres, 21, sizeof(float));

    npairs = nres * (nres - 1) / 2;
    
    if (!fgets(buf, 512, cfp))
	fail("Bad column stats file!");
    
    sscanf(buf, "%d", &target.nseq);
    
    if (!fgets(buf, 512, cfp))
	fail("Bad column stats file!");
    
    sscanf(buf, "%f", &target.effnseq);
    
    if (!fgets(buf, 512, cfp))
	fail("Bad column stats file!");
    
    if (sscanf(buf, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
	       &consv[0], &consv[1], &consv[2], &consv[3], &consv[4],
	       &consv[5], &consv[6], &consv[7], &consv[8], &consv[9],
	       &consv[10], &consv[11], &consv[12], &consv[13],
	       &consv[14], &consv[15], &consv[16], &consv[17],
	       &consv[18], &consv[19], &consv[20]) != 21)
	fail("Bad AAcomp records in column stats file!");
    
    for (j = 0; j < 21; j++)
	target.aacomp[j] = consv[j];
    
    for (i=0; i<nres; i++)
    {
	if (!fgets(buf, 512, cfp))
	    fail("Bad column stats file!");
	
	if (sscanf(buf, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
		   &consv[0], &consv[1], &consv[2], &consv[3], &consv[4],
		   &consv[5], &consv[6], &consv[7], &consv[8], &consv[9],
		   &consv[10], &consv[11], &consv[12], &consv[13],
		   &consv[14], &consv[15], &consv[16], &consv[17],
		   &consv[18], &consv[19], &consv[20], &consv[21]) != 22)
	{
	    puts(buf);
	    fail("Bad profile records in column stats file!");
	}
	
	for (aventropy = j = 0; j < 21; j++)
	    target.profile[i][j] = consv[j];
	aventropy += (target.entropy[i] = consv[21]);
    }
    
    fclose(cfp);

    aventropy /= nres;
    
    target.entropmean = aventropy;
    
    if (!(pfp = fopen(pairname, "r")))
	fail("Cannot open pair stats file!");
    
    target.potmat = allocmat(nres, nres, sizeof(float));
    target.mimat = allocmat(nres, nres, sizeof(float));
    target.mipmat = allocmat(nres, nres, sizeof(float));
    
    for (k=0; k<npairs; k++)
    {
	if (!fgets(buf, 512, pfp))
	    fail("Bad pair stats file!");
	
	if (sscanf(buf, "%d%d%f%f%f", &i, &j, &consv[0], &consv[1], &consv[2]) != 5)
	    fail("Bad pair stats file!");
	
	target.potmat[i-1][j-1] = target.potmat[j-1][i-1] = consv[0];
	target.mimat[i-1][j-1] = target.mimat[j-1][i-1] = consv[1];
	target.mipmat[i-1][j-1] = target.mipmat[j-1][i-1] = consv[2];
    }
    
    fclose(pfp);

    target.psicov = allocmat(nres, nres, sizeof(float));
    
    if ((pfp = fopen(psiname, "r")))
    {
	for (;;)
	{
	    if (!fgets(buf, 512, pfp))
		break;
	    
	    if (sscanf(buf, "%d%d%*s%*s%f", &i, &j, &consv[0]) != 3)
		fail("Bad PSICOV file!");
	    
	    target.psicov[i-1][j-1] = target.psicov[j-1][i-1] = consv[0];
	}
    
	fclose(pfp);
    }
    
    target.evfold = allocmat(nres, nres, sizeof(float));
    
    if ((pfp = fopen(evfname, "r")))
    {
	for (;;)
	{
	    if (!fgets(buf, 512, pfp))
		break;
	    
	    if (sscanf(buf, "%d%*s%d%*s%*s%f", &i, &j, &consv[0]) != 3)
		fail("Bad evfold file!");
	    
	    target.evfold[i-1][j-1] = target.evfold[j-1][i-1] = consv[0];
	}
    
	fclose(pfp);
    }
    
    target.ccmpred = allocmat(nres, nres, sizeof(float));

    if ((pfp = fopen(ccmpfname, "r")))
    {
	if (fscanf(pfp, "%f", &consv[0]) == 1)
	{
	    rewind(pfp);
	    
	    for (i=0; i<nres; i++)
		for (j=0; j<nres; j++)
		{
		    if (fscanf(pfp, "%f", &consv[0]) != 1)
			fail("Bad CCMPRED file!");
		    
		    /* NOTE: Reading full NxN matrix for CCMPRED! */
		    target.ccmpred[i][j] = consv[0];
		}
	}
    
	fclose(pfp);
    }
    
    if (!(sfp = fopen(ssname, "r")))
	fail("Cannot open SS2 file!");
    
    if (!fgets(buf, 512, sfp) || !fgets(buf, 512, sfp))
	fail("Bad SS2 file!");
    
    i = 0;
    while (!feof(sfp))
    {
	if (!fgets(buf, 512, sfp))
	    break;
	if (sscanf(buf, "%*s%*s%*s%f%f%f", target.cprob+i, target.hprob+i, target.eprob+i) != 3)
	    fail("Bad SS2 record!");
	target.cmean += target.cprob[i];
	target.hmean += target.hprob[i];
	target.emean += target.eprob[i];
	i++;
    }
    
    fclose(sfp);
    
    nres = i;
    
    if (!(sfp = fopen(solvname, "r")))
	fail("Cannot open solv file!");
    
    i = 0;
    while (!feof(sfp))
    {
	if (!fgets(buf, 512, sfp))
	    break;
	if (sscanf(buf, "%*s%*s%f", target.psolv+i) != 1)
	    fail("Bad solv record!");
	target.smean += target.psolv[i];
	i++;
    }
    
    fclose(sfp);
    
    if (i != nres)
	fail("Solv length mismatch!");
    
    target.cmean /= (float)nres;
    target.hmean /= (float)nres;
    target.emean /= (float)nres;
    target.smean /= (float)nres;
}

main(int argc, char **argv)
{
    int             i, j, niters;
    FILE *ifp;

    /* malloc_debug(3); */
    if (argc < 10)
	fail("usage : deepmetapsicov colstats-file pairstats-file psicov-file evfold-file ccmpred-file ss2-file solv-file mapoutputfile fixedparamoutfile");
    
    read_dat(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);

    gen_output(argv[8], argv[9]);

    return 0;
}
