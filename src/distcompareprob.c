/* Compare PDB structures to distograms */

/* David Jones, Jan 2019 */

/* Copyright (C) 2019 University College London */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAXNSTRUC 100000

#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

#define BIG ((float)1.0e10)
#define SMALL ((float)1.0e-3)
#define ZERO ((float)0.0)
#define HALF ((float)0.5)
#define ONE ((float)1.0)
#define TWO ((float)2.0)
#define THREE ((float)3.0)
#define FOUR ((float)4.0)
#define PI ((float)3.1415927)

#define NTOKENS 26
#define DEFR 1.80

/* A list of common PDB record types... */

#define HEADER             1
#define COMPND             2
#define SOURCE             3
#define AUTHOR             4
#define REVDAT             5
#define REMARK             6
#define SEQRES             7
#define FTNOTE             8
#define HET                9
#define FORMUL             10
#define HELIX              11
#define CRYST1             12
#define ORIGX1             13
#define ORIGX2             14
#define ORIGX3             15
#define SCALE1             16
#define SCALE2             17
#define SCALE3             18
#define ATOM               19
#define TER                20
#define HETATM             21
#define CONECT             22
#define ENDENT             23
#define JRNL               24
#define TURN               25
#define ENDMDL             26

#define MEMGRAIN	   100

#define dotprod(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define vecprod(a,b,c) (a[0]=b[1]*c[2]-b[2]*c[1],a[1]=b[2]*c[0]-b[0]*c[2],a[2]=b[0]*c[1]-b[1]*c[0])
#define veczero(v) memset(v, 0, sizeof(v))
#define veccopy(a,b) (a[0]=b[0],a[1]=b[1],a[2]=b[2])
#define vecadd(a,b,c) (a[0]=b[0]+c[0],a[1]=b[1]+c[1],a[2]=b[2]+c[2])
#define vecsub(a,b,c) (a[0]=b[0]-c[0],a[1]=b[1]-c[1],a[2]=b[2]-c[2])
#define vecscale(a,b,c) (a[0]=b[0]*c,a[1]=b[1]*c,a[2]=b[2]*c)
#define ran0() ((rand()&32767)/(float)32768.0)
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define SQR(x) ((x)*(x))
#define distsq(a,b) (SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))
#define dist(a,b) sqrt(distsq(a,b))

typedef float   Transform[4][4];
typedef float   Point[3];

struct pdbatm
{
    Point           pos;
    float           radius;	/* Radius of sphere */
    int             resnum, aac;
    char            ispolar, isback, donors, acceptors;
    char            ndon, nacc, donflag, accflag;
    char            atmnam[5];
}              *atoms[MAXNSTRUC];

char            pdbfn[80], csdfn[80], logfn[80], keyword[40], buf[512];

/* Record names for decoding record types */
char           *tokstr[] =
{
 "HEADER", "COMPND", "SOURCE", "AUTHOR", "REVDAT",
 "REMARK", "SEQRES", "FTNOTE", "HET", "FORMUL",
 "HELIX", "CRYST1", "ORIGX1", "ORIGX2", "ORIGX3",
 "SCALE1", "SCALE2", "SCALE3", "ATOM", "TER",
 "HETATM", "CONECT", "END", "JRNL", "TURN", "ENDMDL"
};

/* Residue name to allow conversion of a.a. name into numeric code */
char           *rnames[] =
{
 "ALA", "ARG", "ASN", "ASP", "CYS",
 "GLN", "GLU", "GLY", "HIS", "ILE",
 "LEU", "LYS", "MET", "PHE", "PRO",
 "SER", "THR", "TRP", "TYR", "VAL",
 "ASX", "GLX", "UNK"
};

enum RESCODE
{
    Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys,
    Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val,
    Asx, Glx, Unk
};

int             nmodels, natoms[MAXNSTRUC], cbflg;

/***************************************************************************/

void
                fail(errstr)
    char           *errstr;
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

/* Allocate matrix */
void           *allocmat(int rows, int columns, int size, int clrflg)
{
    int             i;
    void          **p;

    p = (void **) malloc(rows * sizeof(void *));

    if (p == NULL)
	fail("allocmat: malloc [] failed!");
    if (clrflg)
    {
	for (i = 0; i < rows; i++)
	    if ((p[i] = calloc(columns, size)) == NULL)
		fail("allocmat: calloc [][] failed!");
    }
    else
	for (i = 0; i < rows; i++)
	    if ((p[i] = malloc(columns * size)) == NULL)
		fail("allocmat: malloc [][] failed!");

    return p;
}

void            freemat(void *p, int rows)
{
    int             i;

    for (i = rows - 1; i >= 0; i--)
	free(((void **) p)[i]);
    free(p);
}

/* Get atomic coordinates from ATOM records */
void
                getcoord(x, y, z, chain, n, aacode, atmnam)
    float          *x, *y, *z;
    int            *n, *aacode;
    char           *atmnam, *chain;
{
    int             i, atomn;

    strncpy(atmnam, buf + 12, 4);
    atmnam[4] = '\0';
    if (atmnam[2] == ' ')
	atmnam[3] = ' ';
    sscanf(buf + 6, "%d", &atomn);
    *chain = buf[21];

    sscanf(buf + 22, "%d%f%f%f", n, x, y, z);
    for (i = 0; i < 22; i++)
	if (!strncmp(rnames[i], buf + 17, 3))
	    break;
    *aacode = i;
}

void
                read_atoms(ifp, natoms, atmptr)
    FILE           *ifp;
    int            *natoms;
    struct pdbatm **atmptr;
{
    int             i, token=ENDENT, namino, aac, resnum = 0, blksize;
    float           x, y, z, u[3][3];
    char            atmnam[5], chain = '?';
    struct pdbatm  *atom = *atmptr;

    blksize = 10 * MEMGRAIN;
    if (!(atom = malloc(blksize * sizeof(struct pdbatm))))
	fail("read_atoms : Out of memory!");

    while (!feof(ifp))
    {
	if (!fgets(buf, 512, ifp))
	    break;

	/* Read the record name */
	if (sscanf(buf, "%s", keyword) != 1)
	    break;

	if (!keyword[0])	/* Odd - there isn't a record name! Exit. */
	    break;

	token = 0;
	for (i = 1; i <= NTOKENS; i++)	/* Decode record type */
	    if (!strcmp(keyword, tokstr[i - 1]))
		token = i;

	if (token == ENDENT || token == ENDMDL)
	    break;

	switch (token)
	{
	case ATOM:
	    if (buf[16] != ' ' && buf[16] != 'A')
		continue;
	    buf[26] = ' ';
	    if (*natoms >= blksize)
	    {
		blksize += MEMGRAIN;
		atom = realloc(atom, blksize * sizeof(struct pdbatm));
		if (!atom)
		    fail("read_atoms : Out of Memory!");
	    }
	    getcoord(&x, &y, &z, &chain, &namino, &aac, atmnam);

	    if (!cbflg && strncmp(atmnam, " CA", 3))
		continue;

	    if (cbflg && aac == Gly && strncmp(atmnam, " CA", 3))
		continue;
	    if (cbflg && aac != Gly && strncmp(atmnam, " CB", 3))
		continue;

	    resnum++;
	    atom[*natoms].pos[0] = x;
	    atom[*natoms].pos[1] = y;
	    atom[*natoms].pos[2] = z;
	    atom[*natoms].resnum = resnum;
	    atom[*natoms].aac = aac;
	    atom[*natoms].isback = i > 0 && i < 4;
	    strcpy(atom[*natoms].atmnam, atmnam);
	    ++*natoms;
	    break;
	default:		/* Ignore all other types in this version */
	    break;
	}
    }

    *atmptr = atom;
}

int main(int argc, char **argv)
{
    int             i, j, k, ii, jj, kk, l, n, nmax = 0, at1, at2, topomin=23, cseqsep;
    int             blksize, hashval, moda, modb, bestmod;
    unsigned int    nmatch[MAXNSTRUC], clustflg[MAXNSTRUC], filepos[MAXNSTRUC], modnum = 1;
    float           x, y, z, d, r, mean, sd, rmsd, u[3][3], **rmsdmat, matchsum, wt2 = 0.01, cutoff = 0.0, dcut, (**dpmat)[20], diff, probs[20], bestmatch = 1e32;
    char            fname[160], buf[4096];
    FILE           *ifp, *ofp, *rfp, *cfp;
    Point           new, CG_a, CG_b;
    Transform       fr_xf;

    if (argc < 2)
	fail("usage : distcompareprob ensemble.pdb distogramlist {param1}");

    ifp = fopen(argv[2], "r");	/* Open contact file in TEXT mode */
    if (!ifp)
    {
	printf("Contact file error!\n");
	exit(-1);
    }

    while (!feof(ifp))
    {
	if (!fgets(buf, 256, ifp))
	    break;

	if (!strncmp(buf, "PFRMAT RR", 9))
	{
	    cbflg = 1;
	    break;
	}
    }
    
    fclose(ifp);

    ifp = fopen(argv[1], "r");	/* Open PDB file in TEXT mode */
    if (!ifp)
    {
	printf("PDB file error!\n");
	exit(-1);
    }

    if (argc > 3)
	topomin = atoi(argv[3]);

    if (argc > 4)
	wt2 = cutoff = atof(argv[4]);

    /* Read models */
    for (nmodels=0; nmodels<MAXNSTRUC; nmodels++)
    {
	filepos[nmodels] = ftell(ifp);
	read_atoms(ifp, natoms+nmodels, atoms+nmodels);

	if (!natoms[nmodels])
	    break;
    }

    fprintf(stderr, "%d models read from PDB file\n", nmodels);

    dpmat = allocmat(natoms[0], natoms[0], sizeof(float[20]), TRUE);

    cfp = fopen(argv[2], "r");	/* Open contact file in TEXT mode */
    if (!cfp)
    {
	printf("Contact file error!\n");
	exit(-1);
    }

    n = cseqsep = 0;
    
    while (!feof(cfp))
    {
	if (!fgets(buf, 1024, cfp))
	    break;

	if (sscanf(buf, "%d%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f", &i, &j, probs, probs+1, probs+2, probs+3, probs+4, probs+5, probs+6, probs+7, probs+8, probs+9, probs+10, probs+11, probs+12, probs+13, probs+14, probs+15, probs+16, probs+17, probs+18, probs+19) != 22)
	    break;

	for (k=0; k<20; k++)
	    dpmat[i-1][j-1][k] = dpmat[j-1][i-1][k] = probs[k];

	cseqsep += abs(i-j);
	n++;
    }
    
    fclose(cfp);

    cseqsep /= n;

    topomin = MIN(cseqsep, topomin);
    
    for (moda=0; moda<nmodels; moda++)
    {
/*	printf("%d %d\n", moda, natoms[0]); */
	for (n = rmsd = matchsum = i = 0; i < natoms[moda]; i++)
	    for (j = i+topomin; j < natoms[moda]; j++)
	    {
/*		printf("%d %d %f\n", i, j, dist(atoms[moda][i].pos, atoms[moda][j].pos)); */

		if (dpmat[i][j][0])
		{
		    d = dist(atoms[moda][i].pos, atoms[moda][j].pos);

		    if (d < 8.0F)
		    {
			k = MAX(0, (int)(0.5F * (d - 4.0F)));
			matchsum -= dpmat[i][j][k];
			n++;
		    }
		    else if (d < 13.0F)
		    {
			k = MIN(19, (int)d);
			matchsum -= dpmat[i][j][k];
			n++;
		    }
		}
	    }

	printf("%d %f %f\n", moda+1, matchsum, matchsum/n);
	if (matchsum/n < bestmatch)
	{
	    bestmod = moda;
	    bestmatch = matchsum/n;
	}
    }

    if (nmodels > 1)
    {
	ofp = fopen("best_dist.pdb", "w");
	if (!ofp)
	    fail("Cannot open PDB file writing!");
	
	fprintf(ofp, "REMARK  ORIGINAL MODEL NUMBER    %4d\n", bestmod+1);
	fseek(ifp, filepos[bestmod], SEEK_SET);
	
	while (!feof(ifp))
	{
	    if (!fgets(buf, 160, ifp))
		break;
	    fputs(buf, ofp);
	    if (!strncmp(buf, "END", 3))
		break;
	}
	
	fclose(ifp);
	fclose(ofp);
    }

    return 0;
}
