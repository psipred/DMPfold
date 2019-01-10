/* Protein Chain Superposition - by David Jones, February 1992 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define MAXSEQLEN 5000

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

char            pdbfn[80], csdfn[80], logfn[80], keyword[40], buf[160];

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

/***************************************************************************/

void
                fail(errstr)
    char           *errstr;
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

/* This routine will read in one sequence from a database file. The
   sequence can be in any of the supported formats. Returns length
   of sequence.
*/
int
  getseq(char *dbname, char *dseq, FILE * lfil)
{
    int i, j, len;
    short badln, fformat;
    enum
    {
	unknown, embl, genbank, staden, fastp, codata, owl, intelgen, gcg
    };
    char temp[8192], split;
    int offset;

    offset = j = 0;

    if (!fgets(temp, 8192, lfil))
	return (-1);
    if (strstr(temp, "of:") != NULL && strstr(temp, "check:") != NULL)
	fformat = gcg;
    else if ((temp[0] == '<') && (temp[19] == '>'))
	fformat = staden;
    else if (strncmp(temp, "ID   ", 5) == 0)
	fformat = embl;
    else if (strncmp(temp, "LOCUS     ", 10) == 0)
	fformat = genbank;
    else if (strncmp(temp, "ENTRY", 5) == 0)
	fformat = codata;
    else if (temp[0] == ';')
	fformat = intelgen;
    else if (temp[0] == '>' && (temp[1] == '>' || temp[3] == ';'))
	fformat = owl;
    else if (temp[0] == '>')
	fformat = fastp;
    else
	fformat = unknown;

    switch (fformat)
    { 
     case gcg:
	sscanf(strstr(temp, "of:")+3, "%s", dbname);
	while (strstr(temp, "..") == NULL)
	    fgets(temp, 8192, lfil);
	fgets(temp, 8192, lfil);
	break;
   case embl:
	strncpy(dbname, temp + 5, 70);
	while (temp[0] != ' ')
	    fgets(temp, 8192, lfil);
	break;

    case genbank:
	while (strncmp(temp, "ORIGIN", 6) != 0)
	{
	    fgets(temp, 8192, lfil);
	    if (strncmp(temp, "DEFINITION", 10) == 0)
		strncpy(dbname, temp + 12, 70);
	}
	fgets(temp, 8192, lfil);
	break;

    case codata:
	strncpy(dbname, temp + 6, 70);
	while (strncmp(temp, "SEQUENCE", 8) != 0)
	    fgets(temp, 8192, lfil);
	fgets(temp, 8192, lfil);
	break;

    case owl:
	fgets(temp, 8192, lfil);
	strncpy(dbname, temp, 70);
	fgets(temp, 8192, lfil);
	break;

    case fastp:
	strncpy(dbname, temp + 1, 70);
	fgets(temp, 8192, lfil);
	break;

    case staden:
	strncpy(dbname, temp + 1, 18);
	offset = 20;
	break;

    case intelgen:
	while (*temp == ';')
	    fgets(temp, 8192, lfil);
	fgets(temp, 8192, lfil);
	break;

    default:
	do
	{
	    len = strlen(temp);
	    for (badln = i = 0; i < len; i++)
		if (islower(temp[i]) || temp[i] == 'J' || temp[i] == 'O' || temp[i] == 'U')
		{
		    badln = TRUE;
		    break;
		}
	    if (badln && !fgets(temp, 8192, lfil))
		return (-1);
	}
	while (badln);
	strcpy(dbname, "<NO NAME>");
	break;
    }

    if (dbname[(len = strlen(dbname)) - 1] == '\n')
	dbname[--len] = '\0';
    if (len >= 70)
	dbname[70] = '\0';

    for (;;)
    {
	if (!strncmp(temp, "//", 2))
	    break;
	len = strlen(temp);
	for (i = offset; i < len && j < MAXSEQLEN; i++)
	{
	    split = islower(temp[i]) ? toupper(temp[i]) : temp[i];
	    if (split == '@' || (fformat == owl && split == '*'))
	    {
		dseq[j] = '\0';
		while (fgets(temp, 8192, lfil));
		return (j);
	    }
	    if (isalpha(split))
		dseq[j++] = split;
	    else if (temp[i] == '\n')
		break;
	}
	if (staden)
	    offset = 0;
	if (!fgets(temp, 8192, lfil))
	    break;
    }

    if (j == MAXSEQLEN)
	printf("\nWARNING: sequence %s over %d long; truncated!\n",
	       dbname, MAXSEQLEN);

    dseq[j] = '\0';
    return (j);
}

/* Convert AA letter to numeric code (0-22) */
int
                aanum(int ch)
{
    const static int aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
	22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 22);
}    

int main(int argc, char **argv)
{
    int i, j, seqlen;
    float d1, d2, r;
    char            fname[160], buf[65536], seq[MAXSEQLEN+1];
    FILE           *ifp, *ofp, *rfp;
    
    if (argc != 3)
	fail("usage : contact2noe seq contactlist");

    ifp = fopen(argv[1], "r");	/* Open seq file in TEXT mode */
    if (!ifp)
    {
	printf("Seq file error!\n");
	exit(-1);
    }

    seqlen = getseq(fname, seq, ifp);

    fclose(ifp);

    for (i=0; i<seqlen; i++)
	seq[i] = aanum(seq[i]);

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

	if (isalpha(buf[0]))
	    continue;
	
	if (sscanf(buf, "%d%d%f%f%f", &i, &j, &d1, &d2, &r) != 5)
	    break;

	if (d1 < 3.5)
	    d1 = 3.5;

	if (seq[i-1] == 7 && seq[j-1] == 7)
	    printf("assign (resid %d and name ca) (resid %d and name ca) %f 0.1 %f\n", i, j, d1 + 0.1, d2 - d1 - 0.1);
	else if (seq[i-1] == 7)
	    printf("assign (resid %d and name ca) (resid %d and name cb) %f 0.1 %f\n", i, j, d1 + 0.1, d2 - d1 - 0.1);
	else if (seq[j-1] == 7)
	    printf("assign (resid %d and name cb) (resid %d and name ca) %f 0.1 %f\n", i, j, d1 + 0.1, d2 - d1 - 0.1);
	else
	    printf("assign (resid %d and name cb) (resid %d and name cb) %f 0.1 %f\n", i, j, d1 + 0.1, d2 - d1 - 0.1);
    }

    fclose(ifp);
}
