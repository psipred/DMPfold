/* QMODOPE - Run QMODCHECK + DOPE on protein structure ensemble - by David T. Jones, February 2008 */

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <sys/wait.h>
#include <math.h>

#define MAXNSTRUC 10000

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
#define dist(a,b) sqrt(SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))
#define distsq(a,b) (SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))

extern void    *
                calloc(), *malloc(), *realloc();

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

struct distrec
{
    short           a, b;
}              *distmat, **sortlist;

int             nmodels, natoms[MAXNSTRUC];
float         **dmat;

/***************************************************************************/

void
                fail(errstr)
    char           *errstr;
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
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

int
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

    if (feof(ifp))
    {
	*natoms = 0;
	return 1;
    }
	
    while (!feof(ifp))
    {
	if (!fgets(buf, 160, ifp))
	    return 1;

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
	    if (!strncmp(atmnam, " N ", 3))
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

    return 0;
}


int main(int argc, char **argv)
{
    int             i, j, k, nres;
    int             blksize, hashval, nmatch, moda, modb, ret;
    float           x, y, z, d, r, rmsd, u[3][3], **rmsdmat, matchsum, qm_val, best_qm=1e32, d_val;
    FILE           *ifp, *ofp, *qfp, *dfp;
    char cmdstr[512];
    Point           new, CG_a, CG_b;
    Transform       fr_xf;

    if (argc != 2)
	fail("usage : collect_data ensemble.pdb");

    ifp = fopen(argv[1], "r");	/* Open PDB file in TEXT mode */

    if (!ifp)
    {
	printf("PDB file error!\n");
	exit(-1);
    }

    /* Read models */
    for (nmodels=0; nmodels<MAXNSTRUC; nmodels++)
    {
	if (read_atoms(ifp, natoms+nmodels, atoms+nmodels))
	    break;

	if (natoms[nmodels] < 10)
	    continue;

	ofp = fopen("temp.pdb", "w");
	for (i = 0; i < natoms[nmodels]; i++)
	{
	    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
		    i + 1, atoms[nmodels][i].atmnam, rnames[atoms[nmodels][i].aac], atoms[nmodels][i].resnum, atoms[nmodels][i].pos[0], atoms[nmodels][i].pos[1], atoms[nmodels][i].pos[2], 1.0, 1.0);
	}
	fprintf(ofp, "TER\nEND\n");
	fclose(ofp);

	printf("%d ", nmodels+1);

	fflush(stdout);

	sprintf(cmdstr, "python dope.scr | tail -1");

	dfp = popen(cmdstr, "r");

	if (!dfp)
	    fail("Pipe error!");

	if ((WIFSIGNALED(ret) && (WTERMSIG(ret) == SIGINT || WTERMSIG(ret) == SIGQUIT)))
	    break;

	if (fscanf(dfp, "%f", &d_val) != 1)
	    puts("FAILED!");

	pclose(dfp);

	sprintf(cmdstr, "./qmodcheck temp.pdb | tail -1 > qmod.tmp");
//	sprintf(cmdstr, "./qmodcheck qmod_temp.pdb A | tail -1");
	
	ret = system(cmdstr);

	if (ret < 0 || (WIFSIGNALED(ret) && (WTERMSIG(ret) == SIGINT || WTERMSIG(ret) == SIGQUIT)))
	    break;

	qfp = fopen("qmod.tmp", "r");
	if (!qfp)
	    fail("Cannot read qmod.tmp");
	if (fscanf(qfp, "%f", &qm_val) == 1)
	{
	    qm_val = d_val * -qm_val;
	    printf("%f\n", qm_val);
	    if (qm_val < best_qm)
	    {
		best_qm = qm_val;
		system("/bin/cp -f temp.pdb best_qdope.pdb");
	    }
	}
	else
	    puts("FAILED!");

	fclose(qfp);

	fflush(stdout);
    }

    fclose(ifp);
}
