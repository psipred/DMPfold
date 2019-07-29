/* TM-Score-based Protein Chain Clustering Program by David Jones, May 2004 */

/* Copyright (C) 2004 University College London */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define MAXNSTRUC 50000

#ifndef FALSE
#define TRUE 1
#define FALSE 0
#endif

#define BIG ((float)1.0e10)
#define SMALL ((float)1.0e-5)
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

/* Apply a Transform matrix to a point */
void
                transform_point(transform, p, tp)
    Transform       transform;	/* transform to apply to the point */
    Point           p;		/* the point to transform */
    Point           tp;		/* the returned point after transformation */
{
    int             i, j;
    Point           temp;

    temp[0] = p[0] + transform[0][3];
    temp[1] = p[1] + transform[1][3];
    temp[2] = p[2] + transform[2][3];
    tp[0] = dotprod(transform[0], temp);
    tp[1] = dotprod(transform[1], temp);
    tp[2] = dotprod(transform[2], temp);

#if 0
    printf("%g %g %g\n%g %g %g\n\n", p[0], p[1], p[2], tp[0], tp[1], tp[2]);
#endif
}



/* Eigen decomposition code for symmetric 6x6 matrix, taken/modified from the public
   domain Java Matrix library JAMA. */

#define NDIM 6

#define hypot2(x, y) (sqrt((x)*(x)+(y)*(y)))

/* Symmetric Householder reduction to tridiagonal form. */

void tred2(double V[NDIM][NDIM], double d[NDIM], double e[NDIM]) 
{
    int i, j, k;
    
    for (j = 0; j < NDIM; j++)
	d[j] = V[NDIM-1][j];

    /* Householder reduction to tridiagonal form. */

    for (i = NDIM-1; i > 0; i--)
    {
	/* Scale to avoid under/overflow. */

	double scale = 0.0;
	double h = 0.0;
	
	for (k = 0; k < i; k++)
	    scale = scale + fabs(d[k]);
	
	if (scale == 0.0)
	{
	    e[i] = d[i-1];
	    
	    for (j = 0; j < i; j++)
	    {
		d[j] = V[i-1][j];
		V[i][j] = 0.0;
		V[j][i] = 0.0;
	    }
	}
	else
	{
	    /* Generate Householder vector. */

	    for (k = 0; k < i; k++)
	    {
		d[k] /= scale;
		h += d[k] * d[k];
	    }
	    
	    double f = d[i-1];
	    double g = sqrt(h);
	    
	    if (f > 0)
		g = -g;
	    
	    e[i] = scale * g;
	    h = h - f * g;
	    d[i-1] = f - g;
	    
	    for (j = 0; j < i; j++)
		e[j] = 0.0;

	    /* Apply similarity transformation to remaining columns. */

	    for (j = 0; j < i; j++)
	    {
		f = d[j];
		V[j][i] = f;
		g = e[j] + V[j][j] * f;

		for (k = j+1; k <= i-1; k++)
		{
		    g += V[k][j] * d[k];
		    e[k] += V[k][j] * f;
		}
		
		e[j] = g;
	    }
	    
	    f = 0.0;
	    
	    for (j = 0; j < i; j++)
	    {
		e[j] /= h;
		f += e[j] * d[j];
	    }
	    
	    double hh = f / (h + h);
	    
	    for (j = 0; j < i; j++)
		e[j] -= hh * d[j];
	    
	    for (j = 0; j < i; j++)
	    {
		f = d[j];
		g = e[j];
		
		for (k = j; k <= i-1; k++)
		    V[k][j] -= (f * e[k] + g * d[k]);
		
		d[j] = V[i-1][j];
		
		V[i][j] = 0.0;
	    }
	}
	
	d[i] = h;
    }
    

    /* Accumulate transformations. */

    for (i = 0; i < NDIM-1; i++)
    {
	V[NDIM-1][i] = V[i][i];
	
	V[i][i] = 1.0;
	
	double h = d[i+1];
	
	if (h != 0.0)
	{
	    for (k = 0; k <= i; k++)
		d[k] = V[k][i+1] / h;
	    
	    for (j = 0; j <= i; j++)
	    {
		double g = 0.0;
		
		for (k = 0; k <= i; k++)
		    g += V[k][i+1] * V[k][j];
		
		for (k = 0; k <= i; k++)
		    V[k][j] -= g * d[k];
	    }
	}
	
	for (k = 0; k <= i; k++)
	    V[k][i+1] = 0.0;
    }
    
    for (j = 0; j < NDIM; j++)
    {
	d[j] = V[NDIM-1][j];
	V[NDIM-1][j] = 0.0;
    }
    
    V[NDIM-1][NDIM-1] = 1.0;
    
    e[0] = 0.0;
}
 

/* Symmetric tridiagonal QL algorithm. */

void tql2(double V[NDIM][NDIM], double d[NDIM], double e[NDIM]) 
{
    int i, j, k, l, m;
    double f = 0.0;
    double tst1 = 0.0;
    double eps = DBL_EPSILON;
    
    for (i = 1; i < NDIM; i++)
	e[i-1] = e[i];
    
    e[NDIM-1] = 0.0;
    
    for (l = 0; l < NDIM; l++)
    {
	/* Find small subdiagonal element */

	tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
	
	int m = l;
	
	while (m < NDIM)
	{
	    if (fabs(e[m]) <= eps*tst1)
		break;
	    
	    m++;
	}
	

	/* If m == l, d[l] is an eigenvalue, otherwise, iterate. */

	if (m > l)
	{
	    int iter = 0;
	    
	    do
	    {
		iter = iter + 1;

		/* Compute implicit shift */

		double g = d[l];
		double p = (d[l+1] - g) / (2.0 * e[l]);
		double r = hypot2(p,1.0);
		
		if (p < 0)
		    r = -r;
		
		d[l] = e[l] / (p + r);
		d[l+1] = e[l] * (p + r);
		
		double dl1 = d[l+1];
		double h = g - d[l];
		
		for (i = l+2; i < NDIM; i++)
		    d[i] -= h;
		
		f = f + h;

		/* Implicit QL transformation. */

		p = d[m];
		
		double c = 1.0;
		double c2 = c;
		double c3 = c;
		double el1 = e[l+1];
		double s = 0.0;
		double s2 = 0.0;
		
		for (i = m-1; i >= l; i--)
		{
		    c3 = c2;
		    c2 = c;
		    s2 = s;
		    g = c * e[i];
		    h = c * p;
		    r = hypot2(p,e[i]);
		    e[i+1] = s * r;
		    s = e[i] / r;
		    c = p / r;
		    p = c * d[i] - s * g;
		    d[i+1] = h + s * (c * g + s * d[i]);

		    /* Accumulate transformation. */

		    for (k = 0; k < NDIM; k++)
		    {
			h = V[k][i+1];
			V[k][i+1] = s * V[k][i] + c * h;
			V[k][i] = c * V[k][i] - s * h;
		    }
		}
		
		p = -s * s2 * c3 * el1 * e[l] / dl1;
		e[l] = s * p;
		d[l] = c * p;
		
		/* Check for convergence. */
		
	    } while (fabs(e[l]) > eps*tst1);
	}
	
	d[l] = d[l] + f;
	e[l] = 0.0;
    }
    
    /* Sort eigenvalues and corresponding vectors. */

    for (i = 0; i < NDIM-1; i++)
    {
	double p = d[i];

	k = i;
	
	for (j = i+1; j < NDIM; j++)
	    if (d[j] > p)
	    {
		k = j;
		p = d[j];
	    }
	
	if (k != i)
	{
	    d[k] = d[i];
	    d[i] = p;
	    
	    for (j = 0; j < NDIM; j++)
	    {
		p = V[j][i];
		V[j][i] = V[j][k];
		V[j][k] = p;
	    }
	}
    }
}


void eigen_decomposition(double A[NDIM][NDIM], double d[NDIM], double V[NDIM][NDIM]) 
{
    int i, j;
    double e[NDIM];
    
    for (i = 0; i < NDIM; i++)
	for (j = 0; j < NDIM; j++)
	    V[i][j] = A[i][j];
    
    tred2(V, d, e);
    tql2(V, d, e);
}


int
                lsq_fit(double u[3][3], Transform R)
{
    double           du, omega[6][6], vom[6][6];
    double           dom[6], root2, h[3][3], k[3][3], sign;
    int             i, j, l, rot;

    /* Constant */
    root2 = sqrt(2.0);

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    R[i][j] = ZERO;

    for (i = 0; i < 6; i++)
	for (j = 0; j < 6; j++)
	    omega[i][j] = 0;

    /* Calculate determinant of U */
    du = u[0][0] * u[1][1] * u[2][2] - u[0][0] * u[1][2] * u[2][1]
	- u[0][1] * u[1][0] * u[2][2] + u[0][1] * u[1][2] * u[2][0]
	+ u[0][2] * u[1][0] * u[2][1] - u[0][2] * u[1][1] * u[2][0];

    /* If determinant is zero return */
    if (fabs(du) < SMALL)
	return TRUE;

    /* Make metric matrix omega */
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	{
	    omega[i][j] = ZERO;
	    omega[i + 3][j + 3] = ZERO;
	    omega[i + 3][j] = u[j][i];
	    omega[i][j + 3] = u[i][j];
	}

    /* Diagonalise matrix */
    eigen_decomposition(omega, dom, vom);

    /* Check for degeneracy */
    if (du <= ZERO && fabs(dom[2] - dom[5]) < SMALL)
	return TRUE;

    /* Determine h and k */
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	{
	    h[i][j] = root2 * vom[i][j];
	    k[i][j] = root2 * vom[i + 3][j];
	}

    sign = h[0][0] * h[1][1] * h[2][2] - h[0][0] * h[1][2] * h[2][1]
	- h[0][1] * h[1][0] * h[2][2] + h[0][1] * h[1][2] * h[2][0]
	+ h[0][2] * h[1][0] * h[2][1] - h[0][2] * h[1][1] * h[2][0];

    if (sign <= ZERO)
	for (i = 0; i < 3; i++)
	{
	    h[i][2] = -h[i][2];
	    k[i][2] = -k[i][2];
	}

    /* Determine rotation matrix */
    du /= fabs(du);
    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    R[i][j] = k[j][0] * h[i][0] + k[j][1] * h[i][1] + du * k[j][2] * h[i][2];

    R[0][3] = R[1][3] = R[2][3] = R[3][0] = R[3][1] = R[3][2] = ZERO;
    R[3][3] = 1.0;

    return FALSE;
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
    float           x, y, z;
    char            atmnam[5], chain = '?';
    struct pdbatm  *atom = *atmptr;

    blksize = 10 * MEMGRAIN;
    if (!(atom = malloc(blksize * sizeof(struct pdbatm))))
	fail("read_atoms : Out of memory!");

    while (!feof(ifp))
    {
	if (!fgets(buf, 160, ifp))
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
	
	if (token == ATOM)
	{
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
	    if (strncmp(atmnam, " CA", 3))
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
	}
    }

    *atmptr = atom;
}

int main(int argc, char **argv)
{
    int             i, j, k, ii, jj, kk, l, n, nmax = 0, at1, at2, nequiv, eqlen;
    int             first, last, blksize, hashval, moda, modb, modc, repnum, nclust, modnum;
    int    clustflg[MAXNSTRUC];
    unsigned int    filepos[MAXNSTRUC];
    char            fname[80];
    float           x, y, z, d, r, rsq, rmsd, **tmscmat, cutoff, maxrmsd = 6.0, matchsum[MAXNSTRUC];
    float           maxclusc, tmsc, besttm, d0sq, tmsccut = 0.5, tmsccut2 = 0.6;
    double          u[3][3];
    FILE           *ifp, *ofp, *rfp;
    Point           new, CG_a, CG_b;
    Transform       fr_xf;

    if (argc < 2)
	fail("usage : tmclust ensemble.pdb {tmsccutoff(default=0.65)} {tmsccutoff2(default=0.4)}");

    ifp = fopen(argv[1], "r");	/* Open PDB file in TEXT mode */
    if (!ifp)
    {
	printf("PDB file error!\n");
	exit(-1);
    }

    /* Cutoff for selecting best cluster centroids */
    if (argc > 2)
	tmsccut = atof(argv[2]);

    /* Cutoff for single-linkage clustering around centroid */
    if (argc > 3)
	tmsccut2 = atof(argv[3]);

/*     printf("MAXRMSD = %f\n", maxrmsd); */

    /* Read models */
    for (nmodels=0; nmodels<MAXNSTRUC; nmodels++)
    {
	filepos[nmodels] = ftell(ifp);
	read_atoms(ifp, natoms+nmodels, atoms+nmodels);

	if (!natoms[nmodels])
	    break;

	if (nmodels && natoms[nmodels] != natoms[0])
	{
	    printf("Mismatching number of atoms in model %d (%d vs %d)\n", nmodels+1, natoms[nmodels], natoms[0]);
	    exit(1);
	}
    }

    printf("%d models read from PDB file\n", nmodels);

    tmscmat = allocmat(nmodels, nmodels, sizeof(float), TRUE);

    d0sq = SQR(1.24 * pow(natoms[0]-15.0, 1.0/3.0) - 1.8);
    
    for (moda=0; moda<nmodels; moda++)
    {
	for (modb=moda+1; modb<nmodels; modb++)
	{
	    /* Calculate centroids */
	    veczero(CG_a);
	    veczero(CG_b);
	    for (j = 0; j < natoms[0]; j++)
	    {
		vecadd(CG_a, CG_a, atoms[moda][j].pos);
		vecadd(CG_b, CG_b, atoms[modb][j].pos);
	    }
	    
	    vecscale(CG_a, CG_a, ONE / natoms[0]);
	    vecscale(CG_b, CG_b, ONE / natoms[0]);
	    
	    for (i = 0; i < natoms[0]; i++)
		vecsub(atoms[moda][i].pos, atoms[moda][i].pos, CG_a);
	    for (j = 0; j < natoms[0]; j++)
		vecsub(atoms[modb][j].pos, atoms[modb][j].pos, CG_b);
	    
	    /* Calculate U */
	    for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		    u[i][j] = ZERO;
	    for (j = 0; j < natoms[0]; j++)
		for (k = 0; k < 3; k++)
		    for (l = 0; l < 3; l++)
			u[k][l] += atoms[moda][j].pos[k] * atoms[modb][j].pos[l];
		    
	    if (lsq_fit(u, fr_xf))
	    {
		printf("%d vs %d : LSQFIT failed!\n", moda+1, modb+1);
		tmsc = 0.0;
	    }
	    else
	    {
		for (j = 0; j < natoms[0]; j++)
		{
		    transform_point(fr_xf, atoms[modb][j].pos, new);
		    veccopy(atoms[modb][j].pos, new);
		}
		
		for (tmsc = i = 0; i < natoms[0]; i++)
		{
		    rsq = distsq(atoms[moda][i].pos, atoms[modb][i].pos);
		    tmsc += 1.0F / (1.0F + rsq/d0sq);
		}

		tmsc /= natoms[0];
	    }
	    
	    tmscmat[moda][modb] = tmscmat[modb][moda] = tmsc;
	}
    }
    
    for (moda=0; moda<nmodels; moda++)
	clustflg[moda] = 0;

    /* Eliminate trivial repeat structures */
    for (moda=0; moda<nmodels; moda++)
	if (!clustflg[moda])
	    for (modb=moda+1; modb<nmodels; modb++)
		if (!clustflg[modb] && tmscmat[moda][modb] > 0.999F)
		{
		    printf("Repeated structure removed: %d\n", modb+1);
		    clustflg[modb] = 3;
		}

    modnum = 1;
    rfp = fopen("clustreps.pdb", "w");
    if (!rfp)
	fail("Cannot open clustreps.pdb for writing!");

    for (nclust=1; nclust<=200; nclust++)
    {
	maxclusc = repnum = 0;
	for (moda=0; moda<nmodels; moda++)
	    if (!clustflg[moda])
	    {
		matchsum[moda] = 0;
		for (modb=0; modb<nmodels; modb++)
		    if (modb != moda && !clustflg[modb])
		    {
			if (tmscmat[moda][modb] >= tmsccut)
			    matchsum[moda] += tmscmat[moda][modb];
			if (matchsum[moda] > maxclusc)
			{
			    maxclusc = matchsum[moda];
			    repnum = moda;
			}
		    }
	    }

	if (maxclusc == 0.0)
	{
	    printf("Outliers:\n");
	    for (moda=0; moda<nmodels; moda++)
		if (!clustflg[moda])
		    printf(" %d", moda+1);
	    puts("\n");
	    break;
	}

	sprintf(fname, "CLUSTER_%03d.pdb", nclust);

	printf("Cluster %d:\n", nclust);
	
	ofp = fopen(fname, "w");
	if (!ofp)
	    fail("Cannot open cluster file writing!");

	modc = 0;
	for (moda=0; moda<nmodels; moda++)
	    if (!clustflg[moda] && (moda == repnum || tmscmat[moda][repnum] >= tmsccut))
		clustflg[moda] = 1;
	for (moda=0; moda<nmodels; moda++)
	    if (clustflg[moda] == 1)
		for (modb=0; modb<nmodels; modb++)
		    if (!clustflg[modb] && tmscmat[moda][modb] >= tmsccut2)
			clustflg[modb] = 2;
	
	fprintf(ofp, "REMARK  ORIGINAL REPRESENTATIVE MODEL    %4d\n", repnum+1);
	
	fprintf(rfp, "MODEL    %4d\n", modnum++);
	fprintf(rfp, "REMARK   Original model number: %d\n", repnum+1);
	fseek(ifp, filepos[repnum], SEEK_SET);
	
	while (!feof(ifp))
	{
	    if (!fgets(buf, 160, ifp))
		break;
	    if (!strncmp(buf, "MODEL", 5))
		continue;
	    fputs(buf, rfp);
	    if (!strncmp(buf, "END", 3))
		break;
	}
	
	for (modb=0; modb<nmodels; modb++)
	    if (clustflg[modb] == 1 || clustflg[modb] == 2)
	    {
		printf(" %d", modb+1);
		clustflg[modb] = 3;
		
		fprintf(ofp, "MODEL    %4d\n", ++modc);
		fprintf(ofp, "REMARK   Original model number: %d\n", modb+1);
		fseek(ifp, filepos[modb], SEEK_SET);
		
		while (!feof(ifp))
		{
		    if (!fgets(buf, 160, ifp))
			break;
		    if (!strncmp(buf, "MODEL", 5))
			continue;
		    fputs(buf, ofp);
		    if (!strncmp(buf, "END", 3))
			break;
		}
	    }
	
	fclose(ofp);
	
	printf("\nRepresentative: %d\n\n", repnum+1);
    }
    
    fclose(rfp);
    fclose(ifp);

    return 0;
}
