/* Protein Chain Clustering Program by David Jones, May 2004 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXNSTRUC 50000

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


#define ROTATE(a,i,j,k,l) (g=a[i][j],h=a[k][l],a[i][j]=g-s*(h+g*tau),a[k][l]=h+s*(g-h*tau))

int             jacobi(float a[6][6], float d[6], float v[6][6], int *nrot)
{
    int             j, iq, ip, i;
    float           tresh, theta, tau, t, sm, s, h, g, c, b[6], z[6];

    for (ip = 0; ip < 6; ip++)
    {
	for (iq = 0; iq < 6; iq++)
	    v[ip][iq] = ZERO;
	v[ip][ip] = 1.0;
    }

    for (ip = 0; ip < 6; ip++)
    {
	b[ip] = d[ip] = a[ip][ip];
	z[ip] = ZERO;
    }

    *nrot = 0;
    for (i = 0; i < 50; i++)
    {
	sm = ZERO;
	for (ip = 0; ip < 5; ip++)
	{
	    for (iq = ip + 1; iq < 6; iq++)
		sm += fabs(a[ip][iq]);
	}
	if (sm == ZERO)
	    return 0;
	if (i < 3)
	    tresh = 0.2 * sm / 36;
	else
	    tresh = ZERO;
	for (ip = 0; ip < 5; ip++)
	{
	    for (iq = ip + 1; iq < 6; iq++)
	    {
		g = 100.0F * fabs(a[ip][iq]);
		if (i > 3 && fabs(d[ip]) + g == fabs(d[ip])
		    && fabs(d[iq]) + g == fabs(d[iq]))
		    a[ip][iq] = ZERO;
		else
		    if (fabs(a[ip][iq]) > tresh)
		    {
			h = d[iq] - d[ip];
			if (fabs(h) + g == fabs(h))
			    t = (a[ip][iq]) / h;
			else
			{
			    theta = 0.5 * h / (a[ip][iq]);
			    t = 1.0F / (fabs(theta) + sqrtf(1.0F + theta * theta));
			    if (theta < ZERO)
				t = -t;
			}
			c = 1.0F / sqrtf(1.0F + t * t);
			s = t * c;
			tau = s / (1.0F + c);
			h = t * a[ip][iq];
			z[ip] -= h;
			z[iq] += h;
			d[ip] -= h;
			d[iq] += h;
			a[ip][iq] = ZERO;
			for (j = 0; j <= ip - 1; j++)
			{
			    ROTATE(a, j, ip, j, iq);
			}
			for (j = ip + 1; j <= iq - 1; j++)
			{
			    ROTATE(a, ip, j, j, iq);
			}
			for (j = iq + 1; j < 6; j++)
			{
			    ROTATE(a, ip, j, iq, j);
			}
			for (j = 0; j < 6; j++)
			{
			    ROTATE(v, j, ip, j, iq);
			}
			++(*nrot);
		    }
	    }
	}
	for (ip = 0; ip < 6; ip++)
	{
	    b[ip] += z[ip];
	    d[ip] = b[ip];
	    z[ip] = ZERO;
	}
    }
    return 1;
}

#undef ROTATE

/*
 * function eigsrt 
 *
 * Given the eigenvalues d[n] and eignevectors v[n][n] as output from jacobi
 * this routine sourts the eigenvalues into descending order and rearranges
 * the columns of v correspondingly 
 */
void            eigsrt(float d[6], float v[6][6])
{
    int             k, j, i;
    float           p;

    for (i = 0; i < 5; i++)
    {
	p = d[k = i];
	for (j = i + 1; j < 6; j++)
	    if (d[j] >= p)
		p = d[k = j];
	if (k != i)
	{
	    d[k] = d[i];
	    d[i] = p;
	    for (j = 0; j < 6; j++)
	    {
		p = v[j][i];
		v[j][i] = v[j][k];
		v[j][k] = p;
	    }
	}
    }
}

int
                lsq_fit(float u[3][3], Transform R)
{
    float du;
    float           omega[6][6], vom[6][6];
    float           dom[6], h[3][3], k[3][3], sign;
    const float root2 = sqrt(2.0);
    int             i, j, l, rot;

    for (i = 0; i < 3; i++)
	for (j = 0; j < 3; j++)
	    R[i][j] = ZERO;

    for (i = 0; i < 6; i++)
	for (j = 0; j < 6; j++)
	    omega[i][j] = 0;

    /* Calculate determinant of U */
    du = (double) u[0][0] * u[1][1] * u[2][2] - u[0][0] * u[1][2] * u[2][1]
	- u[0][1] * u[1][0] * u[2][2] + u[0][1] * u[1][2] * u[2][0]
	+ u[0][2] * u[1][0] * u[2][1] - u[0][2] * u[1][1] * u[2][0];

    /* If determinant is zero return */
    if (fabs(du) < 1.0e-5)
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
    if (jacobi(omega, dom, vom, &rot))
	return TRUE;

    /* Sort by eigenvalues */
    eigsrt(dom, vom);

    /* Check for degeneracy */
    if (du <= ZERO && fabs(dom[2] - dom[5]) < 1.0e-5)
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
    R[3][3] = ONE;

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
    float           x, y, z, u[3][3];
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

main(int argc, char **argv)
{
    int             i, j, k, ii, jj, kk, l, n, nmax = 0, at1, at2, nequiv, eqlen;
    int             first, last, blksize, hashval, moda, modb, modc, repnum, nclust, modnum;
    int    clustflg[MAXNSTRUC];
    unsigned int    filepos[MAXNSTRUC];
    char            fname[80];
    float           x, y, z, d, r, rsq, rmsd, u[3][3], **tmscmat, cutoff, maxrmsd = 6.0, matchsum[MAXNSTRUC];
    float           maxclusc, tmsc, besttm, d0sq, tmsccut = 0.5, tmsccut2 = 0.6;
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
		if (!clustflg[modb] && tmscmat[moda][modb] > 0.99F)
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
