/* QMODCHECK - evaluate protein model quality using statistical potentials */

/* David Jones, April 2002 */

/* Copyright (C) 2002 University College London */

#define bestDCUTOFF 12.0F

#define DCUTOFF 15.0F

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/times.h>

#define RTCONST (0.582F)

#define noGA
#define MAXSTEPS 100000
#define noELITIST
#define unWEIGHTED
#define MAXFRAGS 5

#define Version		"0.9"
#define Last_Edit	"21st Aug 2004"

/*
 * This program attempts to evaluate the plausibility of a tertiary
 * structure by a sequence shuffling test.
 */

#define FALSE 0
#define TRUE 1

#define BIG (1000000)
#define VBIG (1e32F)

#define MAXSEQLEN 5000

#define MAXNMOT 20

#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define CH alloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);
#define veczero(v) memset(v, 0, sizeof(v))
#define vecadd(a,b,c) (a[0]=b[0]+c[0],a[1]=b[1]+c[1],a[2]=b[2]+c[2])
#define vecsub(a,b,c) (a[0]=b[0]-c[0],a[1]=b[1]-c[1],a[2]=b[2]-c[2])
#define vecscale(a,b,c) ((a[0]=b[0]*(c)),(a[1]=b[1]*(c)),(a[2]=b[2]*(c)))
#define vecprod(a,b,c) ((a[0]=b[1]*c[2]-b[2]*c[1]),(a[1]=b[2]*c[0]-b[0]*c[2]),(a[2]=b[0]*c[1]-b[1]*c[0]))
#define dotprod(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define veccopy(a,b) ((a[0]=b[0]),(a[1]=b[1]),(a[2]=b[2]))
#define distsq(a,b) (SQR(a[0]-b[0])+SQR(a[1]-b[1])+SQR(a[2]-b[2]))
#define dist(a,b) (sqrtf(distsq(a,b)))

/* logistic 'squashing' function (+/- 1.0) */
#define logistic(x) (1.0 / (1.0 + exp(-(x))))

#define TOPOMAX 11
#define INTERVALS 20
#define NACC 5
#define OOICUTOFF 10.0F
#define OOIMIN 6
#define OOIMAX 30
#define OOIDIV 1
#define NRR 40
#define RRWIDTH (1.0F)
#define RTCONST (0.582F)

float     sr_wt = 1.0, lr_wt = 1.0, solvwt = 1.0,
                hbwt = 1.0, compactwt = 0*1.0, sterwt = 1.0, dswt = 0*1.0,
                patwt = 0*1.0, supwt = 0.0;

/* Constants for calculating CB coords (Prot. Eng. Vol.2 p.121) */
#define	TETH_ANG 0.9128
#define CACBDIST 1.538

static char     buf[MAXSEQLEN], seq[MAXSEQLEN], segnum[MAXSEQLEN];
static char     tplt_ss[MAXSEQLEN];
static short    tcbooi[MAXSEQLEN], relacc[MAXSEQLEN], tplt_occ[MAXSEQLEN],
                var_lst[MAXSEQLEN], gaps[MAXSEQLEN], nvart;
static int      firstid, seqlen, nseqs, npats;
static float    cn_energy, econtrib[MAXSEQLEN], e_min = VBIG;
static unsigned short chcount[MAXSEQLEN];

static float cg_x, cg_y, cg_z;

static float    param[10], exprad, expsurf, expmean;

const char     *rnames[] =
{
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "UNK", "UNK"
};

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    GAP, UNK
};

const char     *atmnames[] =
{
    "CA", "CB", "O ", "N ", "C "
};

enum atmcodes
{
    CAATOM, CBATOM, OATOM, NATOM, CATOM
};

enum sstypes
{
  COIL = 1, HELIX = 2, STRAND = 4
};

const char     *sscodes = ".CH.E";

#define NPAIRS 8

enum PAIRS
{
    CA_CA, CB_CB, N_O, O_N, CB_N, N_CB, CB_O, O_CB
};

const int       pairs[NPAIRS][2] =
{
    {CAATOM, CAATOM},
    {CBATOM, CBATOM},
    {CBATOM, NATOM},
    {NATOM, CBATOM},
    {CBATOM, OATOM},
    {OATOM, CBATOM}
};

const float mindist[TOPOMAX][4][4] =
{
  {
    { 2.71, 2.28, 3.14, 1.91 },
    { 2.56, 2.73, 3.15, 2.05 },
    { 1.86, 1.73, 1.65, 1.86 },
    { 2.55, 1.39, 3.16, 2.29 },
  },
  {
    { 4.38, 3.62, 3.38, 3.32 },
    { 4.12, 3.68, 2.66, 3.26 },
    { 3.08, 2.61, 2.30, 2.37 },
    { 4.18, 3.55, 2.71, 3.47 },
  },
  {
    { 3.89, 2.90, 3.04, 3.84 },
    { 2.61, 2.32, 2.26, 3.03 },
    { 2.79, 2.06, 2.53, 2.63 },
    { 3.58, 2.69, 2.73, 3.51 },
  },
  {
    { 3.80, 2.77, 2.68, 3.61 },
    { 2.98, 2.49, 1.97, 3.53 },
    { 2.67, 1.73, 2.67, 2.49 },
    { 3.45, 2.33, 2.63, 3.46 },
  },
  {
    { 3.50, 2.51, 3.06, 3.22 },
    { 2.54, 2.50, 1.69, 2.79 },
    { 2.80, 2.03, 2.76, 2.53 },
    { 3.16, 2.94, 2.61, 3.37 },
  },
  {
    { 3.75, 2.78, 3.07, 3.66 },
    { 2.45, 2.22, 2.02, 2.61 },
    { 2.98, 1.83, 3.07, 2.66 },
    { 3.33, 2.90, 2.79, 3.72 },
  },
  {
    { 3.78, 2.45, 3.16, 3.51 },
    { 2.81, 2.83, 2.57, 2.89 },
    { 2.89, 1.83, 3.08, 2.72 },
    { 3.57, 3.25, 2.69, 3.80 },
  },
  {
    { 3.96, 3.12, 3.00, 3.71 },
    { 3.09, 2.51, 2.90, 2.87 },
    { 3.01, 1.86, 3.09, 2.66 },
    { 3.87, 3.41, 2.70, 3.75 },
  },
  {
    { 3.71, 3.30, 3.01, 3.80 },
    { 3.49, 2.64, 1.88, 3.36 },
    { 3.09, 2.62, 3.01, 2.69 },
    { 3.37, 2.52, 2.68, 3.56 },
  },
  {
    { 3.91, 2.63, 3.10, 3.81 },
    { 2.72, 2.62, 2.13, 2.84 },
    { 2.98, 1.93, 3.17, 2.65 },
    { 3.53, 2.34, 2.74, 3.73 },
  },
  {
    { 3.96, 2.86, 3.05, 3.69 },
    { 2.62, 2.80, 2.11, 2.90 },
    { 3.09, 1.85, 3.05, 2.72 },
    { 3.96, 3.24, 2.70, 3.60 },
  }
};

const float maxdist[TOPOMAX][4][4] =
{
  {
    { 4.11, 5.27, 6.18, 2.86 },
    { 5.11, 6.42, 7.40, 3.92 },
    { 3.68, 4.92, 5.81, 2.57 },
    { 5.25, 6.40, 7.24, 3.97 },
  },
  {
    { 7.71, 8.60, 9.55, 6.37 },
    { 8.63, 9.80, 10.78, 7.38 },
    { 7.06, 8.13, 9.03, 5.84 },
    { 8.87, 9.75, 10.64, 7.46 },
  },
  {
    { 11.00, 12.22, 12.86, 9.80 },
    { 12.15, 13.23, 13.86, 10.91 },
    { 10.46, 11.48, 11.98, 9.16 },
    { 12.10, 13.28, 14.00, 10.89 },
  },
  {
    { 14.52, 15.42, 16.19, 13.15 },
    { 15.52, 16.53, 17.22, 14.21 },
    { 13.65, 14.75, 15.29, 12.52 },
    { 15.60, 16.53, 17.41, 14.33 },
  },
  {
    { 17.75, 18.79, 19.68, 16.50 },
    { 18.96, 20.08, 20.66, 17.63 },
    { 16.98, 18.00, 18.86, 15.73 },
    { 19.00, 20.06, 20.78, 17.68 },
  },
  {
    { 21.15, 22.38, 23.02, 19.93 },
    { 22.09, 23.00, 23.59, 20.83 },
    { 20.18, 21.25, 21.92, 18.91 },
    { 22.32, 23.55, 24.07, 21.09 },
  },
  {
    { 24.59, 25.49, 26.34, 23.27 },
    { 25.54, 26.60, 27.03, 24.28 },
    { 23.59, 24.67, 25.39, 22.42 },
    { 25.69, 26.81, 27.46, 24.39 },
  },
  {
    { 27.98, 28.73, 29.73, 26.64 },
    { 28.76, 29.69, 30.34, 27.52 },
    { 26.89, 27.92, 28.55, 25.51 },
    { 28.97, 29.96, 30.88, 27.75 },
  },
  {
    { 31.29, 32.26, 33.10, 29.99 },
    { 32.35, 33.51, 33.44, 31.08 },
    { 30.22, 31.15, 32.11, 29.00 },
    { 32.38, 33.34, 34.30, 31.19 },
  },
  {
    { 34.39, 35.43, 35.84, 33.20 },
    { 34.98, 35.86, 36.51, 34.06 },
    { 33.60, 34.57, 34.84, 32.36 },
    { 35.76, 36.78, 36.96, 34.55 },
  },
  {
    { 37.75, 38.83, 39.41, 36.52 },
    { 38.74, 39.60, 40.25, 37.35 },
    { 36.54, 37.68, 38.52, 35.32 },
    { 38.50, 39.55, 39.50, 37.26 },
  }
};

const float cbmin[20][20] =
{
    {
	3.99, 3.83, 3.62, 3.76, 3.83, 3.83, 3.65, 2.43, 4.46, 3.85, 3.89, 3.84, 3.72, 4.26, 3.53, 3.61, 3.52, 4.71, 3.61, 3.79 },
    {
	3.83, 4.84, 4.13, 3.78, 4.97, 4.35, 3.86, 2.61, 4.26, 4.28, 3.90, 4.18, 5.51, 3.83, 4.03, 4.19, 4.07, 4.86, 3.79, 4.90 },
    {
	3.62, 4.13, 4.53, 3.77, 4.19, 4.34, 4.18, 2.81, 4.35, 5.29, 3.67, 4.31, 4.43, 4.25, 3.65, 4.04, 4.09, 4.18, 3.85, 4.14 },
    {
	3.76, 3.78, 3.77, 4.02, 3.76, 3.90, 4.54, 2.71, 4.10, 4.36, 3.88, 4.13, 4.36, 3.84, 3.70, 3.63, 3.85, 3.47, 4.13, 4.07 },
    {
	3.83, 4.97, 4.19, 3.76, 3.51, 4.18, 5.02, 2.34, 4.73, 4.20, 4.53, 5.38, 3.94, 4.30, 4.18, 3.93, 3.97, 5.54, 4.00, 4.49 },
    {
	3.83, 4.35, 4.34, 3.90, 4.18, 3.78, 3.68, 2.95, 4.17, 4.52, 4.31, 3.92, 4.20, 4.20, 4.20, 4.09, 3.87, 4.89, 4.28, 4.41 },
    {
	3.65, 3.86, 4.18, 4.54, 5.02, 3.68, 4.00, 3.43, 3.89, 5.04, 4.06, 4.15, 5.49, 4.08, 3.48, 3.87, 4.50, 4.12, 3.86, 4.19 },
    {
	2.43, 2.61, 2.81, 2.71, 2.34, 2.95, 3.43, 1.92, 3.29, 3.61, 3.24, 3.00, 4.17, 3.26, 2.85, 2.53, 2.61, 3.32, 2.87, 3.32 },
    {
	4.46, 4.26, 4.35, 4.10, 4.73, 4.17, 3.89, 3.29, 5.07, 5.32, 3.70, 3.96, 5.23, 4.37, 3.86, 4.07, 4.34, 5.27, 5.03, 4.30 },
    {
	3.85, 4.28, 5.29, 4.36, 4.20, 4.52, 5.04, 3.61, 5.32, 4.04, 3.89, 4.17, 4.11, 4.65, 3.97, 4.06, 4.71, 4.45, 4.74, 4.26 },
    {
	3.89, 3.90, 3.67, 3.88, 4.53, 4.31, 4.06, 3.24, 3.70, 3.89, 4.35, 3.83, 3.85, 4.21, 4.15, 3.59, 4.29, 3.81, 3.92, 4.26 },
    {
	3.84, 4.18, 4.31, 4.13, 5.38, 3.92, 4.15, 3.00, 3.96, 4.17, 3.83, 4.24, 5.20, 4.48, 4.11, 3.62, 4.03, 4.58, 3.75, 4.19 },
    {
	3.72, 5.51, 4.43, 4.36, 3.94, 4.20, 5.49, 4.17, 5.23, 4.11, 3.85, 5.20, 6.40, 4.10, 4.20, 4.53, 3.92, 5.36, 3.99, 4.07 },
    {
	4.26, 3.83, 4.25, 3.84, 4.30, 4.20, 4.08, 3.26, 4.37, 4.65, 4.21, 4.48, 4.10, 4.13, 3.96, 3.87, 4.17, 4.00, 4.68, 4.07 },
    {
	3.53, 4.03, 3.65, 3.70, 4.18, 4.20, 3.48, 2.85, 3.86, 3.97, 4.15, 4.11, 4.20, 3.96, 3.74, 3.85, 3.97, 4.40, 3.94, 4.26 },
    {
	3.61, 4.19, 4.04, 3.63, 3.93, 4.09, 3.87, 2.53, 4.07, 4.06, 3.59, 3.62, 4.53, 3.87, 3.85, 4.21, 3.76, 5.51, 3.61, 3.47 },
    {
	3.52, 4.07, 4.09, 3.85, 3.97, 3.87, 4.50, 2.61, 4.34, 4.71, 4.29, 4.03, 3.92, 4.17, 3.97, 3.76, 4.19, 4.72, 4.32, 4.31 },
    {
	4.71, 4.86, 4.18, 3.47, 5.54, 4.89, 4.12, 3.32, 5.27, 4.45, 3.81, 4.58, 5.36, 4.00, 4.40, 5.51, 4.72, 6.43, 3.74, 4.93 },
    {
	3.61, 3.79, 3.85, 4.13, 4.00, 4.28, 3.86, 2.87, 5.03, 4.74, 3.92, 3.75, 3.99, 4.68, 3.94, 3.61, 4.32, 3.74, 3.89, 4.13 },
    {
	3.79, 4.90, 4.14, 4.07, 4.49, 4.41, 4.19, 3.32, 4.30, 4.26, 4.26, 4.19, 4.07, 4.07, 4.26, 3.47, 4.31, 4.93, 4.13, 4.34 }
};

/* DSSP residue solvent accessible area in GGXGG extended pentapeptide */
const float     resacc[22] =
{
    113.0, 253.0, 167.0, 167.0, 140.0, 199.0, 198.0, 88.0, 194.0, 178.0,
    179.0, 215.0, 194.0, 226.0, 151.0, 134.0, 148.0, 268.0, 242.0, 157.0,
    88.0, 88.0
};

/* Amino acid molecular weights */
const float     molwt[22] =
{
89.09, 174.20, 132.12, 133.10, 121.15, 146.15, 147.13, 75.07, 155.16, 131.17,
    131.17, 146.19, 149.21, 165.19, 115.13, 105.09, 119.12, 204.24, 181.19, 117.15,
    89.09, 89.09
};

/* Amino acid solvation energies */
const float     tfenergy[20] =
{
    0.67, -2.1, -0.6, -1.2, 0.38, -0.22, -0.76, 0.0, 0.64, 1.9,
    1.9, -0.57, 2.4, 2.3, 1.2, 0.01, 0.52, 2.6, 1.6, 1.5
};

/*  BLOSUM 62 */
const short           aamat[23][23] =
{
    {4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3,
     -2, 0, -2, -1, 0},
    {-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3,
     -2, -3, -1, 0, -1},
    {-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4,
     -2, -3, 3, 0, -1},
    {-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4,
     -3, -3, 4, 1, -1},
    {0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2,
     -2, -1, -3, -3, -2},
    {-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2,
     -1, -2, 0, 3, -1},
    {-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3,
     -2, -2, 1, 4, -1},
    {0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2,
     -3, -3, -1, -2, -1},
    {-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2,
     2, -3, 0, 0, -1},
    {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3,
     -1, 3, -3, -3, -1},
    {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2,
     -1, 1, -4, -3, -1},
    {-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3,
     -2, -2, 0, 1, -1},
    {-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1,
     -1, 1, -3, -1, -1},
    {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1,
     3, -1, -3, -3, -1},
    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4,
     -3, -2, -2, -1, -2},
    {1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3,
     -2, -2, 0, 0, 0},
    {0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2,
     -2, 0, -1, -1, 0},
    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11,
     2, -3, -4, -3, -2},
    {-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2,
     7, -1, -3, -2, -1},
    {0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3,
     -1, 4, -3, -2, -1},
    {-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4,
     -3, -3, 4, 1, -1},
    {-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3,
     -2, -2, 1, 4, -1},
    {0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2,
     -1, -1, -1, -1, 4}
};

typedef float   Transform[4][4];
typedef float   Vector[3];

typedef struct
{
    Vector           n, h, ca, c, o, cb;
    float           phi, psi, omega;
    short           fcode;
}
RESATM;

static struct chnentry
{
    short           length, *relacc, **psimat;
    char           *seq, *sstruc;
    RESATM         *chain;
} chnlist[1000];

struct chnidx
{
    short           chain, pos, length;
    unsigned int    scount;
    float           weight, rmsd;
};

static struct flist
{
    struct chnidx  *frags;
    int             nfrags, ntot;
    short           nchn, from, bestlen;
    float           totwt, totbest, rmsd, totsum, totsumsq;
    char           *type;
}
fraglist[MAXSEQLEN + 1];

static struct SSTRUC
{
    short           start, length, type;
}
sstlist[100];

static int      nchn, chnres, ncontact;

static float radgyr, contord;

static float ***cb_pmat, **cb_targ, **cb_last, **ooi_pmat, **acc_pmat;

static FILE    *rfp;

typedef struct
{
    RESATM         *genome;
    float           perfval, selval;
    short           evalflg;
}
Schema;

static float    diffdist[TOPOMAX][4][4];

const char     *rescodes = "ARNDCQEGHILKMFPSTWYVX";

/* Statistical Sampling Arrays */

/* Local conformation tables */
static float    sr_de[4][4][21][21][INTERVALS][TOPOMAX];

/* MR Residue-Residue interaction matrices */
static float    mr_de[4][4][21][21][NRR];

/* LR Residue-Residue interaction matrices */
static float    lr_de[4][4][21][21][NRR];

/* Residue accessibility matrices */
static float    acc_de[21][NACC];

/* Residue OOI number matrices */
static float    ooi_de[21][(OOIMAX - OOIMIN) / OOIDIV + 1];

/* Structure generation arrays */
static float    (**distmat)[NPAIRS], **hbdmat;
static char   **hbflg;
static RESATM  *bestchn;


/* Dump a rude message to standard error and exit */
void
                fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

/* Return accessibility range */
int
                accrange(int acc)
{
    if (acc < 12)
	return 0;
    else if (acc < 36)
	return 1;
    else if (acc < 44)
	return 2;
    else if (acc < 87)
	return 3;
    return 4;
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

void           *allocmat(int rows, int columns, int size)
{
    int             i;
    void          **p;

    p = malloc(rows * sizeof(void *));

    if (p == NULL)
	fail("allocmat: malloc [] failed!");
    for (i = 0; i < rows; i++)
	if ((p[i] = calloc(columns, size)) == NULL)
	    fail("allocmat: malloc [][] failed!");

    return p;
}

void
                freemat(void *p, int rows)
{
    int             i;

    for (i = 0; i < rows; i++)
	free(((void **) p)[i]);
    free(p);
}

/* Return potential term for at1(a) -> at2(b) */
float pairpot(const int a, const int b, const int at1, const int at2, const int t, const float d)
{
    int    k;

    if (t > TOPOMAX)
    {
	k = (d-4.0F) / RRWIDTH;
	if (k < 0)
	    k = 0;
	else if (k >= NRR)
	    k = NRR - 1;

	if (t >= 23)
	    return lr_de[at1][at2][a][b][k];
	else
	    return mr_de[at1][at2][a][b][k];
    }

    if (t >= 1)
    {
	k = (int) (INTERVALS * (d - mindist[t-1][at1][at2]) / diffdist[t-1][at1][at2]);
	if (k < 0)
	    k = 0;
	if (k < INTERVALS)
	    return sr_de[at1][at2][a][b][k][t-1];
    }

    return 0.0F;
}

/* Secondary structure hydrogen bond evaluation */
float           calchb(float (**distmat)[NPAIRS], RESATM * chn)
{
    int             i, j, ndonate[MAXSEQLEN], naccept[MAXSEQLEN], donor[MAXSEQLEN][2], acceptor[MAXSEQLEN], nhbs=0, bestj;
    float           energy[MAXSEQLEN], hb_energy = 0.0F, costheta, e, ebest, chi;
    float           dnosq, dho;
    Vector          c_nvec, ca_nvec, n_hvec;

    /* Build backbone amide hydrogens (Pauling & Corey rule) */
    for (i = 0; i < seqlen; i++)
    {
	ndonate[i] = naccept[i] = 0;
	acceptor[i] = -999;
	donor[i][0] = donor[i][1] = -999;
	
	/* Generate main chain amide H position bisecting C-N-CA (mean of C->N and CA->N directions) */
	if (i)
	{
	    vecsub(c_nvec, chn[i].n, chn[i-1].c);
	    vecscale(c_nvec, c_nvec, 1.0F / dist(chn[i-1].c, chn[i].n));
	    vecsub(ca_nvec, chn[i].n, chn[i].ca);
	    vecscale(ca_nvec, ca_nvec, 1.0F / dist(chn[i].n, chn[i].ca));
	    
	    vecadd(n_hvec, c_nvec, ca_nvec);

	    vecscale(n_hvec, n_hvec, 1.0F / sqrtf(dotprod(n_hvec, n_hvec)));

	    vecadd(chn[i].h, chn[i].n, n_hvec);
	}
    }
    
    for (i = 0; i < seqlen; i++)
	if (seq[i] != PRO)
	{
	    /* Find best acceptor for this donor */
	    bestj = -1;
	    ebest = VBIG;
	    for (j = 0; j < seqlen; j++)
		if (abs(j-i) > 2 && distmat[i][j][CA_CA] < 9.0F && naccept[j] < 2)
		{
		    dnosq = distsq(chn[i].n, chn[j].o);

		    if (dnosq > 25.0F)
			continue;

		    dho = dist(chn[i].h, chn[j].o);

		    if (i)
		    {
			/* Dot product N->H . O->H (N.B. N-H vector is already unit length) */
			costheta = ((chn[i].h[0] - chn[i].n[0]) * (chn[i].h[0] - chn[j].o[0]) + (chn[i].h[1] - chn[i].n[1]) * (chn[i].h[1] - chn[j].o[1]) + (chn[i].h[2] - chn[i].n[2]) * (chn[i].h[2] - chn[j].o[2])) / dho;

			/* DREIDING hydrogen bond potential with cos^2 weighting */
			e = 5.5F * (5.0F * powf(8.41F / dnosq, 6.0F) - 6.0F * powf(8.41F / dnosq, 5.0F)) * SQR(costheta);
		    }
		    else
		    {
			/* Non-directional potential for first N-H */
			e = 5.5F * (5.0F * powf(8.41F / dnosq, 6.0F) - 6.0F * powf(8.41F / dnosq, 5.0F));
		    }

		    if (e < ebest)
		    {
			ebest = e;
			bestj = j;
		    }
		}
		    
	    if (ebest < -0.5F)
	    {
		ndonate[i]++;
		acceptor[i] = bestj;
		energy[i] = ebest;
		donor[bestj][naccept[bestj]++] = i;
	    }
	}

    for (i = 0; i < seqlen; i++)
    {
	/* Parallel bridges */
	if (acceptor[i] >= 0)
	    if (acceptor[i] == donor[i][0]-2 || acceptor[i] == donor[i][1]-2)
		hb_energy += energy[i];
	    else if (acceptor[i] == donor[i][0]+2 || acceptor[i] == donor[i][1]+2)
		hb_energy += energy[i];
	
	/* Antiparallel bridges */
	if (acceptor[i] >= 0 && (acceptor[i] == donor[i][0] || acceptor[i] == donor[i][1]))
	    hb_energy += energy[i];
	else if (i > 1 && acceptor[i] >= 0 && (acceptor[i] == donor[i-2][0]-2 || acceptor[i] == donor[i-2][1]-2))
	    hb_energy += energy[i];
    }

    /* printf("NHBs = %d\n", nhbs); */
    
    return hb_energy;
}

/* Old structure hydrogen bond evaluation */
float           calchb2(float (**distmat)[NPAIRS], RESATM * chn)
{
    int             i, j, ndonate[MAXSEQLEN], naccept[MAXSEQLEN], donor[MAXSEQLEN][2], acceptor[MAXSEQLEN], nhbs=0, bestj;
    float           energy[MAXSEQLEN], hb_energy = 0.0F, costheta, e, ebest, chi;
    float           dnosq, dho;
    Vector          c_nvec, ca_nvec, n_hvec;

    /* Build backbone amide hydrogens (Pauling & Corey rule) */
    for (i = 0; i < seqlen; i++)
    {
	ndonate[i] = naccept[i] = 0;
	acceptor[i] = -999;
	donor[i][0] = donor[i][1] = -999;
	
	/* Generate main chain amide H position bisecting C-N-CA (mean of C->N and CA->N directions) */
	if (i)
	{
	    vecsub(c_nvec, chn[i].n, chn[i-1].c);
	    vecscale(c_nvec, c_nvec, 1.0F / dist(chn[i-1].c, chn[i].n));
	    vecsub(ca_nvec, chn[i].n, chn[i].ca);
	    vecscale(ca_nvec, ca_nvec, 1.0F / dist(chn[i].n, chn[i].ca));
	    
	    vecadd(n_hvec, c_nvec, ca_nvec);

	    vecscale(n_hvec, n_hvec, 1.0F / sqrtf(dotprod(n_hvec, n_hvec)));

	    vecadd(chn[i].h, chn[i].n, n_hvec);
	}
    }
    
    for (i = 0; i < seqlen; i++)
	if (seq[i] != PRO)
	{
	    /* Find best acceptor for this donor */
	    bestj = -1;
	    ebest = VBIG;
	    for (j = 0; j < seqlen; j++)
		if (abs(j-i) > 2 && distmat[i][j][CA_CA] < 9.0F && naccept[j] < 2)
		{
		    dnosq = distsq(chn[i].n, chn[j].o);

		    if (dnosq > 25.0F)
			continue;

		    dho = dist(chn[i].h, chn[j].o);

		    if (i)
		    {
			/* Dot product N->H . O->H (N.B. N-H vector is already unit length) */
			costheta = ((chn[i].h[0] - chn[i].n[0]) * (chn[i].h[0] - chn[j].o[0]) + (chn[i].h[1] - chn[i].n[1]) * (chn[i].h[1] - chn[j].o[1]) + (chn[i].h[2] - chn[i].n[2]) * (chn[i].h[2] - chn[j].o[2])) / dho;

			/* DREIDING hydrogen bond potential with cos^2 weighting */
			e = 5.5F * (5.0F * powf(8.41F / dnosq, 6.0F) - 6.0F * powf(8.41F / dnosq, 5.0F)) * SQR(costheta);
		    }
		    else
		    {
			/* Non-directional potential for first N-H */
			e = 5.5F * (5.0F * powf(8.41F / dnosq, 6.0F) - 6.0F * powf(8.41F / dnosq, 5.0F));
		    }

		    if (e < ebest)
		    {
			ebest = e;
			bestj = j;
		    }
		}
		    
	    if (ebest < -0.5F)
	    {
		ndonate[i]++;
		acceptor[i] = bestj;
		energy[i] = ebest;
		donor[bestj][naccept[bestj]++] = i;
		if (bestj == i-4 || abs(i - bestj) > 9)
		    hb_energy += ebest;
	    }
	}
    
    return hb_energy;
}


int calcdisu(RESATM *chn, float (**distmat)[NPAIRS])
{
    int             i, j, ndisulph;

    for (ndisulph = i = 0; i < seqlen; i++)
	if (seq[i] == CYS)
	    for (j = 0; j < seqlen; j++)
		if (abs(j-i) > 2 && seq[j] == CYS && distmat[i][j][CA_CA] < 6.5F && distmat[i][j][CB_CB] < 4.5)
		    ndisulph++;
    
    return ndisulph;
}


float           torsion(Vector a, Vector b, Vector c, Vector d)
{
    int             i;
    float           a_b[3], b_c[3], c_d[3], len_a_b, len_b_c, len_c_d;
    float           ab_bc, ab_cd, bc_cd, s_ab, s_bc, costor, costsq, tor;
    float           sintor, conv = 0.01745329, sign;

    /* Calculate vectors and lengths a-b, b-c, c-d */
    for (i = 0; i < 3; i++)
    {
	a_b[i] = b[i] - a[i];
	b_c[i] = c[i] - b[i];
	c_d[i] = d[i] - c[i];
    }

    len_a_b = sqrt(dotprod(a_b, a_b));
    len_b_c = sqrt(dotprod(b_c, b_c));
    len_c_d = sqrt(dotprod(c_d, c_d));

    /* Error check, are any vectors of zero length ? */
    if ((len_a_b == 0.0F) || (len_b_c == 0.0F) || (len_c_d == 0.0F))
	return (-999.0F);

    /* Calculate dot products to form cosines */
    ab_bc = dotprod(a_b, b_c) / (len_a_b * len_b_c);
    ab_cd = dotprod(a_b, c_d) / (len_a_b * len_c_d);
    bc_cd = dotprod(b_c, c_d) / (len_b_c * len_c_d);

    /* Form sines */
    s_ab = sqrt(1.0F - SQR(ab_bc));
    s_bc = sqrt(1.0F - SQR(bc_cd));

    costor = (ab_bc * bc_cd - ab_cd) / (s_ab * s_bc);
    costsq = SQR(costor);
    if ((costsq >= 1.0F) && (costor < 0.0F))
	tor = 180.0F;

    /* If the angle is not == 180 degs calculate sign using sine */
    if (costsq < 1.0F)
    {
	sintor = sqrt(1.0F - costsq);
	tor = atan2(sintor, costor);
	tor /= conv;

	/* Find unit vectors */
	for (i = 0; i < 3; i++)
	{
	    a_b[i] /= len_a_b;
	    b_c[i] /= len_b_c;
	    c_d[i] /= len_c_d;
	}

	/* Find determinant */
	sign = a_b[0] * (b_c[1] * c_d[2] - b_c[2] * c_d[1]) + a_b[1] * (b_c[2] * c_d[0] - b_c[0] * c_d[2]) + a_b[2] * (b_c[0] * c_d[1] - b_c[1] * c_d[0]);

	/* Change sign if necessary */
	if (sign < 0.0F)
	    tor = -tor;
    }

    return (tor);
}

void calctor(RESATM *chn, int nres)
{
  int i;

  for (i = 0; i < nres; i++)
    {
      if (i < nres-1)
	{
	  chn[i].psi = torsion(chn[i].n, chn[i].ca, chn[i].c, chn[i+1].n);
	  chn[i].omega = torsion(chn[i].ca, chn[i].c, chn[i+1].n, chn[i+1].ca);
	}
      else
	chn[i].psi = chn[i].omega = 999.999F;
      if (i)
	chn[i].phi = torsion(chn[i-1].c, chn[i].n, chn[i].ca, chn[i].c);
      else
	chn[i].phi = 999.999F;
      /*      printf("%3d %f %f %f\n", i, chn[i].phi, chn[i].psi, chn[i].omega); */
    }
}


unsigned long rng_x=123456789, rng_y=362436069, rng_z=521288629, rng_w=88675123;

/* Fast Marsaglia XOR-shift RNG with period 2^128-1 */
unsigned long xor128(void)
{
    unsigned long t;

    t = (rng_x^(rng_x<<11));
    rng_x = rng_y;
    rng_y = rng_z;
    rng_z = rng_w;

    return rng_w = (rng_w^(rng_w>>19))^(t^(t>>8));
}

/* Generate random number 0<=x<1 */
#define ran0()  (xor128()*(1.0/4294967296.0))

/* randint(a,b) : return random integer a <= n <= b */
#define randint(low,high) ((int) ((low) + ((high)-(low)+1) * ran0()))

void
                randomise(void)
{
    char dum;

    /* Attempt to generate a random state unique to this process/host/time */
    rng_x = (unsigned long)time(NULL);
    rng_y = (unsigned long)&dum;
    rng_z = (unsigned long)gethostid();
    rng_w = (unsigned long)getpid();
}

/* Perform random shuffle of array */
void            shuffle(int *s, int len)
{
    int             i, r, temp;

    for (i = len - 1; i >= 1; i--)
    {
	r = ran0() * (i + 1);
	temp = s[i];
	s[i] = s[r];
	s[r] = temp;
    }
}

/* Byte-swap for little-endian systems */
void byteswap(void *f, register int n)
{
    register unsigned int *vp;

    vp = (unsigned int *)f;

    while (n--)
    {
        *vp = ((*vp >> 24) & 0xFF) | ((*vp >> 8) & 0x0000FF00) | ((*vp << 8) & 0x00FF0000) | ((*vp << 24) & 0xFF000000);
	vp++;
    }
}

void
                readxyz(char *buf, float *x, float *y, float *z)
{
    char            temp[9];

    temp[8] = '\0';
    strncpy(temp, buf, 8);
    *x = atof(temp);
    strncpy(temp, buf + 8, 8);
    *y = atof(temp);
    strncpy(temp, buf + 16, 8);
    *z = atof(temp);
}

/* Read target structure coords */
void
                maketemplate(char *pdbname, char chainid, int model, RESATM *chn)
{
    FILE           *pfp, *dfp, *xfp;
    char            buf[160], whichatm[MAXSEQLEN], inscode;
    int             acc, atc, i, j, k, nres, resid, lastid, nallat, atmresn[MAXSEQLEN*20], nseg;
    float           f, dv, cca[3], nca[3], xx[3], yy[3], sx, sy, r;
    short           at1, at2;
    float           x[MAXSEQLEN][5], y[MAXSEQLEN][5], z[MAXSEQLEN][5];
    float           allats_x[MAXSEQLEN*20], allats_y[MAXSEQLEN*20], allats_z[MAXSEQLEN*20];

    printf("Reading %s chain %c ...\n", pdbname, chainid);
    pfp = fopen(pdbname, "r");
    if (!pfp)
	fail("maketemplate: Cannot open PDB file!");

    for (i = 0; i < MAXSEQLEN; i++)
	whichatm[i] = 0;

    if (model >= 1)
    {
	while (!feof(pfp))
	{
	    if (fgets(buf, 160, pfp) == NULL)
		fail("Cannot find MODEL!");
	    if (strncmp(buf, "MODEL", 5))
		continue;
	    sscanf(buf + 5, "%d", &i);
	    if (i == model)
		break;
	}
    }

    i = -1;
    nallat = nseg = 0;
    while (!feof(pfp))
    {
	if (fgets(buf, 160, pfp) == NULL)
	    break;
	if (!strncmp(buf, "END", 3))
	    break;
	if (!strncmp(buf, "TER", 3))
	    nseg++;
	/* printf("%d %s\n",i,buf); */
	if (strncmp(buf, "ATOM", 4) || buf[21] != chainid || buf[16] != ' ' && buf[16] != 'A')
	    continue;
	for (atc = 0; atc <= CATOM; atc++)
	    if (!strncmp(buf + 13, atmnames[atc], 2))
	    {
		inscode = buf[26];
		buf[26] = ' ';
		sscanf(buf + 22, "%d", &resid);
		if (atc == NATOM)
		{
		    ++i;
		    readxyz(buf + 30, &x[i][atc], &y[i][atc], &z[i][atc]);
		    if (!i)
			firstid = resid;
		    else if (inscode == ' ' && resid - lastid > 1)
		    {
			if (SQR(x[i][NATOM] - x[i - 1][NATOM]) + SQR(y[i][NATOM] - y[i - 1][NATOM]) + SQR(z[i][NATOM] - z[i - 1][NATOM]) > 15.0)
			{
			    printf("WARNING: Skipping %d missing residues!\n", resid - lastid - 1);
			    for (k = 0; k < resid - lastid - 1; k++)
			    {
				seq[i++] = GAP;
			    }
			}
		    }
		    for (k = 0; k < 20; k++)
			if (!strncmp(rnames[k], buf + 17, 3))
			    break;
		    seq[i] = k;
		    lastid = resid;
		    segnum[i] = nseg;
		}
		else if (i >= 0)
		    readxyz(buf + 30, &x[i][atc], &y[i][atc], &z[i][atc]);
		whichatm[i] |= 1 << atc;
	    }
	
	if (buf[13] != 'H' && buf[13] != 'D')
	{
	    readxyz(buf + 30, &allats_x[nallat], &allats_y[nallat], &allats_z[nallat]);
	    atmresn[nallat] = i+1;
	    nallat++;
	}
    }
    fclose(pfp);

    nres = i + 1;

    /* Check atoms */
    for (i = 0; i < nres; i++)
	if (seq[i] != GAP)
	{
	    if (!(whichatm[i] & (1 << NATOM)))
	    {
		printf("FATAL: Missing N atom in %d!\n", i + 1);
		exit(1);
	    }
	    if (!(whichatm[i] & (1 << CAATOM)))
	    {
		printf("WARNING: Missing CA atom in %d!\n", i + 1);
	    }
	    if (!(whichatm[i] & (1 << CBATOM)))
	    {
		if (!(whichatm[i] & (1 << CAATOM)) || !(whichatm[i] & (1 << CATOM)) || !(whichatm[i] & (1 << NATOM)))
		{
		    /* Not much left of this residue! */
		    printf("WARNING: Missing main-chain atom in %d!\n", i + 1);
		    seq[i] = GAP;
		    continue;
		}
		
		/* Reconstruct CB atom */
		nca[0] = x[i][CAATOM] - x[i][NATOM];
		nca[1] = y[i][CAATOM] - y[i][NATOM];
		nca[2] = z[i][CAATOM] - z[i][NATOM];
		cca[0] = x[i][CAATOM] - x[i][CATOM];
		cca[1] = y[i][CAATOM] - y[i][CATOM];
		cca[2] = z[i][CAATOM] - z[i][CATOM];
		for (k = 0; k < 3; k++)
		    xx[k] = nca[k] + cca[k];
		vecprod(yy, nca, cca);
		sx = CACBDIST * cos(TETH_ANG) / sqrt(dotprod(xx, xx));
		sy = CACBDIST * sin(TETH_ANG) / sqrt(dotprod(yy, yy));
		x[i][CBATOM] = x[i][CAATOM] + xx[0] * sx + yy[0] * sy;
		y[i][CBATOM] = y[i][CAATOM] + xx[1] * sx + yy[1] * sy;
		z[i][CBATOM] = z[i][CAATOM] + xx[2] * sx + yy[2] * sy;
		whichatm[i] |= 1 << CBATOM;
	    }
	}
    
    for (i = 0; i < nres; i++)
      {
	chn[i].n[0] = x[i][NATOM];
	chn[i].n[1] = y[i][NATOM];
	chn[i].n[2] = z[i][NATOM];
	chn[i].ca[0] = x[i][CAATOM];
	chn[i].ca[1] = y[i][CAATOM];
	chn[i].ca[2] = z[i][CAATOM];
	chn[i].c[0] = x[i][CATOM];
	chn[i].c[1] = y[i][CATOM];
	chn[i].c[2] = z[i][CATOM];
	chn[i].o[0] = x[i][OATOM];
	chn[i].o[1] = y[i][OATOM];
	chn[i].o[2] = z[i][OATOM];
	chn[i].cb[0] = x[i][CBATOM];
	chn[i].cb[1] = y[i][CBATOM];
	chn[i].cb[2] = z[i][CBATOM];
      }

    /* Calculate radius of gyration */
    for (i = 0; i < nallat; i++)
    {
	cg_x += allats_x[i];
	cg_y += allats_y[i];
	cg_z += allats_z[i];
    }

    cg_x /= nallat;
    cg_y /= nallat;
    cg_z /= nallat;

    for (radgyr = i = 0; i < nallat; i++)
    {
	radgyr += SQR(cg_x - allats_x[i]) + SQR(cg_y - allats_y[i]) + SQR(cg_z - allats_z[i]);
    }

    radgyr = sqrt(radgyr / nallat);

    /* Calculate contact order */
    for (i = 0; i < nallat; i++)
	for (j = i+1; j < nallat; j++)
	    if (abs(atmresn[j] - atmresn[i]) > 0 && SQR(allats_x[i] - allats_x[j]) + SQR(allats_y[i] - allats_y[j]) + SQR(allats_z[i] - allats_z[j]) < 36.0F)
	    {
		contord += abs(atmresn[j] - atmresn[i]);
		ncontact++;
	    }

    contord /= ncontact * nres;
    
    /* Calculate interatomic distance templates */
    at1 = at2 = CBATOM;
    for (i = 0; i < nres; i++)
	for (j = i + 1; j < nres; j++)
	{
	    dv = sqrt(SQR(x[i][at1] - x[j][at2]) + SQR(y[i][at1] - y[j][at2]) + SQR(z[i][at1] - z[j][at2]));
	    distmat[i][j][CB_CB] = distmat[j][i][CB_CB] = dv;
	}

    /* Calculate interatomic distance templates */
    at1 = at2 = CAATOM;
    for (i = 0; i < nres; i++)
	for (j = i + 1; j < nres; j++)
	{
	    dv = sqrt(SQR(x[i][at1] - x[j][at2]) + SQR(y[i][at1] - y[j][at2]) + SQR(z[i][at1] - z[j][at2]));
	    distmat[i][j][CA_CA] = distmat[j][i][CA_CA] = dv;
	}

    for (i = 0; i < nres; i++)
	tcbooi[i] = 0;

    /* Compute 10 A CB-Ooi numbers */
    for (i = 0; i < nres; i++)
	for (j = i + 1; j < nres; j++)
	    if (distmat[i][j][CB_CB] <= OOICUTOFF && segnum[i] == segnum[j])
	    {
		tcbooi[i]++;
		tcbooi[j]++;
	    }

    seqlen = nres;
}

int main(int argc, char **argv)
{
    int             a, aa, b, i, j, k, l, n, t, at1, at2, seqidx[MAXSEQLEN], aafreq[20];
    float           p, pair, solv, natpair, natsolv, pairsum, pairsumsq, solvsum, solvsumsq, zpair, zsolv, dcut, ehb, min_ncd, rpair[MAXSEQLEN], rsolv[MAXSEQLEN], rnatpair[MAXSEQLEN], rnatsolv[MAXSEQLEN], scsum[MAXSEQLEN], aaprob[20], avpair, avsolv;
    double sum, sumsq, sdpair, sdsolv, av, sd, rpsum[MAXSEQLEN], rpsumsq[MAXSEQLEN], rssum[MAXSEQLEN], rssumsq[MAXSEQLEN];
    FILE           *ifp;
    char            chainid, id[512];

    printf("Model Checker\n");
    printf("Build date : %s\n",__DATE__);
    printf("Copyright (C) 1994/2000/2002 David T. Jones\n\n");

    if (argc < 2)
    {
	fprintf(stderr, "usage : %s pdbfname [chainid]\n", argv[0]);
	exit(1);
    }

    randomise();

    puts("\nReading mean force potential tables...");
/*    ifp = fopen("/data/threader/potentials.dat", "rb"); */
    ifp = fopen("modcheckpot.dat", "rb");
    if (!ifp)
	fail("Cannot open modcheckpot.dat!");
    fread(sr_de, 1, sizeof(sr_de), ifp);
    fread(mr_de, 1, sizeof(mr_de), ifp);
    fread(lr_de, 1, sizeof(lr_de), ifp);
    fread(acc_de, 1, sizeof(acc_de), ifp);
    fread(ooi_de, 1, sizeof(ooi_de), ifp);
    fclose(ifp);

#if defined(sun) || defined(sgi)
    byteswap(sr_de, sizeof(sr_de)/4);
    byteswap(mr_de, sizeof(mr_de)/4);
    byteswap(lr_de, sizeof(lr_de)/4);
    byteswap(acc_de, sizeof(acc_de)/4);
    byteswap(ooi_de, sizeof(ooi_de)/4);
#endif

    for (i = 0; i < TOPOMAX; i++)
	for (a = 0; a < 4; a++)
	    for (b = 0; b < 4; b++)
		diffdist[i][a][b] = maxdist[i][a][b] - mindist[i][a][b];

    distmat = allocmat(MAXSEQLEN, MAXSEQLEN, sizeof(**distmat));

    bestchn = malloc(MAXSEQLEN * sizeof(RESATM));

    if (argc == 2)
    {
	maketemplate(argv[1], ' ', -1, bestchn);
    }
    else if (argc == 3)
    {
	if (strlen(argv[2]) != 1)
	    sscanf(argv[2] + 1, "%d", &a);
	else
	    a = -1;
	maketemplate(argv[1], argv[2][0], a, bestchn);
    }

    if (seqlen < 20)
    {
	printf(" 999\n");
	return 0;
    }

    for (i=0; i<seqlen; i++)
	seqidx[i] = i;

    for (a=0; a<20; a++)
	aafreq[a] = 0;
    
    pair = solv = 0.0;
    for (i=0; i<seqlen; i++)
    {
	if (seq[i] < 20)
	    aafreq[seq[i]]++;
	rnatpair[i] = rnatsolv[i] = 0.0;
    }

    for (a=0; a<20; a++)
	aaprob[a] = (float) aafreq[a] / seqlen;
    
    for (i=0; i<seqlen; i++)
	if (seq[i] < 20)
	{
	    for (j=i+1; j<seqlen; j++)
		if (seq[j] < 20)
		{
		    if (distmat[i][j][CB_CB] < DCUTOFF && segnum[i] == segnum[j])
		    {
			rnatpair[i] += (p = pairpot(seq[i], seq[j], CBATOM, CBATOM, j-i, distmat[i][j][CB_CB]));
			rnatpair[j] += p;
			pair += p;
		    }
		}
	    rnatsolv[i] += (p = ooi_de[seq[i]][MAX(0, MIN(OOIMAX, tcbooi[i]) - OOIMIN) / OOIDIV]);
	    solv += p;
	}
    
    natpair = pair;
    natsolv = solv;

    avpair = avsolv = sdpair = sdsolv = 0.0;
    for (i=0; i<seqlen; i++)
    {
	for (j=i+1; j<seqlen; j++)
	{
	    if (distmat[i][j][CB_CB] < DCUTOFF && segnum[i] == segnum[j])
	    {
		for (sum=a=0; a<20; a++)
		    for (b=0; b<20; b++)
		    {
			avpair += (p = pairpot(a, b, CBATOM, CBATOM, j-i, distmat[i][j][CB_CB])) * aaprob[a] * aaprob[b];
			sum += p * aaprob[a] * aaprob[b];
		    }
		for (a=0; a<20; a++)
		    for (b=0; b<20; b++)
		    {
			p = pairpot(a, b, CBATOM, CBATOM, j-i, distmat[i][j][CB_CB]) - sum;
			sdpair += SQR(p) * aaprob[a] * aaprob[b];
		    }
	    }
	}
	
	for (sum=a=0; a<20; a++)
	    avsolv += ooi_de[a][MAX(0, MIN(OOIMAX, tcbooi[i]) - OOIMIN) / OOIDIV] * aaprob[a];
	for (sum=a=0; a<20; a++)
	    sdsolv += SQR(ooi_de[a][MAX(0, MIN(OOIMAX, tcbooi[i]) - OOIMIN) / OOIDIV] - sum) * aaprob[a];
    }

    printf("avpair = %f avsolv = %f\n", avpair, avsolv);

    sdpair = sqrt(sdpair) * 1.8;
    sdsolv = sqrt(sdsolv);

    printf("sdpair = %f sdsolv = %f\n", sdpair, sdsolv);
    
    min_ncd = 99999.0;
    for (i=0; i<5; i++)
	for (j=0; j<5; j++)
	    if (distmat[i][seqlen-j-1][CB_CB] < min_ncd)
		min_ncd = distmat[i][seqlen-j-1][CB_CB];
    
    zpair = (natpair - avpair) / sdpair;
    
    zsolv = (natsolv - avsolv) / sdsolv;

    printf("Epair = %f (Zpair = %f)\n", -natpair * RTCONST, zpair);
    printf("Esolv = %f (Zsolv = %f)\n", -natsolv * RTCONST, zsolv);
    
    printf("Hydrogen bond score = %f\n", ehb = calchb(distmat, bestchn));

    printf("MQAP Score %f\n", logistic(0.2*(zsolv + zpair - 10.0)));

    printf(" %f %f %f %f %f %f %d %f %f %f %d\n", -natpair * RTCONST, zpair, -natsolv * RTCONST, zsolv, zpair+zsolv, ehb, calcdisu(bestchn, distmat), radgyr, contord, min_ncd, seqlen);

//    printf(" %f\n", -(zpair+zsolv));
//    printf(" %f\n", -zpair);
    printf(" %f\n", -(0.29*zpair+0.17*zsolv)-0.2*contord);
//    printf(" %f\n", ehb);

    return 0;
}
