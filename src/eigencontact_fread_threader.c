/****************************************************************************
 *
 * EiGenTHREADER version 0.4 (Sep 2015)
 * written by David T. Jones,
 * Department of Computer Science,
 * University College London,
 * Gower Street,
 * London
 *
 * Email: d.t.jones@ucl.ac.uk
 *
 ****************************************************************************/

/*
 * This program attempts to align a protein contact map with structures in a fold library
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/times.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>

#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

#define BIG 1000000000
#define VBIG 1e32F

/* "Tweakable" Parameters */

/* Scaling factor for fixed-point arithmetic */
#define SCALEFAC 1000

/* Max-size-of-problem parameters */
#define MAXSEQLEN 5000
#define MAXATOMS 10000

#define SQR(x) ((x)*(x))
#define ABS(x) ((x)>=0?(x):-(x))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define CH malloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);
#define vecprod(a,b,c) (a[0]=b[1]*c[2]-b[2]*c[1],a[1]=b[2]*c[0]-b[0]*c[2],a[2]=b[0]*c[1]-b[1]*c[0])
#define dotprod(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])


/* Constants for calculating CB coords (Prot. Eng. Vol.2 p.121) */
#define	TETH_ANG 0.9128
#define CACBDIST 1.538

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

enum atmcodes
{
    CAATOM, CBATOM, OATOM, NATOM, CATOM
};

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    GAP, UNK
};

const float     ZERO = 0.0F, ONE = 1.0F, TWO = 2.0F, THREE = 3.0F;

char            brkid[10], seq1[MAXSEQLEN], seq2[10000];
char            modseq[MAXSEQLEN], tsstruc[MAXSEQLEN], coreflg[MAXSEQLEN];
char            ssstruc[10000], ssstrel[10000], *resid[MAXSEQLEN];
short           modsdx[MAXSEQLEN], modsdx2[MAXSEQLEN], gaps[MAXSEQLEN];
short           sstidx[MAXSEQLEN], tooi[MAXSEQLEN];
int             maxgplen[MAXSEQLEN];
float           targf1[MAXSEQLEN][20], targf2[MAXSEQLEN][20];
int             tpltsc1[MAXSEQLEN][20], tpltsc2[MAXSEQLEN][20], firstid,
                seq1len, seq2len, nseqs;
float           rmax;

float          **eigenvec1, *eigenval1, **eigenvec2, *eigenval2;

int             lc1, lc2;

int           **pmat;

int             psidata;

struct hash
{
    short           resid, acc;
    char            inscode, sstruc;
}
hashtbl[MAXSEQLEN];

struct SSTRUC
{
    short           start, length, type;
}
sstlist[100];

float           e_contrib[MAXSEQLEN], e_cav[MAXSEQLEN], contrib_min, contrib_max;

struct HIT
{
    char            *brkid;
    float           contactsc, nwsc, perco_a, perco_b;
} hits[100000];

typedef struct
{
    short           length;
    short           posn_a[MAXSEQLEN], posn_b[MAXSEQLEN];
}
ALNS;

int             sstlen;
FILE           *efp;

const char     *rnames[] =
{
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "GAP", "UNK"
};

const char     *ssnames[] =
{
    "COIL", "HELIX", "STRAND", "TURN"
};

enum sscodes
{
    COIL, HELIX, STRAND, TURN
};

const char     *atmnames[] =
{
    "CA", "CB", "O ", "N ", "C "
};

#define NPAIRS 5

const short     atompair[][2] =
{
    {CBATOM, CBATOM},
    {CBATOM, NATOM},
    {CBATOM, OATOM},
    {NATOM, CBATOM},
    {OATOM, CBATOM},
};

const char     *rescodes = "ARNDCQEGHILKMFPSTWYV-X";
const char     *sscodes = "CHE";

/* Amino acid composition of SWISS-PROT 47.0 */
float dbaaf[20] =
{
    0.078558, 0.053457, 0.041873, 0.053116, 0.015533, 0.039327, 0.066123, 0.069598, 0.022845, 0.059217,
    0.096416, 0.059123, 0.023859, 0.040095, 0.048484, 0.068456, 0.054442, 0.011580, 0.030628, 0.067271
};

short           mutcode[22] =
{
    16, 11, 3, 6, 19, 8, 3, 0, 18, 12, 12, 1, 19, 18, 15, 16, 15, 18, 13, 9, 20, 20
};

short ssmat[3][3] = {
    {  135, -46, -472 },
    { -200,  90, -374 },
    { -162, -398, 215 }
};

/*  BLOSUM 50 */
short           aamat2[23][23] =
{
    {500, -200, -100, -200, -100, -100, -100, 000, -200, -100, -200, -100, -100, -300, -100, 100, 000, -300, -200, 000, -200, -100, -100},
    {-200, 700, -100, -200, -400, 100, 000, -300, 000, -400, -300, 300, -200, -300, -300, -100, -100, -300, -100, -300, -100, 000, -100},
    {-100, -100, 700, 200, -200, 000, 000, 000, 100, -300, -400, 000, -200, -400, -200, 100, 000, -400, -200, -300, 400, 000, -100},
    {-200, -200, 200, 800, -400, 000, 200, -100, -100, -400, -400, -100, -400, -500, -100, 000, -100, -500, -300, -400, 500, 100, -100},
    {-100, -400, -200, -400, 1300, -300, -300, -300, -300, -200, -200, -300, -200, -200, -400, -100, -100, -500, -300, -100, -300, -300, -200},
    {-100, 100, 000, 000, -300, 700, 200, -200, 100, -300, -200, 200, 000, -400, -100, 000, -100, -100, -100, -300, 000, 400, -100},
    {-100, 000, 000, 200, -300, 200, 600, -300, 000, -400, -300, 100, -200, -300, -100, -100, -100, -300, -200, -300, 100, 500, -100},
    {000, -300, 000, -100, -300, -200, -300, 800, -200, -400, -400, -200, -300, -400, -200, 000, -200, -300, -300, -400, -100, -200, -200},
    {-200, 000, 100, -100, -300, 100, 000, -200, 1000, -400, -300, 000, -100, -100, -200, -100, -200, -300, 200, -400, 000, 000, -100},
    {-100, -400, -300, -400, -200, -300, -400, -400, -400, 500, 200, -300, 200, 000, -300, -300, -100, -300, -100, 400, -400, -300, -100},
    {-200, -300, -400, -400, -200, -200, -300, -400, -300, 200, 500, -300, 300, 100, -400, -300, -100, -200, -100, 100, -400, -300, -100},
    {-100, 300, 000, -100, -300, 200, 100, -200, 000, -300, -300, 600, -200, -400, -100, 000, -100, -300, -200, -300, 000, 100, -100},
    {-100, -200, -200, -400, -200, 000, -200, -300, -100, 200, 300, -200, 700, 000, -300, -200, -100, -100, 000, 100, -300, -100, -100},
    {-300, -300, -400, -500, -200, -400, -300, -400, -100, 000, 100, -400, 000, 800, -400, -300, -200, 100, 400, -100, -400, -400, -200},
    {-100, -300, -200, -100, -400, -100, -100, -200, -200, -300, -400, -100, -300, -400, 1000, -100, -100, -400, -300, -300, -200, -100, -200},
    {100, -100, 100, 000, -100, 000, -100, 000, -100, -300, -300, 000, -200, -300, -100, 500, 200, -400, -200, -200, 000, 000, -100},
    {000, -100, 000, -100, -100, -100, -100, -200, -200, -100, -100, -100, -100, -200, -100, 200, 500, -300, -200, 000, 000, -100, 000},
    {-300, -300, -400, -500, -500, -100, -300, -300, -300, -300, -200, -300, -100, 100, -400, -400, -300, 1500, 200, -300, -500, -200, -300},
    {-200, -100, -200, -300, -300, -100, -200, -300, 200, -100, -100, -200, 000, 400, -300, -200, -200, 200, 800, -100, -300, -200, -100},
    {000, -300, -300, -400, -100, -300, -300, -400, -400, 400, 100, -300, 100, -100, -300, -200, 000, -300, -100, 500, -400, -300, -100},
    {-200, -100, 400, 500, -300, 000, 100, -100, 000, -400, -400, 000, -300, -400, -200, 000, 000, -500, -300, -400, 500, 200, -100},
    {-100, 000, 000, 100, -300, 400, 500, -200, 000, -300, -300, 100, -100, -400, -100, 000, -100, -200, -200, -300, 200, 500, -100},
    {-100, -100, -100, -100, -200, -100, -100, -200, -100, -100, -100, -100, -100, -200, -200, -100, 000, -300,
     -100, -100, -100, -100, -100}
};

/* BLOSUM 62 */
short           aamat[23][23] =
{
    {400, -100, -200, -200, 000, -100, -100, 000, -200, -100, -100, -100, -100, -200, -100, 100, 000, -300, -200, 000, -200, -100, -100},
    {-100, 500, 000, -200, -300, 100, 000, -200, 000, -300, -200, 200, -100, -300, -200, -100, -100, -300, -200, -300, -100, 000, -100},
    {-200, 000, 600, 100, -300, 000, 000, 000, 100, -300, -300, 000, -200, -300, -200, 100, 000, -400, -200, -300, 300, 000, -100},
    {-200, -200, 100, 600, -300, 000, 200, -100, -100, -300, -400, -100, -300, -300, -100, 000, -100, -400, -300, -300, 400, 100, -100},
    {000, -300, -300, -300,1000, -300, -400, -300, -300, -100, -100, -300, -100, -200, -300, -100, -100, -200, -200, -100, -300, -300, -100},
    {-100, 100, 000, 000, -300, 500, 200, -200, 000, -300, -200, 100, 000, -300, -100, 000, -100, -200, -100, -200, 000, 300, -100},
    {-100, 000, 000, 200, -400, 200, 500, -200, 000, -300, -300, 100, -200, -300, -100, 000, -100, -300, -200, -200, 100, 400, -100},
    {000, -200, 000, -100, -300, -200, -200, 600, -200, -400, -400, -200, -300, -300, -200, 000, -200, -200, -300, -300, -100, -200, -100},
    {-200, 000, 100, -100, -300, 000, 000, -200, 800, -300, -300, -100, -200, -100, -200, -100, -200, -200, 200, -300, 000, 000, -100},
    {-100, -300, -300, -300, -100, -300, -300, -400, -300, 400, 200, -300, 100, 000, -300, -200, -100, -300, -100, 300, -300, -300, -100},
    {-100, -200, -300, -400, -100, -200, -300, -400, -300, 200, 400, -200, 200, 000, -300, -200, -100, -200, -100, 100, -400, -300, -100},
    {-100, 200, 000, -100, -300, 100, 100, -200, -100, -300, -200, 500, -100, -300, -100, 000, -100, -300, -200, -200, 000, 100, -100},
    {-100, -100, -200, -300, -100, 000, -200, -300, -200, 100, 200, -100, 500, 000, -200, -100, -100, -100, -100, 100, -300, -100, -100},
    {-200, -300, -300, -300, -200, -300, -300, -300, -100, 000, 000, -300, 000, 600, -400, -200, -200, 100, 300, -100, -300, -300, -100},
    {-100, -200, -200, -100, -300, -100, -100, -200, -200, -300, -300, -100, -200, -400, 700, -100, -100, -400, -300, -200, -200, -100, -100},
    {100, -100, 100, 000, -100, 000, 000, 000, -100, -200, -200, 000, -100, -200, -100, 400, 100, -300, -200, -200, 000, 000, -100},
    {000, -100, 000, -100, -100, -100, -100, -200, -200, -100, -100, -100, -100, -200, -100, 100, 500, -200, -200, 000, -100, -100, -100},
    {-300, -300, -400, -400, -200, -200, -300, -200, -200, -300, -200, -300, -100, 100, -400, -300, -200, 1100, 200, -300, -400, -300, -100},
    {-200, -200, -200, -300, -200, -100, -200, -300, 200, -100, -100, -200, -100, 300, -300, -200, -200, 200, 700, -100, -300, -200, -100},
    {000, -300, -300, -300, -100, -200, -200, -300, -300, 300, 100, -200, 100, -100, -200, -200, 000, -300, -100, 400, -300, -200, -100},
    {-200, -100, 300, 400, -300, 000, 100, -100, 000, -300, -400, 000, -300, -300, -200, 000, -100, -400, -300, -300, 400, 100, -100},
    {-100, 000, 000, 100, -300, 300, 400, -200, 000, -300, -300, 100, -100, -300, -100, 000, -100, -300, -200, -200, 100, 400, -100},
    {-100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, 000}
};

/* Structure generation arrays */
float           (**distmat)[4][4];
float           x[MAXSEQLEN][5], y[MAXSEQLEN][5], z[MAXSEQLEN][5];
float           phi[MAXSEQLEN], psi[MAXSEQLEN], omega[MAXSEQLEN];
short           trelacc[MAXSEQLEN];
float           tmat[MAXSEQLEN][22];

/* Vars for Gotoh-N-W alignment routine */
int           **pat;

int             nsst;

float **contmap;

/* Control parameters */
int             SHUFCOUNT = 0;
int             MINSSTREL = 30;
float           SC_PARAM[20];	/* TEST parameters */
int             GP_OPEN[3] = { 748, 1263, 711 };
int             GP_EXT[3] = { 64, 162, 284 };
float           SCUTOFF = 20.0;
int             NEIGEN = 7;
float           SCALE = 1250.0F;
float           DCUTOFF = 12.0F;
int             ssscore = 100;
int             prtalnflg = FALSE, verbflg = FALSE, motifflg = FALSE, filtgapflg = FALSE;
int             mod3flg = FALSE;
int             htmlflg = FALSE;
int             seqssflg = FALSE;
int             quickalnflg = FALSE;
int             ooisolv = FALSE;
int             ssstflg = FALSE;
char            *ssfname, *htmname;
char            motpref[160];

/* Dump a rude message to standard error and exit */
void
                fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

#define freevec(x) (free(x))

/* Allocate vector */
void           *allocvec(int columns, int size, int clrflg)
{
    void          *p;

    if (clrflg)
	p = calloc(columns, size);
    else
	p = malloc(columns * size);

    if (p == NULL)
	fail("allocvec: calloc failed!");

    return p;
}

void           *allocmat(int rows, int columns, int size, int clrflg)
{
    int             i;
    void          **p;

    p = malloc(rows * sizeof(void *));

    if (p == NULL)
	fail("allocmat: malloc [] failed!");
    if (clrflg)
    {
	if ((p[0] = calloc(columns * rows, size)) == NULL)
	    fail("allocmat: calloc [][] failed!");
    }
    else
	if ((p[0] = malloc(columns * rows * size)) == NULL)
	    fail("allocmat: malloc [][] failed!");

    for (i=1; i<rows; i++)
	p[i] = (char *)p[i-1] + columns * size;

    return p;
}

void            freemat(void *p)
{
    free(*(void **)p);
    free(p);
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

/* Convert string to lower case */
void
                lc_str(char *str)
{
    int             i;

    for (i = strlen(str) - 1; i >= 0; i--)
	if (isupper(str[i]))
	    str[i] = tolower(str[i]);
}


unsigned int rng_x=123456789, rng_y=362436069, rng_z=521288629, rng_w=88675123;

/* Fast Marsaglia XOR-shift RNG with period 2^128-1 */
unsigned int xor128(void)
{
    unsigned int t;

    t = (rng_x^(rng_x<<11));
    rng_x = rng_y;
    rng_y = rng_z;
    rng_z = rng_w;

    return rng_w = (rng_w^(rng_w>>19))^(t^(t>>8));
}

/* Generate random number 0<=x<1 */
#define ran0()  (xor128()*(1.0/4294967296.0))

/* randint(a,b) : return random integer a <= n <= b */
#define randint(low,high) ((low) + (int)(((high)-(low)+1) * ran0()))


/* Attempt to generate a random state unique to this process/host/time */
void
                randomise(void)
{
    int i;

    rng_x = (unsigned int)time(NULL);
    rng_y = (unsigned int)getppid();
    rng_z = (unsigned int)gethostid();
    rng_w = (unsigned int)getpid();

    if (verbflg)
	printf("Random seeds: %u %u %u %u\n", rng_x, rng_y, rng_z, rng_w);

    /* Warm up the generator */
    for (i=0; i<100; i++)
	xor128();
}


/* Perform random shuffle of sequence skipping GAPS/UNKs */
void            shufseqng(char *s, int len)
{
    int             i, r;
    char            temp;

    for (i = len - 1; i >= 1; i--)
	if (s[i] < 20)
	{
	    r = ran0() * (i + 1);
	    if (s[r] < 20)
	    {
		temp = s[i];
		s[i] = s[r];
		s[r] = temp;
	    }
	}
}

/* Perform random shuffle of sequence */
void            shufseq(char *s, int len)
{
    int             i, r;
    char            temp;

    for (i = len - 1; i >= 1; i--)
    {
	r = ran0() * (i + 1);
	temp = s[i];
	s[i] = s[r];
	s[r] = temp;
    }
}

/* Perform random rotation of sequence */
void            rotseq(char *s, int len)
{
    int             i, n;
    char            temp;

    n = randint(1, len - 1);

    while (n--)
    {
	temp = s[0];
	for (i = 0; i < len - 1; i++)
	    s[i] = s[i + 1];
	s[len - 1] = temp;
    }
}

/* Print alignment */
void
                prtalign(ALNS * newtem)
{
    int             i, b, id = 0, nb, p1, p2, seqflg;
    char            r1, r2, sim;

    nb = newtem->length / 60 + 1;

    for (b = 0; b < nb; b++)
    {
	printf("         ");
	for (i = 0; i < 58; i++)
	{
	    if (b * 60 + i + 3 > newtem->length)
		break;
	    p1 = newtem->posn_a[b * 60 + i + 3];
	    if (p1 && !(p1 % 10))
	    {
		printf("%3d", p1);
		i += 2;
	    }
	    else
		printf(" ");
	}
	putchar('\n');
	for (seqflg = 0; seqflg < 2; seqflg++)
	{
	    if (seqflg == 1)
		printf("%-8.8s ", brkid);
	    else
		printf("         ");
	    for (i = 0; i < 60; i++)
	    {
		if (b * 60 + i >= newtem->length)
		    break;
		p1 = newtem->posn_a[b * 60 + i + 1];

		switch (seqflg)
		{
		case 0:
		    r1 = (p1) ? sscodes[tsstruc[p1 - 1]] : '-';
		    break;
		case 1:
		    r1 = (p1) ? rescodes[seq1[p1 - 1]] : '-';
		    break;
		}

		putchar(r1);
	    }
	    putchar('\n');
	}
	printf("         ");
	for (i = 0; i < 60; i++)
	{
	    if (b * 60 + i >= newtem->length)
		break;
	    p1 = newtem->posn_a[b * 60 + i + 1];
	    p2 = newtem->posn_b[b * 60 + i + 1];
	    r1 = (p1) ? rescodes[seq1[p1 - 1]] : '-';
	    r2 = (p2) ? rescodes[seq2[p2 - 1]] : '-';
	    if (r1 == r2 && r1 != '-')
	    {
		id++;
		sim = '|';
	    }
	    else
		sim = ' ';
	    putchar(sim);
	}
	putchar('\n');
	printf("%-8.8s ", motpref[0] ? motpref : "Query");
	for (i = 0; i < 60; i++)
	{
	    if (b * 60 + i >= newtem->length)
		break;
	    p2 = newtem->posn_b[b * 60 + i + 1];
	    r2 = (p2) ? rescodes[seq2[p2 - 1]] : '-';
	    putchar(r2);
	}
	putchar('\n');
	if (ssstflg && seqssflg)
	{
	    printf("         ");
	    for (i = 0; i < 60; i++)
	    {
		if (b * 60 + i >= newtem->length)
		    break;
		p2 = newtem->posn_b[b * 60 + i + 1];
		r2 = (p2) ? sscodes[ssstruc[p2 - 1]] : '-';
		putchar(r2);
	    }
	    putchar('\n');
	}
	printf("         ");
	for (i = 0; i < 58; i++)
	{
	    if (b * 60 + i + 3 > newtem->length)
		break;
	    p2 = newtem->posn_b[b * 60 + i + 3];
	    if (p2 && !(p2 % 10))
	    {
		printf("%3d", p2);
		i += 2;
	    }
	    else
		printf(" ");
	}
	puts("\n\n");
    }
    printf("Percentage Identity = %3.1f.\n\n", 100.0F * id / MIN(seq1len, seq2len));
}

/* Print alignment in HTML format */
void
                htmlalign(ALNS * newtem)
{
    int             i, b, id = 0, nb, p1, p2, seqflg, fontcolour;
    char            r1, r2;
    FILE *hfp;

    hfp = fopen(htmname, "a");
    if (!hfp)
	fail("Cannot append to HTML output file!");

    nb = newtem->length / 60 + 1;

    fprintf(hfp, "<hr><body text=\"#ffffff\" bgcolor=\"#000000\"><tt><pre><br><a name=\"%s\"></a>", brkid);

    for (b = 0; b < nb; b++)
    {
	fprintf(hfp, "         ");
	for (i = 0; i < 58; i++)
	{
	    if (b * 60 + i + 3 > newtem->length)
		break;
	    p1 = newtem->posn_a[b * 60 + i + 3];
	    if (p1 && !(p1 % 10))
	    {
		fprintf(hfp, "<font color=\"#ffffff\">%3d</font>", p1);
		i += 2;
	    }
	    else
		fprintf(hfp, " ");
	}
	fputc('\n', hfp);
	for (seqflg = 0; seqflg < 2; seqflg++)
	{
	    if (seqflg == 1)
		fprintf(hfp, "<font color=\"#ffffff\">%-8.8s</font> ", brkid);
	    else
		fprintf(hfp, "         ");
	    fontcolour = 0;
	    for (i = 0; i < 60; i++)
	    {
		if (b * 60 + i >= newtem->length)
		    break;
		p1 = newtem->posn_a[b * 60 + i + 1];
		p2 = newtem->posn_b[b * 60 + i + 1];

		switch (seqflg)
		{
		case 0:
		    r1 = (p1) ? sscodes[tsstruc[p1 - 1]] : '-';
		    switch(r1)
		    {
		    case 'H':
			if (fontcolour != 1)
			{
			    if (fontcolour)
				fprintf(hfp, "</font>");
			    fprintf(hfp, "<font color=\"#ff00ff\">");
			    fontcolour = 1;
			}
			break;
		    case 'E':
			if (fontcolour != 2)
			{
			    if (fontcolour)
				fprintf(hfp, "</font>");
			    fprintf(hfp, "<font color=\"#ffff00\">");
			    fontcolour = 2;
			}
			break;
		    default:
			if (fontcolour != 3)
			{
			    if (fontcolour)
				fprintf(hfp, "</font>");
			    fprintf(hfp, "<font color=\"#ffffff\">");
			    fontcolour = 3;
			}
			break;
		    }
		    fprintf(hfp, "%c", r1);
		    break;
		case 1:
		    r1 = (p1) ? rescodes[seq1[p1 - 1]] : '-';
		    r2 = (p2) ? rescodes[seq2[p2 - 1]] : '-';
		    if (r1 == r2 && r1 != '-')
		    {
			if (fontcolour != 4)
			{
			    if (fontcolour)
				fprintf(hfp, "</font>");
			    fprintf(hfp, "<font color=\"#00ff00\">");
			    fontcolour = 4;
			}
			fprintf(hfp, "%c", r1);
		    }
		    else if (p1 && p2 && aamat[seq1[p1 - 1]][seq2[p2 - 1]] > 100)
		    {
			if (fontcolour != 5)
			{
			    if (fontcolour)
				fprintf(hfp, "</font>");
			    fprintf(hfp, "<font color=\"#ff0000\">");
			    fontcolour = 5;
			}
			fprintf(hfp, "%c", r1);
		    }
		    else
		    {
			if (fontcolour != 6)
			{
			    if (fontcolour)
				fprintf(hfp, "</font>");
			    fprintf(hfp, "<font color=\"#6060ff\">");
			    fontcolour = 6;
			}
			fprintf(hfp, "%c", r1);
		    }
		    break;
		}
	    }
	    fprintf(hfp, "</font>\n");
	}
	fprintf(hfp, "<font color=\"#ffffff\">%-8.8s</font> ", motpref[0] ? motpref : "Query");

	fontcolour = 0;

	for (i = 0; i < 60; i++)
	{
	    if (b * 60 + i >= newtem->length)
		break;
	    p1 = newtem->posn_a[b * 60 + i + 1];
	    p2 = newtem->posn_b[b * 60 + i + 1];
	    r1 = (p1) ? rescodes[seq1[p1 - 1]] : '-';
	    r2 = (p2) ? rescodes[seq2[p2 - 1]] : '-';
	    if (r1 == r2 && r1 != '-')
	    {
		id++;
		if (fontcolour != 4)
		{
		    if (fontcolour)
			fprintf(hfp, "</font>");
		    fprintf(hfp, "<font color=\"#00ff00\">");
		    fontcolour = 4;
		}
	    }
	    else if (p1 && p2 && aamat[seq1[p1 - 1]][seq2[p2 - 1]] > 100)
	    {
		if (fontcolour != 5)
		{
		    if (fontcolour)
			fprintf(hfp, "</font>");
		    fprintf(hfp, "<font color=\"#ff0000\">");
		    fontcolour = 5;
		}
	    }
	    else
	    {
		if (fontcolour != 6)
		{
		    if (fontcolour)
			fprintf(hfp, "</font>");
		    fprintf(hfp, "<font color=\"#6060ff\">");
		    fontcolour = 6;
		}
	    }
	    fprintf(hfp, "%c", r2);
	}
	fprintf(hfp, "</font>\n");

	fontcolour = 0;

	if (ssstflg && seqssflg)
	{
	    fprintf(hfp, "         ");
	    for (i = 0; i < 60; i++)
	    {
		if (b * 60 + i >= newtem->length)
		    break;
		p2 = newtem->posn_b[b * 60 + i + 1];
		r2 = (p2) ? sscodes[ssstruc[p2 - 1]] : '-';
		switch(r2)
		{
		case 'H':
		    if (fontcolour != 1)
		    {
			if (fontcolour)
			    fprintf(hfp, "</font>");
			fprintf(hfp, "<font color=\"#ff00ff\">");
			fontcolour = 1;
		    }
		    break;
		case 'E':
		    if (fontcolour != 2)
		    {
			if (fontcolour)
			    fprintf(hfp, "</font>");
			fprintf(hfp, "<font color=\"#ffff00\">");
			fontcolour = 2;
		    }
		    break;
		default:
		    if (fontcolour != 3)
		    {
			if (fontcolour)
			    fprintf(hfp, "</font>");
			fprintf(hfp, "<font color=\"#ffffff\">");
			fontcolour = 3;
		    }
		    break;
		}
		fprintf(hfp, "%c", r2);
	    }
	}
	fprintf(hfp, "</font>\n");
	fprintf(hfp, "         ");
	for (i = 0; i < 58; i++)
	{
	    if (b * 60 + i + 3 > newtem->length)
		break;
	    p2 = newtem->posn_b[b * 60 + i + 3];
	    if (p2 && !(p2 % 10))
	    {
		fprintf(hfp, "<font color=\"#ffffff\">%3d</font>", p2);
		i += 2;
	    }
	    else
		fprintf(hfp, " ");
	}
	fprintf(hfp, "\n\n");
    }
    fprintf(hfp, "</pre><font color=\"#ffffff\">Percentage Identity = %3.1f%%</font></body><P>\n", 100.0F * id / MIN(seq1len, seq2len));

    fclose(hfp);
}

/* Calculate point-biserial score */
float pointbiserial()
{
    int i, j, n=0, n0=0, n1=0;
    float m0=0.0, m1=0.0, delta, mean=0.0, variance=0.0, sd, rpb;

    for (i = 0; i < seq2len; i++)
	for (j = i+2; j < seq2len; j++)
	    if (modsdx2[i] >= 0 && modsdx2[j] >= 0)
	    {
		n++;

		if (distmat[modsdx2[i]][modsdx2[j]][CBATOM][CBATOM] < DCUTOFF)
		{
		    n1++;
		    m1 += contmap[i][j];
		}
		else
		{
		    n0++;
		    m0 += contmap[i][j];
		}
		delta = contmap[i][j] - mean;
		mean += delta / n;
		variance += delta * (contmap[i][j] - mean);
	    }

    if (!n0 || !n1 || n < 2)
	return 0.0F;

    m0 /= n0;
    m1 /= n1;
    variance /= n;

    sd = sqrtf(variance);

    rpb = (m1 - m0) * sqrtf((float)n1 * n0 / SQR((float)n)) / sd;

    return rpb * sqrtf((n1 + n0 - 2.0F)/(1.0F - SQR(rpb)));
}

/* Calculate Jones point-biserial score */
float jonespointbiserial()
{
    int i, j, n;
    float frac, n01, n0=0.0, n1=0.0, m0=0.0, m1=0.0, delta, mean=0.0, variance=0.0, sd, rpb, viol;

    for (n = i = 0; i < seq2len; i++)
	for (j = i+2; j < seq2len; j++)
	    if (modsdx2[i] >= 0 && modsdx2[j] >= 0)
	    {
		viol = distmat[modsdx2[i]][modsdx2[j]][CBATOM][CBATOM] - DCUTOFF;

		frac = 0.5F * expf(-viol*viol);

		if (viol < 0.0F)
		{
		    n0 += frac;
		    m0 += frac * contmap[i][j];
		    n1 += 1.0F - frac;
		    m1 += (1.0F - frac) * contmap[i][j];
		}
		else
		{
		    n1 += frac;
		    m1 += frac * contmap[i][j];
		    n0 += 1.0F - frac;
		    m0 += (1.0F - frac) * contmap[i][j];
		}

		n++;
		delta = contmap[i][j] - mean;
		mean += delta / n;
		variance += delta * (contmap[i][j] - mean);
	    }

    n01 = n0 + n1;

    if (n01 < 6.0F)
	return 0.0F;

    m0 /= n0;
    m1 /= n1;
    variance /= n;

    sd = sqrtf(variance);

    rpb = (m1 - m0) * sqrtf(n1 * n0 / SQR(n01)) / sd;

    return rpb * sqrtf((n01 - 2.0F)/(1.0F - SQR(rpb)));
}

/* Calculate correl t-score */
float correlation()
{
    int i, j, n=0;
    float pc, dviol, delta, mean=0.0, variance=0.0, sd, rpb, *x, *y;
    float correl, xav=0.0, yav=0.0, sxy=0.0, sxx=0.0, syy=0.0;

    x = allocvec(seq2len * (seq2len - 1) / 2, sizeof(float), FALSE);
    y = allocvec(seq2len * (seq2len - 1) / 2, sizeof(float), FALSE);

    for (i = 0; i < seq2len; i++)
	for (j = i+2; j < seq2len; j++)
	    if (modsdx2[i] >= 0 && modsdx2[j] >= 0 && modsdx2[j] - modsdx2[i] > 2)
	    {
		dviol = distmat[modsdx2[i]][modsdx2[j]][CBATOM][CBATOM] - DCUTOFF;

		if (dviol < 0.0F)
		    pc = 0.5F + 0.5F * (1.0F - expf(-SQR(dviol)));
		else
		    pc = 0.5F * expf(-SQR(dviol));

		x[n] = pc;
		y[n] = contmap[i][j];
		n++;
	    }

    for (i=0; i<n; i++)
    {
	xav += x[i];
	yav += y[i];
    }

    xav /= n;
    yav /= n;

    for (i=0; i<n; i++)
    {
	sxy += (x[i] - xav) * (y[i] - yav);
	sxx += SQR(x[i] - xav);
	syy += SQR(y[i] - yav);
    }

    correl = sxy / sqrtf(sxx*syy);

    printf("correl = %f\n", correl);

    freevec(x);
    freevec(y);

    return correl / sqrtf((1.0F-SQR(correl)) / (n - 2.0F));
}

/* Compute energy sums for final model */
void
sc_final(float *contactsc, float *perco_a, float *perco_b)
{
    int             i, j, n, nc, t, wt;
    short           resa, resb;
    float           e, pc, dviol;

    for (n = i = 0; i < seq1len; i++)
	if (modsdx[i] >= 0)
	    n++;

    *perco_a = 100.0 * n / seq1len;
    *perco_b = 100.0 * n / seq2len;

    *contactsc = ZERO;

    for (nc = i = 0; i < seq1len; i++)
	if (seq1[i] < GAP && modsdx[i] >= 0)
	    for (j = i + 1; j < seq1len; j++)
		if (seq1[j] < GAP && modsdx[j] >= 0 && contmap[modsdx[i]][modsdx[j]] > 0.0)
		{
		    resa = seq2[modsdx[i]];
		    resb = seq2[modsdx[j]];

		    t = modsdx[j] - modsdx[i];

		    if (t < 2)
			continue;

		    dviol = distmat[i][j][CBATOM][CBATOM] - DCUTOFF;

#if 1
		    if (dviol > 0.0F)
			*contactsc += contmap[modsdx[i]][modsdx[j]] * expf(-dviol);
		    else
			*contactsc += contmap[modsdx[i]][modsdx[j]];
#else
		    if (dviol > 0.0F)
			*contactsc -= contmap[modsdx[i]][modsdx[j]] * dviol * 0.005F;
		    else
			*contactsc += contmap[modsdx[i]][modsdx[j]];
#endif
		    nc++;
		}

    *contactsc /= seq2len;

    *contactsc = jonespointbiserial();
//    *contactsc = pointbiserial();
//    *contactsc = correlation();
}


/*
 * Trace back highest scoring path
 */
void
                trace(short *posa, short *posb, int mati, int matj,
		      int pati, int patj, int lasti, int lastj, short *n)
{
    int             pij = pat[pati][patj], i, j;

    for (i = lasti + 1; i < mati; i++)
    {
	*(++posa) = i;
	*(++posb) = 0;
	(*n)++;
    }
    for (j = lastj + 1; j < matj; j++)
    {
	*(++posa) = 0;
	*(++posb) = j;
	(*n)++;
    }
    *(++posa) = mati;
    *(++posb) = matj;
    (*n)++;

    if (!pij)
	return;

    if (pij == 1)
	trace(posa, posb, mati + 1, matj + 1, pati + 1, patj + 1,
	      mati, matj, n);
    if (pij < 1)
	trace(posa, posb, mati + 1, matj - pij, pati + 1, patj - pij,
	      mati, matj, n);
    if (pij > 1)
	trace(posa, posb, mati + pij, matj + 1, pati + pij, patj + 1,
	      mati, matj, n);
}

int
                seqscore(const char *seq1, const char *seq2, ALNS * aln, const int seq1len, const int seq2len)
{
    short          *posa, *posb;
    int             trace_back = (aln != NULL);
    int             now = 0, last = 1;
    int             pati, patj, mati, matj, i, j, k, l;
    int             toprows[MAXSEQLEN + 1], topcol, toprow;
    int             maxrows[MAXSEQLEN + 1], maxscore, maxflg = FALSE, maxcol,
	maxrow, diag, row, col, envclass, gap_open, gap_extend;
    int             mat[2][MAXSEQLEN + 1];

#if 0
    for (k=2,i=0; i<3; i++)
	for (j=0; j<3; j++)
	    ssmat[i][j] = SC_PARAM[k++];
#endif

    if (trace_back)
	for (i = 1; i <= seq1len; i++)
	    pat[i][seq2len] = 0;

    for (j = seq2len; j > 0; j--)
    {
	if (trace_back)
	    pat[seq1len][j] = 0;

	for (i = seq1len; i > 0; i--)
	{
	    /* Get matrix element */

	    mat[now][i] = pmat[j-1][i-1];

	    if (j != seq2len && i != seq1len)
	    {
		diag = mat[last][i + 1];

		maxrow = maxrows[i];
		toprow = toprows[i];

		envclass = tsstruc[i - 1];

		gap_open = GP_OPEN[envclass] * 0.4;
		gap_extend = GP_EXT[envclass] * 0.4;

		if (toprow)
		    row = maxrow - gap_open - (toprow - j) * gap_extend + gap_extend;
		else
		    row = maxrow - gap_open;
		if (topcol)
		    col = maxcol - gap_open - (topcol - i) * gap_extend + gap_extend;
		else
		    col = maxcol - gap_open;

		if (diag > col && diag > row)
		{
		    mat[now][i] += diag;
		    if (trace_back)
			pat[i][j] = 1;
		}
		else
		{
		    if (row > col)
		    {
			mat[now][i] += row;
			if (trace_back)
			    pat[i][j] = -(toprow - j);
		    }
		    else
		    {
			mat[now][i] += col;
			if (trace_back)
			    pat[i][j] = topcol - i;
		    }
		}

/*		printf("i=%d j=%d toprow=%d maxrow=%d mat=%d\n", i, j, toprow, maxrow, mat[now][i]); */

		if (diag > maxrows[i])
		{
		    maxrows[i] = diag;
		    toprows[i] = j + 1;
		}
		if (diag > maxcol)
		{
		    maxcol = diag;
		    topcol = i + 1;
		}

		if (i == 1 || j == 1)
		{
		    if (!maxflg || mat[now][i] > maxscore)
		    {
			maxflg = TRUE;
			maxscore = mat[now][i];
/*			printf("MAXSCORE = %d\n"); */
			if (trace_back)
			{
			    pati = i;
			    patj = j;
			    mati = matj = 1;
			    if (i == 1)
				matj = j;
			    if (j == 1)
				mati = i;
			}
		    }
		}
	    }
	    else if (j == seq2len)
	    {
		maxrows[i] = mat[now][i];
		toprows[i] = j;
	    }
	    else if (i == seq1len)
	    {
		maxcol = mat[now][i];
		topcol = i;
	    }
	}
	now = !now;
	last = !last;
    }

    if (!trace_back)
	return (maxscore/100);

    posa = aln->posn_a;
    posb = aln->posn_b;
    aln->length = 0;
    trace(posa, posb, mati, matj, pati, patj, 0, 0, &aln->length);
    posa += aln->length;
    posb += aln->length;
    if (*posa == seq1len)
    {
	for (i = *(posb) + 1; i <= seq2len; i++)
	{
	    *(++posa) = 0;
	    *(++posb) = i;
	    (aln->length)++;
	}
	if (aln->length > MAXSEQLEN)
	    fail("score : max. align length exceeded!");
	return (maxscore/100);
    }
    if (*posb == seq2len)
	for (i = *(posa) + 1; i <= seq1len; i++)
	{
	    *(++posb) = 0;
	    *(++posa) = i;
	    (aln->length)++;
	}
    if (aln->length > MAXSEQLEN)
	fail("score : max. align length exceeded!");

    return (maxscore/100);
}

/* Eigen decomposition code for symmetric matrix, taken/modified from the public
   domain Java Matrix library JAMA. */

#define hypot2(x, y) (sqrtf((x)*(x)+(y)*(y)))

/* Symmetric Householder reduction to tridiagonal form. */

void tred2(float **V, float *d, float *e, int NDIM)
{
    int i, j, k;

    for (j = 0; j < NDIM; j++)
	d[j] = V[NDIM-1][j];

    /* Householder reduction to tridiagonal form. */

    for (i = NDIM-1; i > 0; i--)
    {
	/* Scale to avoid under/overflow. */

	float scale = 0.0F;
	float h = 0.0F;

	for (k = 0; k < i; k++)
	    scale = scale + fabsf(d[k]);

	if (scale == 0.0F)
	{
	    e[i] = d[i-1];

	    for (j = 0; j < i; j++)
	    {
		d[j] = V[i-1][j];
		V[i][j] = 0.0F;
		V[j][i] = 0.0F;
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

	    float f = d[i-1];
	    float g = sqrtf(h);

	    if (f > 0)
		g = -g;

	    e[i] = scale * g;
	    h = h - f * g;
	    d[i-1] = f - g;

	    for (j = 0; j < i; j++)
		e[j] = 0.0F;

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

	    f = 0.0F;

	    for (j = 0; j < i; j++)
	    {
		e[j] /= h;
		f += e[j] * d[j];
	    }

	    float hh = f / (h + h);

	    for (j = 0; j < i; j++)
		e[j] -= hh * d[j];

	    for (j = 0; j < i; j++)
	    {
		f = d[j];
		g = e[j];

		for (k = j; k <= i-1; k++)
		    V[k][j] -= (f * e[k] + g * d[k]);

		d[j] = V[i-1][j];

		V[i][j] = 0.0F;
	    }
	}

	d[i] = h;
    }


    /* Accumulate transformations. */

    for (i = 0; i < NDIM-1; i++)
    {
	V[NDIM-1][i] = V[i][i];

	V[i][i] = 1.0F;

	float h = d[i+1];

	if (h != 0.0F)
	{
	    for (k = 0; k <= i; k++)
		d[k] = V[k][i+1] / h;

	    for (j = 0; j <= i; j++)
	    {
		float g = 0.0F;

		for (k = 0; k <= i; k++)
		    g += V[k][i+1] * V[k][j];

		for (k = 0; k <= i; k++)
		    V[k][j] -= g * d[k];
	    }
	}

	for (k = 0; k <= i; k++)
	    V[k][i+1] = 0.0F;
    }

    for (j = 0; j < NDIM; j++)
    {
	d[j] = V[NDIM-1][j];
	V[NDIM-1][j] = 0.0F;
    }

    V[NDIM-1][NDIM-1] = 1.0F;

    e[0] = 0.0F;
}


/* Symmetric tridiagonal QL algorithm. */

void tql2(float **V, float *d, float *e, int NDIM)
{
    int i, j, k, l, m;
    float f = 0.0F;
    float tst1 = 0.0F;
    float eps = FLT_EPSILON;

    for (i = 1; i < NDIM; i++)
	e[i-1] = e[i];

    e[NDIM-1] = 0.0F;

    for (l = 0; l < NDIM; l++)
    {
	/* Find small subdiagonal element */

	tst1 = MAX(tst1,fabsf(d[l]) + fabsf(e[l]));

	int m = l;

	while (m < NDIM)
	{
	    if (fabsf(e[m]) <= eps*tst1)
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

		float g = d[l];
		float p = (d[l+1] - g) / (2.0F * e[l]);
		float r = hypot2(p,1.0F);

		if (p < 0)
		    r = -r;

		d[l] = e[l] / (p + r);
		d[l+1] = e[l] * (p + r);

		float dl1 = d[l+1];
		float h = g - d[l];

		for (i = l+2; i < NDIM; i++)
		    d[i] -= h;

		f = f + h;

		/* Implicit QL transformation. */

		p = d[m];

		float c = 1.0F;
		float c2 = c;
		float c3 = c;
		float el1 = e[l+1];
		float s = 0.0F;
		float s2 = 0.0F;

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

	    } while (fabsf(e[l]) > eps*tst1);
	}

	d[l] = d[l] + f;
	e[l] = 0.0F;
    }

    /* Sort eigenvalues and corresponding vectors. */

    for (i = 0; i < NDIM-1; i++)
    {
	float p = d[i];

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

/* Return eigenvectors (in V) and eigenvalues (in D) for symmetric matrix A */
void eigen_decomposition(float **A, float *d, float **V, int NDIM)
{
    int i, j;
    float e[NDIM];

    for (i = 0; i < NDIM; i++)
	for (j = 0; j < NDIM; j++)
	    V[i][j] = A[i][j];

    tred2(V, d, e, NDIM);
    tql2(V, d, e, NDIM);
}

/* Align sequence to structural template */
void
                seqfit(float *contactsc, float *perco_a, float *perco_b, ALNS * finaltplt)
{
    int             i, j, k, l, t, tt, neig;
    float score, bestscore, bestsignsc, bestsign;
    float sign[MAXSEQLEN], sum;
    ALNS            tplt;

    /* Build structure gap penalty array */
    for (i = 0; i < nsst; i++)
	if (!i)
	    for (j = 0; j < sstlist[0].start; j++)
		maxgplen[j] = sstlist[0].start - j;
	else
	    for (j = sstlist[i - 1].start + sstlist[i - 1].length; j < sstlist[i].start; j++)
		maxgplen[j] = sstlist[i].start - j;
    for (j = sstlist[nsst - 1].start + sstlist[nsst - 1].length; j < seq1len; j++)
	maxgplen[j] = BIG;

    bestscore = -VBIG;

    neig = MIN(MIN(NEIGEN, seq1len), seq2len);

    for (t=0; t < neig; t++)
	sign[t] = 0;

    for (tt=1; tt <= neig; tt++)
    for (t=0; t < tt; t++)
    {
	bestsignsc = -VBIG;

	for (sign[t] = -1.0F; sign[t] <= 1.0F; sign[t] += 2.0F)
	{
	    for (j = 0; j < seq2len; j++)
		for (i = 0; i < seq1len; i++)
		{
		    for (sum=k=0; k < neig; k++)
			sum += eigenvec1[i][k] * eigenvec2[j][k] * eigenval1[k] * eigenval2[k] * sign[k];
//		printf("%d %d %f\n", i, j, sum);
		    pmat[j][i] = sum * SCALE;
		}

	    if (ssstflg)
		for (j = 0; j < seq2len; j++)
		    if (ssstrel[j] >= MINSSTREL)
			for (i = 0; i < seq1len; i++)
			    pmat[j][i] += ssmat[ssstruc[j-1]][tsstruc[i-1]];

	    seqscore(seq1, seq2, &tplt, seq1len, seq2len);

//	printf("%s score = %d\n", brkid, seqscore(seq1, seq2, &tplt, seq1len, seq2len));

	    for (k = 0; k < seq1len; k++)
		modsdx[k] = -1;
	    for (k = 0; k < seq2len; k++)
		modsdx2[k] = -1;
	    for (k = 1; k <= tplt.length; k++)
		if (tplt.posn_a[k] > 0 && tplt.posn_b[k] > 0)
		{
		    modsdx[tplt.posn_a[k] - 1] = tplt.posn_b[k] - 1;
		    modsdx2[tplt.posn_b[k] - 1] = tplt.posn_a[k] - 1;
		}

	    for (score = i = 0; i < seq2len; i++)
		for (j = i+1; j < seq2len; j++)
		    if (modsdx2[i] >= 0 && modsdx2[j] >= 0)
			if (distmat[modsdx2[i]][modsdx2[j]][CBATOM][CBATOM] < DCUTOFF)
			    score += contmap[i][j];

	    if (score > bestsignsc)
	    {
		bestsignsc = score;
		bestsign = sign[t];
	    }

	    if (score > bestscore)
	    {
		/* printf("%d %f contact score = %f\n", t, sign[t], score);*/
		bestscore = score;
		memcpy(finaltplt, &tplt, sizeof(ALNS));
	    }
	}
	sign[t] = bestsign;
    }


    for (k = 0; k < seq1len; k++)
	modsdx[k] = -1;
    for (k = 0; k < seq2len; k++)
	modsdx2[k] = -1;
    for (k = 1; k <= finaltplt->length; k++)
	if (finaltplt->posn_a[k] > 0 && finaltplt->posn_b[k] > 0)
	{
	    modsdx[finaltplt->posn_a[k] - 1] = finaltplt->posn_b[k] - 1;
	    modsdx2[finaltplt->posn_b[k] - 1] = finaltplt->posn_a[k] - 1;
	}

    if (prtalnflg)
    {
	if (htmlflg)
	    htmlalign(finaltplt);
	prtalign(finaltplt);
    }

    sc_final(contactsc, perco_a, perco_b);

#if 0
    putchar('\n');
    for (i = 0; i < seq1len; i++)
    {
	putchar(rescodes[seq1[i]]);
	if (modsdx[i] >= 0)
	    printf("%c %6.3f\n", rescodes[seq2[modsdx[i]]], e_contrib[i]);
	else
	    printf("-\n");
    }
#endif

    if (verbflg)
	printf("Score_contact = %f\n", *contactsc);
}

/* Change some DSSP Gs to Hs */
void
                helix310(char *struc, int len)
{
    int             i, j, run;

    for (i = 0; i < len; i++)
	if (struc[i] == 'G')
	{
	    for (run = 1; run <= len - i && struc[i + run] == 'G'; run++);
	    if (run >= 4 || (run == 3 && ((i && struc[i - 1] == 'H') || (i < len-run-3 && struc[i + run] == 'H'))))
	    {
		for (j = i; j < i + run; j++)
		    struc[j] = 'H';
		i += run - 1;
	    }
	}
}

/* Smooth DSSP secondary structure definitions */
void
                smooth_sst(void)
{
    int             i;

    for (i = 1; i < seq1len - 1; i++)
	if (tsstruc[i - 1] == tsstruc[i + 1])
	    tsstruc[i] = tsstruc[i - 1];
}

/* Extend DSSP secondary structure definitions */
void
                ext_sst(void)
{
    int             i;

    for (i = 1; i < seq1len - 1; i++)
	if (!tsstruc[i] && !tsstruc[i + 1] && tsstruc[i - 1])
	{
	    tsstruc[i] = tsstruc[i - 1];
	    i++;
	}
	else if (!tsstruc[i - 1] && !tsstruc[i] && tsstruc[i + 1])
	    tsstruc[i] = tsstruc[i + 1];
}

/* Locate secondary structures */
void
                loc_sst(void)
{
    int             i, istart, l;

    nsst = 0;
    for (i = 0; i < seq1len; i++)
	if (tsstruc[i])
	{
	    istart = i;
	    while (i < seq1len && tsstruc[i] == tsstruc[istart])
		i++;
	    l = i - istart;
	    if (l < 2)
		continue;
	    sstlist[nsst].start = istart;
	    sstlist[nsst].length = l;
	    switch (tsstruc[istart])
	    {
	    case HELIX:
		sstlist[nsst].type = HELIX;
		if (verbflg)
		    printf("HELIX  at %3d, length %2d.\n", istart, l);
		break;
	    case STRAND:
		sstlist[nsst].type = STRAND;
		if (verbflg)
		    printf("STRAND at %3d, length %2d.\n", istart, l);
		break;
	    default:
		fail("Unknown secondary structure code!");
	    }
	    sstlen += l;
	    nsst++;
	}
    if (verbflg)
	putchar('\n');

    for (i = 0; i < seq1len; i++)
	sstidx[i] = -1;
    for (i = 0; i < nsst; i++)
	for (l = 0; l < sstlist[i].length; l++)
	    sstidx[sstlist[i].start + l] = i;
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


/* Build structural template */
void
maketemplate(FILE *dfp, int distflg)
{
    char            buf[512], whichatm[MAXSEQLEN], *cp, ri[10];
    int             i, j, k, nres, natoms, n = 0,
                    consv;
    float           dv, cca[3], nca[3], xx[3], yy[3], sx, sy, rmean = ZERO, sum, targf[20];
    short           at1, at2;

    seq1len = i = natoms = 0;

    if (!fgets(buf, 512, dfp))
	return;

    while (!feof(dfp))
    {
	if (!fgets(buf, 512, dfp))
	    break;

	if (distflg)
	    if (sscanf(buf + 12, "%f%f%f", phi + i, psi + i, omega + i) != 3 ||
		sscanf(buf + 39, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
		       &x[i][NATOM], &y[i][NATOM], &z[i][NATOM],
		       &x[i][CAATOM], &y[i][CAATOM], &z[i][CAATOM],
		       &x[i][CATOM], &y[i][CATOM], &z[i][CATOM],
		       &x[i][OATOM], &y[i][OATOM], &z[i][OATOM],
		       &x[i][CBATOM], &y[i][CBATOM], &z[i][CBATOM]) != 15)
		fail("Bad tdb file!");

	sscanf(buf+174, "%s", ri);

	if (resid[i] != NULL)
	    free(resid[i]);

	resid[i] = strdup(ri);

	cp = buf + 181;
	for (j = 0; j < 20; j++,cp+=6)
	{
	  consv = atoi(cp);
	  tpltsc1[i][j] = consv;
	}

	if (buf[5] == '-')
	{
	    seq1[i] = GAP;
	    whichatm[i] = 0;
	    tsstruc[i] = COIL;
	}
	else
	{
	    seq1[i] = aanum(buf[5]);
	    tsstruc[i] = buf[7];
	    sscanf(buf + 9, "%d", &j);
	    trelacc[i] = j;
	    whichatm[i] = 31;
	}
	++i;
    }

    seq1len = nres = i;

#if 0
    for (i = 0; i < nres; i++)
    {
	x[i][CBATOM] = 1.56 * (x[i][CBATOM] - x[i][CAATOM]) + x[i][CAATOM];
	y[i][CBATOM] = 1.56 * (y[i][CBATOM] - y[i][CAATOM]) + y[i][CAATOM];
	z[i][CBATOM] = 1.56 * (z[i][CBATOM] - z[i][CAATOM]) + z[i][CAATOM];
    }
#endif

    helix310(tsstruc, nres);
    for (i = 0; i < nres; i++)
	switch (tsstruc[i])
	{
	case 'H':
	    tsstruc[i] = HELIX;
	    break;
	case 'E':
	case 'A':
	case 'P':
	    tsstruc[i] = STRAND;
	    break;
	default:
	    tsstruc[i] = COIL;
	    break;
	}

    for (i = 0; i < nres; i++)
    {
	for (k=0; k<20; k++)
	    targf[k] = dbaaf[k] * exp(0.00318 * tpltsc1[i][k]);

	for (sum=k=0; k<20; k++)
	    sum += targf[k];

	for (k=0; k<20; k++)
	    targf1[i][k] = targf[k] / sum;
    }

    if (!distflg)
	return;

    /* Check atoms */
    for (i = 0; i < nres; i++)
    {
	if (seq1[i] != GAP && !(whichatm[i] & (1 << NATOM)))
	{
	    printf("FATAL: Missing N atom in %d!\n", i + 1);
	    exit(1);
	}
	if (!(whichatm[i] & (1 << CAATOM)))
	{
	    if (verbflg)
		printf("WARNING: Missing CA atom in %d!\n", i + 1);
	    seq1[i] = GAP;
	}
	if (!(whichatm[i] & (1 << CBATOM)))
	{
	    if (!(whichatm[i] & (1 << CAATOM)) || !(whichatm[i] & (1 << CATOM)) || !(whichatm[i] & (1 << NATOM)))
	    {
		/* Not much left of this residue! */
		if (verbflg)
		    printf("WARNING: Missing main chain atom in %d!\n", i + 1);
		seq1[i] = GAP;
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
	    if (seq1[i] != GLY && verbflg)
		fprintf(stderr, "WARNING: dummy CB atom constructed for %d (%s)!\n", i + 1, rnames[seq1[i]]);
	}
    }

    if (!nres)
	return;

    distmat = allocmat(nres, nres, sizeof(**distmat), TRUE);

    /* Calculate interatomic distance template */
    for (i = 0; i < nres; i++)
	tooi[i] = 0;

    rmax = ZERO;

    for (i = 0; i < nres; i++)
	for (j = i; j < nres; j++)
	    for (at1 = CAATOM; at1 <= NATOM; at1++)
		for (at2 = CAATOM; at2 <= NATOM; at2++)
		    if (i != j || at1 != at2)
			if ((whichatm[i] & (1 << at1)) && (whichatm[j] & (1 << at2)))
			{
			    dv = sqrtf(SQR(x[i][at1] - x[j][at2]) + SQR(y[i][at1] - y[j][at2]) + SQR(z[i][at1] - z[j][at2]));
			    distmat[i][j][at1][at2] = distmat[j][i][at2][at1] = dv;
			}
			else
			    distmat[i][j][at1][at2] = distmat[j][i][at2][at1] = 999.0;
		    else
			distmat[i][j][at1][at2] = ZERO;

    rmean /= n;

    for (i = 0; i < seq1len; i++)
	coreflg[i] = tsstruc[i] != COIL;

    for (i = 0; i < seq1len; i++)
	if (coreflg[i] && (!i || !coreflg[i - 1]) && (i == seq1len - 1 || !coreflg[i + 1]))
	    coreflg[i] = 0;

#if 0
    for (i = 0; i < seq1len; i++)
	printf("%c %c %d\n", rescodes[seq1[i]], sscodes[tsstruc[i]], tooi[i]);
#endif
}

#if 1
/* Write PDB file */
void            writepdb(char *fname, char *brkid)
{
    FILE           *ofp;
    int             atnum, i;
    float coorderr;

    ofp = fopen(fname, "w");
    if (ofp != NULL)
    {
	fprintf(ofp, "HEADER  %s BASED EiGenTHREADER MODEL\n", brkid);
	for (atnum = i = 0; i < seq2len; i++)
	    if (modsdx2[i] < 0 || seq1[modsdx2[i]] == GAP)
	    {
#if 1
		fprintf(ofp, "REMARK  UNALIGNED ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " N  ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		fprintf(ofp, "REMARK  UNALIGNED ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " CA ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		fprintf(ofp, "REMARK  UNALIGNED ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " C  ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		fprintf(ofp, "REMARK  UNALIGNED ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " O  ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		if (seq2[i] != GLY)
		    fprintf(ofp, "REMARK  UNALIGNED ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
			    ++atnum, " CB ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
#endif
#if 0
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " N  ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " CA ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " C  ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " O  ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		if (seq2[i] != GLY)
		    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
			    ++atnum, " CB ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
#endif
	    }
	    else
	    {
		coorderr = 9.9;
		if (coorderr < 1.0)
		    coorderr = 1.0;
		if (coorderr > 10.0)
		    coorderr = 10.0;
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00 %5.2f\n",
			++atnum, " N  ", rnames[seq2[i]], i + 1, x[modsdx2[i]][NATOM], y[modsdx2[i]][NATOM], z[modsdx2[i]][NATOM], coorderr);
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00 %5.2f\n",
			++atnum, " CA ", rnames[seq2[i]], i + 1, x[modsdx2[i]][CAATOM], y[modsdx2[i]][CAATOM], z[modsdx2[i]][CAATOM], coorderr);
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00 %5.2f\n",
			++atnum, " C  ", rnames[seq2[i]], i + 1, x[modsdx2[i]][CATOM], y[modsdx2[i]][CATOM], z[modsdx2[i]][CATOM], coorderr);
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00 %5.2f\n",
			++atnum, " O  ", rnames[seq2[i]], i + 1, x[modsdx2[i]][OATOM], y[modsdx2[i]][OATOM], z[modsdx2[i]][OATOM], coorderr);
		if (seq2[i] != GLY)
		    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00 %5.2f\n",
			++atnum, " CB ", rnames[seq2[i]], i + 1, x[modsdx2[i]][CBATOM], y[modsdx2[i]][CBATOM], z[modsdx2[i]][CBATOM], coorderr);
	    }
	fprintf(ofp, "TER\nEND\n");
	fclose(ofp);
    }
}
#else
/* Write torsions file */
void            writepdb(char *fname)
{
    FILE           *ofp;
    int             i;

    ofp = fopen(fname, "w");
    if (ofp != NULL)
	for (i = 0; i < seq2len; i++)
	    if (modsdx2[i] < 0 || seq1[modsdx2[i]] == GAP)
		fprintf(ofp, "%7.3f %7.3f %7.3f\n", 999.9, 999.9, 999.9);
	    else
		fprintf(ofp, "%7.3f %7.3f %7.3f\n", phi[modsdx2[i]], psi[modsdx2[i]], omega[modsdx2[i]]);
    fclose(ofp);
}
#endif

/* Write profile comparison file */
void            writeprofile(char *fname)
{
    FILE           *ofp;
    int             i, k;

    ofp = fopen(fname, "w");
    if (ofp != NULL)
    {
	for (i = 0; i < seq2len; i++)
	    if (modsdx2[i] >= 0 && seq1[modsdx2[i]] != GAP)
	    {
		fprintf(ofp, "%2d", seq1[modsdx2[i]]);
		for (k=0; k<20; k++)
		    fprintf(ofp, " %d", tpltsc1[modsdx2[i]][k]);
		fprintf(ofp, "\n%2d", seq2[i]);
		for (k=0; k<20; k++)
		    fprintf(ofp, " %d", tpltsc2[i][k]);
		fprintf(ofp, "\n");
	    }

	fclose(ofp);
    }
}

/* Read PSI AA frequency data */
int             getpsi(FILE * lfil)
{
    int             i, j, k, naa, transtab[20];
    float           sum, targf[20];
    char            buf[256];

    if (!fgets(buf, 256, lfil))
	fail("Bad PSI format!");

    for (j = 0, i = 0; i < 80 && buf[i] && j < 20; i++)
	if (isalpha(buf[i]))
	    transtab[j++] = aanum(buf[i]);

    if (j != 20)
	fail("Bad PSI format!");

    naa = 0;
    while (!feof(lfil))
    {
	if (!fgets(buf, 256, lfil))
	    break;
	if (sscanf(buf, "%d", &i) != 1)
	    break;
	if (i != naa + 1)
	    fail("Bad PSI format!");
	seq2[naa] = buf[6];
	if (sscanf(buf + 7, "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d", &tpltsc2[naa][transtab[0]], &tpltsc2[naa][transtab[1]], &tpltsc2[naa][transtab[2]], &tpltsc2[naa][transtab[3]], &tpltsc2[naa][transtab[4]], &tpltsc2[naa][transtab[5]], &tpltsc2[naa][transtab[6]], &tpltsc2[naa][transtab[7]], &tpltsc2[naa][transtab[8]], &tpltsc2[naa][transtab[9]], &tpltsc2[naa][transtab[10]], &tpltsc2[naa][transtab[11]], &tpltsc2[naa][transtab[12]], &tpltsc2[naa][transtab[13]], &tpltsc2[naa][transtab[14]], &tpltsc2[naa][transtab[15]], &tpltsc2[naa][transtab[16]], &tpltsc2[naa][transtab[17]], &tpltsc2[naa][transtab[18]], &tpltsc2[naa][transtab[19]]) != 20)
	    fail("Bad PSI format!");
#if 0
	if (sscanf(buf + 70, "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d", &targf2[naa][transtab[0]], &targf2[naa][transtab[1]], &targf2[naa][transtab[2]], &targf2[naa][transtab[3]], &targf2[naa][transtab[4]], &targf2[naa][transtab[5]], &targf2[naa][transtab[6]], &targf2[naa][transtab[7]], &targf2[naa][transtab[8]], &targf2[naa][transtab[9]], &targf2[naa][transtab[10]], &targf2[naa][transtab[11]], &targf2[naa][transtab[12]], &targf2[naa][transtab[13]], &targf2[naa][transtab[14]], &targf2[naa][transtab[15]], &targf2[naa][transtab[16]], &targf2[naa][transtab[17]], &targf2[naa][transtab[18]], &targf2[naa][transtab[19]]) != 20)
	    fail("Bad PSI format!");
#endif
	naa++;
    }

    for (i=0; i<naa; i++)
    {
	for (k=0; k<20; k++)
	    targf[k] = dbaaf[k] * exp(0.318 * tpltsc2[i][k]);

	for (sum=k=0; k<20; k++)
	    sum += targf[k];

	for (k=0; k<20; k++)
	    targf2[i][k] = targf[k] / sum;

	for (k=0; k<20; k++)
	    tpltsc2[i][k] *= 100;
    }

    psidata = TRUE;

    return naa;
}

/* Read PSI mtx PSSM data */
int             getmtx(FILE *lfil)
{
    int             i, j, k, naa;
    float           sum, targf[20];
    char            buf[256];

    rewind(lfil);

    if (fscanf(lfil, "%d", &naa) != 1)
	fail("Bad mtx file - sequence length not found!");

    if (naa > MAXSEQLEN)
	fail("Input sequence too long!");

    if (fscanf(lfil, "%s", seq2) != 1)
      fail("Bad mtx file - no sequence!");

    while (!feof(lfil))
    {
	if (!fgets(buf, 65536, lfil))
	    fail("Bad mtx file!");
	if (!strncmp(buf, "-32768 ", 7))
	{
	    for (j=0; j<naa; j++)
	    {
		if (sscanf(buf, "%*d%d%*d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%*d%d", &tpltsc2[j][ALA],  &tpltsc2[j][CYS], &tpltsc2[j][ASP],  &tpltsc2[j][GLU],  &tpltsc2[j][PHE],  &tpltsc2[j][GLY],  &tpltsc2[j][HIS],  &tpltsc2[j][ILE],  &tpltsc2[j][LYS],  &tpltsc2[j][LEU],  &tpltsc2[j][MET],  &tpltsc2[j][ASN],  &tpltsc2[j][PRO],  &tpltsc2[j][GLN],  &tpltsc2[j][ARG],  &tpltsc2[j][SER],  &tpltsc2[j][THR],  &tpltsc2[j][VAL],  &tpltsc2[j][TRP],  &tpltsc2[j][TYR]) != 20)
		    fail("Bad mtx format!");
		if (!fgets(buf, 65536, lfil))
		    break;
	    }
	}
    }

    for (i=0; i<naa; i++)
    {
	for (k=0; k<20; k++)
	    targf[k] = dbaaf[k] * exp(0.00318 * tpltsc2[i][k]);

	for (sum=k=0; k<20; k++)
	    sum += targf[k];

	for (k=0; k<20; k++)
	    targf2[i][k] = targf[k] / sum;
    }

    psidata = TRUE;

    return naa;
}


/* Read PSIPRED VFORMAT prediction data */
int             getpsipredv(FILE * lfil)
{
    int             naa;
    float confc, confh, confe;
    char            buf[256];

    if (!fgets(buf, 256, lfil))
	fail("Bad PSIPRED VFORMAT file!");

    naa = 0;
    while (!feof(lfil))
    {
	if (!fgets(buf, 256, lfil))
	    break;
	if (sscanf(buf+10, "%f%f%f", &confc, &confh, &confe) != 3)
	    break;
	seq2[naa] = buf[5];
	switch (buf[7])
	{
	case 'H':
	    ssstruc[naa] = HELIX;
	    break;
	case 'E':
	    ssstruc[naa] = STRAND;
	    break;
	default:
	    ssstruc[naa] = COIL;
	    break;
	}
	ssstrel[naa++] = 100*(2*MAX(MAX(confc, confh),confe)-(confc+confh+confe)-MIN(MIN(confc, confh),confe));
    }

    ssstflg = TRUE;

    if (verbflg)
	puts("Parsed PSIPRED output.");

    return naa;
}

/* Read PSIPRED HFORMAT prediction data */
int             getpsipredh(FILE * lfil)
{
    int             naa, nsst, nrel;
    char            buf[256], *p;

    rewind(lfil);

    while (!feof(lfil))
    {
	if (!fgets(buf, 256, lfil))
	    fail("Bad PSIPRED format!");
	if (!strncmp(buf, "Conf:", 5))
	    break;
    }

    nsst = nrel = naa = 0;
    while (!feof(lfil))
    {
	p = buf+5;
	while (*++p)
	{
	    if (isdigit(*p))
		ssstrel[nrel++] = (*p - '0')*10;
	}
	if (!fgets(buf, 256, lfil))
	    fail("Bad PSIPRED format!");
	p = buf+5;
	while (*++p)
	    if (isalpha(*p))
		switch (*p)
		{
		case 'H':
		    ssstruc[nsst++] = HELIX;
		    break;
		case 'E':
		    ssstruc[nsst++] = STRAND;
		    break;
		default:
		    ssstruc[nsst++] = COIL;
		    break;
		}

	if (!fgets(buf, 256, lfil))
	    fail("Bad PSIPRED format!");
	p = buf+5;
	while (*++p)
	{
	    if (isalpha(*p))
		seq2[naa++] = *p;
	    else if (isdigit(*p))
		seq2[naa++] = -(*p - '0' + 1);
	}
	if (!fgets(buf, 256, lfil))
	    fail("Bad PSIPRED format!");
	while (!feof(lfil))
	{
	    if (!fgets(buf, 256, lfil))
		break;
	    if (!strncmp(buf, "Conf:", 5))
		break;
	}
    }

    if (nsst != nrel || naa != nsst)
    {
	fprintf(stderr, "Incorrect PSIPRED prediction file format!\n");
	return -1;
    }

    ssstflg = TRUE;

    if (verbflg)
	puts("Parsed PSIPRED output.");

    return naa;
}

/* Read contact prediction */
int readcontacts(FILE *ifp)
{
    int i, j, n=0;
    float dcut, r;
    char buf[1024];

    contmap = allocmat(seq2len, seq2len, sizeof(float), TRUE);

    for (;;)
    {
	if (!fgets(buf, 256, ifp))
	    break;

	if (isalpha(buf[0]))
	    continue;

	if (sscanf(buf, "%d%d%*s%f%f", &i, &j, &dcut, &r) != 4)
	    break;

	contmap[i-1][j-1] = contmap[j-1][i-1] = r;
	n++;
    }

    return n;
}


/*
 * This routine will read in one sequence from a database file. The sequence
 * can be in any of the supported formats. Returns length of sequence.
 */
int
                getseq(char *dbname, char *dseq, FILE * lfil)
{
    int             i, j, len;
    short           badln, fformat = -1;
    enum
    {
	unknown, embl, genbank, staden, fastp, codata, owl, intelgen, gcg
    };
    char            temp[8192], split;
    int             offset;

    offset = j = 0;

    if (fgets(temp, 8192, lfil) && !strncmp(temp, "# PSIPRED H", 11))
	return getpsipredh(lfil);

    if (!strncmp(temp, "# PSIPRED V", 11))
	return getpsipredv(lfil);

    rewind(lfil);

    if (fgets(temp, 8192, lfil) && !strncmp(temp, "Conf:", 5))
	return getpsipredh(lfil);

    if (fgets(temp, 8192, lfil) && !strncmp(temp, "Conf:", 5))
	return getpsipredh(lfil);

    rewind(lfil);

    if (fgets(temp, 8192, lfil) && isdigit(temp[0]))
	return getmtx(lfil);

    if (fgets(temp, 8192, lfil) && !strncmp(temp, "Last position-specific", 22))
	return getpsi(lfil);

    rewind(lfil);

    if (fgets(temp, 8192, lfil) && fgets(temp, 8192, lfil))
    {
	if (!strncmp(temp + 14, "999.999", 7))
	{
	    ssstflg = TRUE;
	    do
	    {
		dseq[j] = temp[5];
		switch (temp[7])
		{
		case 'H':
		    ssstruc[j++] = HELIX;
		    break;
		case 'E':
		case 'A':
		case 'P':
		    ssstruc[j++] = STRAND;
		    break;
		default:
		    ssstruc[j++] = COIL;
		    break;
		}
	    }
	    while (fgets(temp, 8192, lfil));
	    dseq[j] = '\0';
	    return j;
	}
    }

    rewind(lfil);

    if (!fgets(temp, 8192, lfil))
	return (-1);

    /* Look for old-style PSI file */
    if (!strncmp(temp, "# PSI", 5))
	return getpsi(lfil);

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
	sscanf(strstr(temp, "of:") + 3, "%s", dbname);
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
	fprintf(stderr, "\nWARNING: sequence %s over %d long; truncated!\n",
	       dbname, MAXSEQLEN);

    dseq[j] = '\0';
    return (j);
}

void
                usage(char *cmdname)
{
    printf("usage: %s {options} sequence-file contacts-file output-file [list-file]\n", cmdname);
    puts("\nwhere options are as follows:");
    puts("-cnnn = set pair potential cutoff distance (default: 10.0)");
    puts("-rnnn = set number of sequence shuffles (default: 0)");
    puts("-ynnn = set secondary structure score weighting (default: 100)");
    puts("-hnnn = set minimum reliability for SST filter (default: 0.5)");
    puts("-Xa,b,c = set CHE gap opening penalties (default: 11)");
    puts("-Od,e,f = set CHE gap extension penalties (default: 1)");
    puts("-Ffnm = read secondary structure data from file fnm");
    puts("-Hfnm = output HTML alignment to file fnm");
    puts("-S = include predicted secondary structure in alignment output");
    puts("-u = use contact number solvation scoring");
    puts("-m{xxx} = generate PDB files (xxx optional prefix)");
    puts("-p{q} = print alignments (q selects quick align mode)");
    puts("-v = select Verbose Mode");
    puts("\nExample: threader -d100 -c12 -v -p globin.seq output.fil");
    exit(1);
}

float           cputime()
{
    float           tot;
    static float    ftime = -1.0, tickspersec = -1.0;
    struct tms      t;

    if (tickspersec < 0.0F)
	tickspersec = sysconf(_SC_CLK_TCK);

    (void) times(&t);

    tot = (float) (t.tms_utime + t.tms_stime) / tickspersec;

    if (ftime < 0.0F)
	ftime = tot;

    return tot - ftime;
}

int             main(int argc, char **argv)
{
    int             aa, i, j, k, l, nn = 0, filtered, start1, start2, end1, end2, alnmode1, alnmode2, nid;
    char           *cmdstr, *cp;
    char            desc[512], wpname[80];
    float           tpair, tsolv, psum, psumsq, ssum, ssumsq, smax, pmax,
	nwsc, nwsc1, nwsc2, nwsc3, nwsc4, dum, perid;
    float contact_av, contact_sd, **inputmat;
    FILE           *ifp, *ofp, *tfp, *efp, *ssfp;
    char            tdbname[512], eigname[512], buf[MAXSEQLEN], templname[512];
    ALNS            tplt;
    struct rlimit limits;

#ifndef DTJCOPY
    if (verbflg)
    {
	printf("EiGenTHREADER - Contact-driven Protein Sequence Threading Program\n");
	printf("Build date : %s\n",__DATE__);
	printf("Program written by David T. Jones\n");
	printf("Copyright (C) 2015 University College London\n");
	printf("Portions Copyright (C) 1997 University of Warwick\n");
	printf("Portions Copyright (C) 1989,1991 David T. Jones\n\n");
    }
#endif

    getrlimit(RLIMIT_STACK, &limits);

    if (limits.rlim_cur < 7340032)
        fail("Stack too small! (type 'limit stacksize 7m' at command prompt)");

    for (cmdstr = *argv++, argc--; argc && **argv == '-'; argv++, argc--)
	switch (*(*argv + 1))
	{
	case 'r':
	    SHUFCOUNT = atoi(*argv + 2);
	    break;
	case 'h':
	    MINSSTREL = 100.0*atof(*argv + 2);
	    if (MINSSTREL > 100)
	      fail("Max. secondary structure reliability should be <= 1.0!");
	    break;
	case 'B':
	    cp = *argv + 2;
	    for (i=0; i<20 && cp; i++)
		SC_PARAM[i] = atof(strsep(&cp, "/:;,"));
	    break;
	case 'O':
	    cp = *argv + 2;
	    for (i=0; i<3 && cp; i++)
		GP_OPEN[i] = 0.5 + 100.0 * atof(strsep(&cp, "/:;,"));
	    break;
	case 'X':
	    cp = *argv + 2;
	    for (i=0; i<3 && cp; i++)
		GP_EXT[i] = 0.5 + 100.0 * atof(strsep(&cp, "/:;,"));
	    break;
	case 'C':
	    SCUTOFF = atof(*argv + 2);
	    break;
	case 'S':
	    seqssflg = TRUE;
	    break;
	case 'g':
	    filtgapflg = TRUE;
	    break;
	case 'y':
	    ssscore = atoi(*argv + 2);
	    break;
	case 'c':
	    DCUTOFF = atof(*argv + 2);
	    break;
	case 'z':
	    SCALE = atof(*argv + 2);
	    break;
	case 't':
	    NEIGEN = atoi(*argv + 2);
	    break;
	case 'm':
	    motifflg = TRUE;
	    strcpy(motpref, *argv + 2);
	    break;
	case 'p':
	    prtalnflg = TRUE;
	    if (argv[0][2] == 'q')
		quickalnflg = TRUE;
	    if (argv[0][2] == 'm')
		mod3flg = TRUE;
	    break;
	case 'H':
	    htmname = strdup(*argv + 2);
	    htmlflg = TRUE;
	    break;
	case 'F':
	    ssfname = strdup(*argv + 2);
	    ssstflg = TRUE;
	    break;
	case 'u':
	    ooisolv = TRUE;
	    break;
	case 'v':
	    verbflg = TRUE;
	    break;
	default:
	    usage(cmdstr);
	}

    if (argc < 2)
	usage(cmdstr);

    randomise();
    printf("Analysing : %s\n", argv[0]);
    ifp = fopen(argv[0], "r");
    if (!ifp)
	fail("Unable to open seq file!");
    seq2len = getseq(desc, seq2, ifp);

    if (ssstflg)
    {
	ssfp = fopen(ssfname, "r");
	if (!ssfp)
	    fail("Cannot open PSIPRED file!");
	ssstflg = FALSE;
	if (fgets(buf, 160, ssfp) && !strncmp(buf, "# PSIPRED H", 11))
	    if (getpsipredh(ssfp) != seq2len)
		fail("Length mismatch between sequence and PSIPRED file!");
	if (!strncmp(buf, "# PSIPRED V", 11))
	    if (getpsipredv(ssfp) != seq2len)
		fail("Length mismatch between sequence and PSIPRED file!");
	if (!strncmp(buf, "Conf:", 5))
	    if (getpsipredh(ssfp) != seq2len)
		fail("Length mismatch between sequence and PSIPRED file!");
	if (fgets(buf, 160, ssfp) && !strncmp(buf, "Conf:", 5))
	    if (getpsipredh(ssfp) != seq2len)
		fail("Length mismatch between sequence and PSIPRED file!");
	fclose(ssfp);
	if (!ssstflg)
	    fail("Cannot parse PSIPRED file!");
    }

    fclose(ifp);

    if (seq2len > MAXSEQLEN - 1)
	fail("Sequence too long!");

    if (seq2len < 1)
	fail("No sequence could be extracted! (the input file may not be correctly formatted)");

    if (verbflg)
	printf("\n%d residues read:\n%s\n\n", seq2len, seq2);

    for (i = 0; i < seq2len; i++)
    {
	aa = aanum(seq2[i]);
	switch (aa)
	{
	case 20:
	    /* ASP is more common than ASN! */
	    seq2[i] = ASP;
	    break;
	case 21:
	    /* GLU is more common than GLN! */
	    seq2[i] = GLU;
	    break;
	case 22:
	    seq2[i] = UNK;
	    break;
	default:
	    seq2[i] = aa;
	    break;
	}
    }

    for (filtered = i = 0; i < seq2len - 3; i++)
	if (seq2[i] == UNK && seq2[i + 1] == UNK && seq2[i + 2] == UNK)
	{
	    filtered = TRUE;
	    break;
	}


    ifp = fopen(argv[1], "r");
    if (!ifp)
	fail("Unable to open contacts file!");

    if (!readcontacts(ifp))
	fail("Bad contacts file!");

    fclose(ifp);

    inputmat = allocmat(seq2len, seq2len, sizeof(float), TRUE);

    eigenval2 = allocvec(seq2len, sizeof(float), FALSE);
    eigenvec2 = allocmat(seq2len, seq2len, sizeof(float), FALSE);

    for (i=0; i<seq2len; i++)
	for (j=i; j<seq2len; j++)
	    if (j-i < 2)
		inputmat[i][j] = inputmat[j][i] = 1.0;
	    else if (contmap[i][j] < 0.0)
		inputmat[i][j] = inputmat[j][i] = 0.0;
	    else
		inputmat[i][j] = inputmat[j][i] = contmap[i][j];
//		inputmat[i][j] = inputmat[j][i] = contmap[i][j] > 0.05;

    eigen_decomposition(inputmat, eigenval2, eigenvec2, seq2len);

    for (i=0; i < seq2len; i++)
	eigenval2[i] = sqrtf(fabsf(eigenval2[i]));

    freemat(inputmat);

    ofp = fopen(argv[2], "w");
    if (!ofp)
	fail("Unable to open output file!");

    pat = allocmat(MAXSEQLEN + 1, MAXSEQLEN + 1, sizeof(int), FALSE);

    if (argc > 3)
	ifp = fopen(argv[3], "r");
    else if ((ifp = fopen("psichain.lst", "r")) == NULL && getenv("THREAD_DIR"))
    {
	strcpy(templname, getenv("THREAD_DIR"));
	if (templname[strlen(templname) - 1] != '/')
	    strcat(templname, "/");
	strcat(templname, "psichain.lst");
	ifp = fopen(templname, "r");
    }

    if (!ifp)
	fail("Cannot open protein list (default: psichain.lst)!");

    while (!feof(ifp))
    {
        if (!fgets(buf, 256, ifp))
            break;
        cp = buf+strlen(buf);
        while (--cp > buf)
            if (*cp == ' ')
                break;
	if (sscanf(cp, "%s", brkid) != 1)
	    break;
	if (brkid[0] == '#')
	    continue;
	if (verbflg)
	    printf("Generating template for : %s ...\n", brkid);
	/* Read coords from PDB file */
	if (getenv("TDB_DIR"))
	{
	    strcpy(tdbname, getenv("TDB_DIR"));
	    if (tdbname[strlen(tdbname) - 1] != '/')
		strcat(tdbname, "/");
	}
	else if (getenv("THREAD_DIR"))
	{
	    strcpy(tdbname, getenv("THREAD_DIR"));
	    if (tdbname[strlen(tdbname) - 1] != '/')
		strcat(tdbname, "/");
	}
	else
#ifdef DTJCOPY
	    strcpy(tdbname, "/ssd/domtdb/");
#else
	    strcpy(tdbname, "domtdb/");
#endif

	strcpy(eigname, tdbname);

	strcat(tdbname, brkid);
	strcat(tdbname, ".tdb");

	strcat(eigname, brkid);
	strcat(eigname, ".eig");

	tfp = fopen(tdbname, "r");
	if (tfp == NULL)
	{
	    fprintf(stderr, "*** main: Cannot open TDB file (%s)!\n", tdbname);
	    fflush(stderr);
	    continue;
	}

	maketemplate(tfp, TRUE);

	fclose(tfp);

	if (!seq1len)
	{
	    fprintf(stderr, "*** main: %s - length error!\n", brkid);
	    fflush(stderr);
	    continue;
	}

#if 0
	ext_sst();
#endif

	loc_sst();

	eigenvec1 = allocmat(seq1len, seq1len, sizeof(float), FALSE);
	eigenval1 = allocvec(seq1len, sizeof(float), FALSE);

	efp = fopen(eigname, "rb");
	if (efp == NULL)
	{
	    fprintf(stderr, "*** main: Cannot open eigen file (%s)!\n", eigname);
	    fflush(stderr);
	    continue;
	}

	if (fread(eigenval1, sizeof(float), seq1len, efp) != seq1len)
	    fail("EOF while reading eigenvalues from file!");

	if (fread(eigenvec1[0], seq1len * sizeof(float), seq1len, efp) != seq1len)
	    fail("EOF while reading eigenvectors from file!");

	fclose(efp);

	for (i=0; i < seq1len; i++)
	    eigenval1[i] = sqrtf(fabsf(eigenval1[i]));

	pmat = allocmat(seq2len, seq1len, sizeof(int), TRUE);

	if (prtalnflg)
	{
	    printf("\n\n>>> Alignment with %s:\n", brkid);
	}

	seqfit(&hits[nn].contactsc, &hits[nn].perco_a, &hits[nn].perco_b, &tplt);

	if (hits[nn].contactsc >= SCUTOFF)
	{
	    if (!(hits[nn].brkid = strdup(brkid)))
		fail("Out of memory!");

	    if (mod3flg)
	    {
		printf(">P1;%c%c%c%c\nstructureX:p%c%c%c%c:%s:%c:%s:%c:PDB::0.00:0.00\n", brkid[0], brkid[1], brkid[2], brkid[3], brkid[0], brkid[1], brkid[2], brkid[3], resid[0], brkid[4] == '0' ? ' ' : brkid[4], resid[seq1len - 1], brkid[4] == '0' ? ' ' : brkid[4]);
		for (i = 1; i <= tplt.length; i++)
		    if (tplt.posn_a[i] <= 0)
			putchar('-');
		    else if (seq1[tplt.posn_a[i] - 1] < GAP)
			putchar(rescodes[seq1[tplt.posn_a[i] - 1]]);
		    else
			putchar('X');
		puts("*");
		puts(">P1;SEQ\nsequence:SEQ:::::SEQ::0.00:0.00");
		for (i = 1; i <= tplt.length; i++)
		    if (tplt.posn_b[i] <= 0)
			putchar('-');
		    else if (seq2[tplt.posn_b[i] - 1] < GAP)
			putchar(rescodes[seq2[tplt.posn_b[i] - 1]]);
		    else
			putchar('X');
		puts("*");
	    }

	    if (motifflg)
	    {
		strcpy(wpname, motpref);
		strcat(wpname, "_");
		strcat(wpname, brkid);
		strcat(wpname, ".model.pdb");
		writepdb(wpname, brkid);
	    }

	    nn++;
	}

	fflush(stdout);

	freemat(pmat);
	freemat(distmat);
	freevec(eigenval1);
	freemat(eigenvec1);
    }

    if (ofp != NULL)
	for (i = 0; i < nn; i++)
	{
	    fprintf(ofp, "%7.3f %5.1f %5.1f %7.7s\n",
		    hits[i].contactsc,
		    hits[i].perco_a,
		    hits[i].perco_b,
		    hits[i].brkid);

	    fflush(ofp);
	}

    return 0;
}
