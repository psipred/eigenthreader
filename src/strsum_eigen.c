#define BREAKS
#define DISCONTINUOUS
/****************************************************************************
 *
 * EIGEN CONTACT THREADER TDB FILE GENERATOR version 1.0 (Sep 2015)
 * written by
 * David T. Jones,
 * University College,
 * Gower Street,
 * London WC1E 6BT
 *
 ****************************************************************************/

/* This program creates the database files for THREADER */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <string.h>

#define FALSE 0
#define TRUE 1

#define BIG 100000000
#define VBIG 1e32F

/* Max-size-of-problem parameters */
#define MAXSEQLEN 2011
#define MAXWINLEN 250
#define MAXATOMS 30000

#define SQR(x) ((x)*(x))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define ran0() ((rand()&32767)/32768.0)
#define CH malloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);
#define vecprod(a,b,c) (a[0]=b[1]*c[2]-b[2]*c[1],a[1]=b[2]*c[0]-b[0]*c[2],a[2]=b[0]*c[1]-b[1]*c[0])
#define dotprod(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/* Constants for calculating CB coords (Prot. Eng. Vol.2 p.121) */
#define	TETH_ANG 0.9128
#define CACBDIST 1.538

const float     ZERO = 0.0F, ONE = 1.0F, TWO = 2.0F, THREE = 3.0F;

char            seq[MAXSEQLEN], ic[MAXSEQLEN];
char            modseq[MAXSEQLEN], tsstruc[MAXSEQLEN];
short           modsdx[MAXSEQLEN], modsdx2[MAXSEQLEN], gaps[MAXSEQLEN], ri[MAXSEQLEN];
short           sstidx[MAXSEQLEN], tooi[MAXSEQLEN];
int             maxgplen[MAXSEQLEN];
int             modtot[MAXSEQLEN], firstid, seqlen, alnlen, nseqs;

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

float           e_contrib[100];

struct HIT
{
    char            brkid[7];
    float           epair, esolv;
}
hits[500];

typedef struct
{
    short           length;
    short           position_a[MAXSEQLEN], position_b[MAXSEQLEN];
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

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    GAP, UNK
};

const char     *ssnames[] =
{
    "COIL", "HELIX", "STRAND", "A-STRAND", "P-STRAND", "TURN"
};

enum ssids
{
    COIL, HELIX, STRAND, ASTRAND, PSTRAND, TURN
};

const char     *atmnames[] =
{
    "CA", "CB", "O ", "N ", "C "
};

enum atmcodes
{
    CAATOM, CBATOM, OATOM, NATOM, CATOM
};

const char     *rescodes = "ARNDCQEGHILKMFPSTWYV-X";
const char     *sscodes = "CHEAPT";

/* DSSP residue solvent accessible area in GGXGG extended pentapeptide */
const float     resacc[22] =
{
    113.0, 253.0, 167.0, 167.0, 140.0, 199.0, 198.0, 88.0, 194.0, 178.0,
    179.0, 215.0, 194.0, 226.0, 151.0, 134.0, 148.0, 268.0, 242.0, 157.0,
    88.0, 88.0
};

/*  BLOSUM 62 */
short           aamat[23][23] = {
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

/* List of atom pairs to evaluate in final model evaluation stage */
#define NPAIRS 5

const short     atompair[][2] =
{
    {CBATOM, CBATOM},
    {CBATOM, NATOM},
    {CBATOM, OATOM},
    {NATOM, CBATOM},
    {OATOM, CBATOM},
};

/* Structure generation arrays */
float           (**distmat)[4][4];
float           x[MAXSEQLEN][5], y[MAXSEQLEN][5], z[MAXSEQLEN][5];
short           trelacc[MAXSEQLEN];
float           tmat[MAXSEQLEN][22];
float           phi[MAXSEQLEN], psi[MAXSEQLEN], omega[MAXSEQLEN];

/* Dump a rude message to standard error and exit */
void
                fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
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
	    if (run >= 4 || (run == 3 && (i && struc[i - 1] == 'H' || struc[i + run] == 'H')))
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
    int             i, istart, l;

    for (i = 1; i < seqlen - 1; i++)
	if (tsstruc[i - 1] == tsstruc[i + 1])
	    tsstruc[i] = tsstruc[i - 1];
}

/* Extend DSSP secondary structure definitions */
void
                ext_sst(void)
{
    int             i, istart, l;

    for (i = 0; i < seqlen; i++)
	if (i > 0 && !tsstruc[i] && tsstruc[i - 1])
	{
	    tsstruc[i] = tsstruc[i - 1];
	    i++;
	}
	else if (i < seqlen - 1 && !tsstruc[i] && tsstruc[i + 1])
	    tsstruc[i] = tsstruc[i + 1];
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

float           torsion(const float *a, const float *b, const float *c, const float *d)
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
    if ((len_a_b == 0.0) || (len_b_c == 0.0) || (len_c_d == 0))
	return (-999.0);

    /* Calculate dot products to form cosines */
    ab_bc = dotprod(a_b, b_c) / (len_a_b * len_b_c);
    ab_cd = dotprod(a_b, c_d) / (len_a_b * len_c_d);
    bc_cd = dotprod(b_c, c_d) / (len_b_c * len_c_d);

    /* Form sines */
    s_ab = sqrt(1.0 - SQR(ab_bc));
    s_bc = sqrt(1.0 - SQR(bc_cd));

    costor = (ab_bc * bc_cd - ab_cd) / (s_ab * s_bc);
    costsq = SQR(costor);
    if ((costsq >= 1.0) && (costor < 0.0))
	tor = 180.0;

    /* If the angle is not == 180 degs calculate sign using sine */
    if (costsq < 1.0)
    {
	sintor = sqrt(1.0 - costsq);
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
	if (sign < 0.0)
	    tor = -tor;
    }

    return (tor);
}


/* Summarize PDB/DSSP file */
int
                summary(FILE * pfp, FILE * dfp, FILE * ofp)
{
    char            buf[160], compnd[160], *p, whichatm[MAXSEQLEN], inscode;
    int             acc, atc, i, j, k, nres, resid, lastid, natoms;
    float           f, dv, cca[3], nca[3], xx[3], yy[3], sx, sy;
    float           a[3], b[3], c[3], d[3], last_nx, last_ny, last_nz;
    short           at1, at2;

    seqlen = 0;
    for (i = 0; i < MAXSEQLEN; i++)
    {
	whichatm[i] = 0;
	tsstruc[i] = ic[i] = ' ';
	ri[i] = 9999;
	hashtbl[i].sstruc = '\0';
	hashtbl[i].acc = 0;
    }
    while (!feof(dfp))
    {
	if (!fgets(buf, 160, dfp))
	{
	    puts("Bad DSSP file!");
	    return;
	}
	if (strstr(buf, "NUMBER OF CHAINS"))
	{
	    if (sscanf(buf, "%d%d%d%d%d", &nres, &i, &i, &i, &i) != 5)
		fail("Bad DSSP file - check number of residues/chains!");
	}
	if (buf[2] == '#')
	    break;
    }
    while (!feof(dfp))
    {
	if (!fgets(buf, 160, dfp))
	    break;
	inscode = buf[10];
	buf[10] = ' ';
	sscanf(buf + 6, "%d", &resid);
	if (++seqlen > MAXSEQLEN)
	    fail("Too many residues!");
	i = (resid + resid) % MAXSEQLEN;
	while (hashtbl[i].sstruc != '\0')
	    i = (i + 1) % MAXSEQLEN;
	hashtbl[i].resid = resid;
	hashtbl[i].inscode = inscode;
	if (buf[16] == 'E' && (isalpha(buf[23]) || isalpha(buf[24])))
	    hashtbl[i].sstruc = isupper(buf[23]) || isupper(buf[24]) ? 'A' : 'P';
	else
	    hashtbl[i].sstruc = buf[16];
	sscanf(buf + 35, "%d", &acc);
	hashtbl[i].acc = acc;
    }

    i = -1;
    natoms = 0;
    compnd[0] = '\0';
    while (!feof(pfp))
    {
	if (fgets(buf, 160, pfp) == NULL)
	    break;
	if (!strncmp(buf, "ENDMDL", 6))
	    break;
	/* printf("%d %s\n",i,buf); */
	if (compnd[0] == '\0' && !strncmp(buf, "COMPND", 6) && !strstr(buf, "MOL_ID:"))
	{
	    buf[strlen(buf)-1] = '\0';
	    strcpy(compnd, buf+10);
	    compnd[72] = '\0';
	    p = strstr(compnd, "  ");
	    if (p)
		*p = '\0';
	}
	if (strncmp(buf, "ATOM", 4) || buf[16] != ' ' && buf[16] != 'A')
	    continue;
	if (!strncmp(buf + 13, "OT1", 3))
	    memcpy(buf+13, "O  ", 3);
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
		    else if (inscode == ' ' && resid - lastid > 1 && SQR(x[i][NATOM] - x[i - 1][CATOM]) + SQR(y[i][NATOM] - y[i - 1][CATOM]) + SQR(z[i][NATOM] - z[i - 1][CATOM]) > 6.25F)
		    {
#ifdef NOBREAKS
			puts("SKIPPED: Chain break - file ignored!");
			seqlen = 0;
			return;
#endif
#ifdef DISCONTINUOUS
			printf("WARNING: %d residues possibly missing (inserted domain?)!\n", resid - lastid - 1);
#else
			printf("SKIPPED: %d residues possibly missing (inserted domain?)!\n", resid - lastid - 1);
			seqlen = 0;
			return;
#endif
		    }
		    else if (inscode != ' ' && SQR(x[i][NATOM] - x[i - 1][NATOM]) + SQR(y[i][NATOM] - y[i - 1][NATOM]) + SQR(z[i][NATOM] - z[i - 1][NATOM]) < 1.0F)
		    {
			puts("WARNING: probable multiple occupancy!");
			x[i - 1][NATOM] = x[i][NATOM];
			y[i - 1][NATOM] = y[i][NATOM];
			z[i - 1][NATOM] = z[i][NATOM];
			i--;
		    }
		    lastid = resid;
		    for (k = 0; k < 20; k++)
			if (!strncmp(rnames[k], buf + 17, 3))
			    break;
		    seq[i] = k;
		    ri[i] = resid;
		    ic[i] = inscode;
		    k = (resid + resid) % MAXSEQLEN;
		    while ((hashtbl[k].resid != resid || hashtbl[k].inscode != inscode) && hashtbl[k].sstruc)
			k = (k + 1) % MAXSEQLEN;
		    if (hashtbl[k].sstruc == '\0')
		    {
#ifdef NOBREAKS
			puts("Missing DSSP entries - chain ignored!");
			seqlen = 0;
			return;
#endif
			printf("WARNING: No DSSP entry found for %d%c!\n", resid, inscode);
			tsstruc[i] = ' ';
			trelacc[i] = 0;
		    }
		    else
		    {
			if (seq[i] <= VAL)
			    trelacc[i] = (short) (100.0F * hashtbl[k].acc / resacc[seq[i]]);
			else
			    trelacc[i] = (short) (100.0F * hashtbl[k].acc / 170.0F);
			tsstruc[i] = hashtbl[k].sstruc;
		    }
		}
		else if (i >= 0)
		    readxyz(buf + 30, &x[i][atc], &y[i][atc], &z[i][atc]);
		whichatm[i] |= 1 << atc;
	    }
    }

    seqlen = nres = i + 1;

    /* Check atoms */
    for (i = 0; i < nres; i++)
	if (seq[i] != GAP)
	{
	    if (!(whichatm[i] & (1 << NATOM)))
	    {
		printf("FATAL: Unexpected missing N atom in %d!\n", i + 1);
		exit(1);
	    }
	    if (!(whichatm[i] & (1 << CBATOM)))
	    {
		if (!(whichatm[i] & (1 << CAATOM)) || !(whichatm[i] & (1 << CATOM)) || !(whichatm[i] & (1 << NATOM)) || !(whichatm[i] & (1 << OATOM)))
		{
		    /* Not much left of this residue! */
		    printf("WARNING: Missing main chain atom in %d!\n", i + 1);
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
		if (seq[i] != GLY)
		    printf("WARNING: dummy CB atom constructed for %d (%s)!\n", i + 1, rnames[seq[i]]);
	    }
	}

    if (!nres)
	return;

    for (i = 0; i < nres; i++)
    {
	if (seq[i] != GAP && i < nres - 1 && seq[i + 1] != GAP)
	{
	    a[0] = x[i][NATOM];
	    a[1] = y[i][NATOM];
	    a[2] = z[i][NATOM];
	    b[0] = x[i][CAATOM];
	    b[1] = y[i][CAATOM];
	    b[2] = z[i][CAATOM];
	    c[0] = x[i][CATOM];
	    c[1] = y[i][CATOM];
	    c[2] = z[i][CATOM];
	    d[0] = x[i + 1][NATOM];
	    d[1] = y[i + 1][NATOM];
	    d[2] = z[i + 1][NATOM];
	    psi[i] = torsion(a, b, c, d);
	}
	else
	    psi[i] = 999.999;

	if (seq[i] != GAP && i < nres - 1 && seq[i + 1] != GAP)
	{
	    a[0] = x[i][CAATOM];
	    a[1] = y[i][CAATOM];
	    a[2] = z[i][CAATOM];
	    b[0] = x[i][CATOM];
	    b[1] = y[i][CATOM];
	    b[2] = z[i][CATOM];
	    c[0] = x[i + 1][NATOM];
	    c[1] = y[i + 1][NATOM];
	    c[2] = z[i + 1][NATOM];
	    d[0] = x[i + 1][CAATOM];
	    d[1] = y[i + 1][CAATOM];
	    d[2] = z[i + 1][CAATOM];
	    omega[i] = torsion(a, b, c, d);
	}
	else
	    omega[i] = 999.999;

	if (seq[i] != GAP && i && seq[i - 1] != GAP)
	{
	    a[0] = x[i - 1][CATOM];
	    a[1] = y[i - 1][CATOM];
	    a[2] = z[i - 1][CATOM];
	    b[0] = x[i][NATOM];
	    b[1] = y[i][NATOM];
	    b[2] = z[i][NATOM];
	    c[0] = x[i][CAATOM];
	    c[1] = y[i][CAATOM];
	    c[2] = z[i][CAATOM];
	    d[0] = x[i][CATOM];
	    d[1] = y[i][CATOM];
	    d[2] = z[i][CATOM];
	    phi[i] = torsion(a, b, c, d);
	}
	else
	    phi[i] = 999.999;
    }

    fprintf(ofp, "#TDB ------");
    fprintf(ofp, " %5d %s\n", nres, compnd);

    for (i = 0; i < nres; i++)
	if (seq[i] < GAP)
    {
	fprintf(ofp, "%4d %c %c %3d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %4d%c ",
		i + 1, rescodes[seq[i]], tsstruc[i], trelacc[i],
		phi[i], psi[i], omega[i],
		x[i][NATOM], y[i][NATOM], z[i][NATOM],
		x[i][CAATOM], y[i][CAATOM], z[i][CAATOM],
		x[i][CATOM], y[i][CATOM], z[i][CATOM],
		x[i][OATOM], y[i][OATOM], z[i][OATOM],
		x[i][CBATOM], y[i][CBATOM], z[i][CBATOM], ri[i], ic[i]);

	for (j = 0; j < 20; j++)
	    fprintf(ofp, "%6d", 100 * aamat[seq[i]][j]);
	fprintf(ofp, "\n");
    }
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
    int             i;
    
    free(*(void **)p);
    free(p);
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

/* Convert string to lower case */
void
                lc_str(char *str)
{
    int             i;

    for (i = strlen(str) - 1; i >= 0; i--)
	if (isupper(str[i]))
	    str[i] = tolower(str[i]);
}

int             main(int argc, char **argv)
{
    int             a, b, aa, i, j, k, l, nn = 0, atpair, t, nres;
    char           *cmdstr, pred;
    char            brkid[10], pdbcode[5], chainid, desc[512], oname[80], *dssppdbn;
    float           **cmat, **eigenvec, *eigenval;
    FILE           *ifp, *ofp, *pfp, *dfp;
    char            pdbname[512], dsspname[512];

    if (argc != 5)
        fail("Usage: strsum pdbfile dsspfile tdbfile eigenfile");

    strcpy(pdbname, argv[1]);

    pfp = fopen(pdbname, "r");
    if (pfp == NULL)
    {
	printf("!!!MISSING %s\n", brkid);
	fprintf(stderr, "*** main: Cannot open PDB file!\n");
	exit(1);
    }

    strcpy(dsspname, argv[2]);
    dfp = fopen(dsspname, "r");
    if (dfp == NULL)
    {
	fclose(pfp);
	fprintf(stderr, "*** main: Cannot open DSSP file!\n");
	exit(1);
    }

    ofp = fopen(argv[3], "w");

    if (!ofp)
	fail("Cannot open TDB file for output!");

    nres = summary(pfp, dfp, ofp);
    
    fclose(dfp);
    fclose(pfp);
    fclose(ofp);
    
    ofp = fopen(argv[4], "wb");
    
    if (!ofp)
	fail("Cannot open eigen file for output!");
    
    cmat = allocmat(nres, nres, sizeof(float), FALSE);
    eigenvec = allocmat(nres, nres, sizeof(float), FALSE);
    eigenval = allocvec(nres, sizeof(float), FALSE);
    
    for (i = 0; i < nres; i++)
	for (j = i; j < nres; j++)
	    cmat[i][j] = cmat[j][i] = (SQR(x[i][CBATOM] - x[j][CBATOM]) + SQR(y[i][CBATOM] - y[j][CBATOM]) + SQR(z[i][CBATOM] - z[j][CBATOM]) < 64.0F);
    
    eigen_decomposition(cmat, eigenval, eigenvec, nres);

    if (fwrite(eigenval, sizeof(float), nres, ofp) != nres)
	fail("Cannot write eigenvalue data!");

    if (fwrite(eigenvec[0], nres * sizeof(float), nres, ofp) != nres)
	fail("Cannot write eigenvector data!");

    fclose(ofp);
    
    return 0;
}
