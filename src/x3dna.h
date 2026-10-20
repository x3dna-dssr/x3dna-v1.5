#ifndef _X3DNA_H
#define _X3DNA_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#define SFACTOR    2L  /* scale factor for increasing memory */
#define PI 3.141592653589793
#define XEPS 1.0e-7
#define XBIG 1.0e+18
#define MFACTOR 10000.0
#define X3DNA_VER "3DNA (v1.6, Feb. 2006) by Dr. Xiang-Jun Lu"
#define NMISC 34  /* number of characters of miscellaneous items */
#define NATOMCOL 11  /* number of atom with defined colors */
#define NBASECOL 7  /* number of base with defined colors */

#define PS_DFTSIZE 500  /* 500 points PS default size */
#define PS_BOUND 10  /* boundary offset */
#define FIG_DFTSIZE 8333  /* 8333 units XFIG default size */
#define FIG_BOUND 166  /* boundary offset: 10/72*1200 */

#define PAR_FILE "misc_3dna.par"  /* miscellaneous parameters */
#define BASE_FILE "baselist.dat"  /* 3-letter to 1-letter base residue */
#define HELP3DNA "help3dna.dat"  /* help file for 3DNA */

#define REF_FILE "ref_frames.dat"  /* reference frames */
#define MREF_FILE "mref_frames.dat"  /* multiplet reference frames in <find_pair> */
#define POC_FILE "poc_haxis.r3d"  /* P, O4' and C1' atom radius in r3d format */
#define MUL_FILE "multiplets.pdb"  /* multiplets */
#define ALLP_FILE "allpairs.pdb"  /* all base-pairs */
#define BESTP_FILE "bestpairs.pdb"  /* best base-pairs */
#define STACK_FILE "stacking.pdb"  /* for stacking diagram */
#define HSTACK_FILE "hstacking.pdb"  /* for stacking w.r.t. to middle helical frame */
#define HLXREG_FILE "hel_regions.pdb"  /* helical regions in <find_pair> */
#define BPORDER_FILE "bp_order.dat"  /* base-pair ordering in <find_pair> */
#define COLCHN_FILE "col_chains.scr"  /* color chains in <find_pair> */
#define COLHLX_FILE "col_helices.scr"  /* color helices in <find_pair> */
#define AUX_FILE "auxiliary.par"  /* auxiliary parameters in <analyze> */
#define BPSTEP_FILE "bp_step.par"  /* base-pair step parameters in <analyze> */
#define HLXSTEP_FILE "bp_helical.par"  /* helical step parameters in <analyze> */
#define SEVEN_FILE "cf_7methods.par"  /* compare seven methods in <analyze> */
#define HB_FILE "hbonds_info.dat"  /* H-bonding information <r3d_atom/stack2img> */
#define LKG_FILE "bonds_lkg.dat"  /* bond linkage, as in HB_FILE */
#define SNUM_FILE "serial_num.pdb"  /* PDB file with serial atom numbers */
#define ATOMALC_FILE "atom_lkg.alc"  /* atomic linkage <r3d_atom/stack2img> */
#define BBLKALC_FILE "bblk_lkg.alc"  /* base block linkage <pdb2img> */
#define PBLKALC_FILE "pblk_lkg.alc"  /* peptide block linkage <pdb2img> */
#define TMP_FILE "tmp_file"  /* temporary multiplets input file <find_pair> */
#define MULBP_FILE "mulbp.inp"  /* multiplets input file <find_pair> */
#define ROTMAT_FILE "rotmat.dat"  /* rotation matrix <rotate_mol> */
#define VIEW1_FILE "pmiview1"  /* PMI view #1 file <rotate_mol> */
#define VIEW2_FILE "pmiview2"  /* PMI view #2 file <rotate_mol> */
#define VIEW3_FILE "pmiview3"  /* PMI view #3 file <rotate_mol> */

#define NP 101L  /* maximum number of pairs per base */
#define BOND_UPPER_LIMIT 2.5  /* for function torsion */
#define HTWIST0 0.05  /* minimum helical twist */
#define BOND_FACTOR 1.15  /* bond distance criterion */
#define NBOND_FNUM  2.0  /* estimated # of bond from # of atoms */
#define NON_WC_IDX  6  /* non-Watson-Crick base index */
#define AXIS_LENGTH 3.5  /* reference axis length */
#define BUF512 512
#define NUM_BASE_ATOMS 50  /* max. no. of base atoms in a residue */
#define NUM_RESIDUE_ATOMS 100  /* max. no. of atoms in a residue */
#define NUM_DINUCLEOTIDE_ATOMS 400  /* max. no. of atoms per dinucleotide */
#define EMPTY_NUMBER -9999.99
#define EMPTY_CRITERION -9999
#define MAXBASE 30000  /* maximum number of bases in regular/fiber */
#define NELE 12  /* 12 elements */
#define WC_SHEAR 2.0  /* maximum SHEAR to be a WC pair: normal distance: ~2.3 A */
#define O3P_UPPER 2.5  /* upper limit for O3-P connection */
#define RTNNUM 37  /* number of returned values from check_pair */
#define PSTNUM 29  /* number of parameters kept in pair_stat */
#define OLCRT 1.2  /* criterion for overlaps */
#define MBASES 50  /* maximum bases per unit */
#define MAXCH 100  /* maximum # of chains in a biological unit */
#define END_STACK_XANG 125.0  /* find_pair; bdl070 (105 deg) */
#define MAXCLEN 52  /* Maximum length of a color name */

#define RA_LIST " C4 ", " N3 ", " C2 ", " N1 ", " C6 ", " C5 ", " N7 ", " C8 ", " N9 "
#define WC_LIST "XX", "AT", "AU", "TA", "UA", "GC", "IC", "CG", "CI"
#define CX_LIST "ACGITUX"
#define CB_LIST "ACGITU"
#define WATER_LIST {"H2O", "HHO", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP"}

/* for calculating overlap area between polygons */
#define MNPOLY 1000  /* maximum size of a polygon */

#define UNUSED_PARAMETER(x) (void)(x)  /* to prevent warning about unused parameters */

typedef struct {
    double x;
    double y;
} point;

typedef struct {
    point ip;
    point rx;
    point ry;
    long inside;
} vertex;

#include "x3dna_header.h"

#include "nrutil.h"

#endif  /* _X3DNA_H */
