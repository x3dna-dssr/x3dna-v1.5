#include "x3dna.h"

int main(int argc, char *argv[])
/* utility for generating 43 + 2 + 2 + 3 fiber DNA/RNA structures */
{
    char FDIR[BUF512], pdbfile[BUF512];
    long ABC[] = { -1, 1, 4, 7 };
    long S1S2[] = { -1, 2, 3, 5, 6, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20,
        21, 22, 23, 24, 25, 26, 27, 28, 29, 37, 38, 39, 40, 41, 43,
        48, 49, 50, 51, 52
    };
    long TRI[] = { -1, 30, 31, 32, 33, 34, 42 };
    long str_id = 4;  /* defaulted to B-DNA */
    double rise, twist;

    fiber_dir(FDIR);
    fiber_cmdline(argc, argv, FDIR, &str_id, pdbfile);
    get_twist_rise(str_id, &twist, &rise);

    ABC[0] = sizeof ABC / sizeof ABC[0] - 1;
    S1S2[0] = sizeof S1S2 / sizeof S1S2[0] - 1;
    TRI[0] = sizeof TRI / sizeof TRI[0] - 1;

    /* Five types of fiber structure generating routines
       % o single residue as repeating unit, including 1, 4, 7, 8 & 12,
       %                                      (44,45), (46,47), (53-55)
       % o base-pair(s) as repeating unit, including 2, 3, 5, 6, 9-11,
       %                                         13-29, 37-41 & 43, 48-52
       % o triplex, including 30-34, 42
       % o quadruplex, including 35
       % o single helix, including 36 */
    if (lval_in_set(str_id, 1, ABC[0], ABC))
        residue(twist, rise, "", 1, FDIR, pdbfile);
    else if (str_id == 8)  /* C-DNA  poly d(GGT) : poly d(ACC) */
        residue(twist, rise, "GGT", 1, FDIR, pdbfile);
    else if (str_id == 12)  /* poly d(AAT) : poly d(ATT) */
        residue(twist, rise, "AAT", 1, FDIR, pdbfile);
    else if (lval_in_range(str_id, 44, 45))  /* poly d(A) : poly d(T) */
        residue(twist, rise, "A", 0, FDIR, pdbfile);  /* T in stnd 2: positive z-axis */
    else if (lval_in_range(str_id, 46, 47) ||  /* BI/BII-DNA: Levitt JMB 2000 */
             lval_in_range(str_id, 53, 55))  /* Premilat & Albiser A-, B- and C-DNA */
        residue(twist, rise, "", -1, FDIR, pdbfile);  /* positive z-axis */
    else if (lval_in_set(str_id, 1, S1S2[0], S1S2))
        multi_bp(str_id, twist, rise, FDIR, pdbfile);
    else if (lval_in_set(str_id, 1, TRI[0], TRI))
        triplex(str_id, twist, rise, FDIR, pdbfile);
    else if (str_id == 35)
        quadruplex(twist, rise, FDIR, pdbfile);
    else if (str_id == 36)
        shelix(twist, rise, FDIR, pdbfile);
    else
        fatal("fiber structure number [%ld] out of 1-55 range\n", str_id);

    return 0;
}

void fiber_dir(char *FDIR)
/* Get the directory containing FIBER models: cf. get_BDIR */
{
    char *temp;

    if ((temp = getenv("X3DNA")) != NULL) {
        strcpy(FDIR, temp);
        check_slash(FDIR);
        strcat(FDIR, "FIBER/");
    } else if ((temp = getenv("HOME")) != NULL ||  /* Unix */
               (temp = getenv("HOMEDRIVE")) != NULL) {  /* PC */
        strcpy(FDIR, temp);
        check_slash(FDIR);
        strcat(FDIR, "X3DNA/FIBER/");
    } else
        fatal("cannot locate fiber structure data files\n");
    fprintf(stderr, "Fiber data in directory: %s\n", FDIR);
}

void fiber_usage(void)
{
    help3dna("fiber");
    contact_msg(0);
}

void fiber_cmdline(int argc, char *argv[], char *FDIR, long *str_id, char *pdbfile)
{
    long i, j;

    if (!lval_in_range(argc, 2, 3))
        fiber_usage();
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        if (sscanf(argv[i], "-%ld", str_id) == 1) {
            if (!lval_in_range(*str_id, 1, 55)) {
                fprintf(stderr, "Structure id# out of range (1-55)\n");
                fiber_usage();
            } else
                continue;
        }
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'M' || argv[i][j] == 'L')
                list_str_info(FDIR);
            else if (argv[i][j] == 'A')
                *str_id = 1;
            else if (argv[i][j] == 'B')
                *str_id = 4;
            else if (argv[i][j] == 'C')
                *str_id = 47;
            else if (argv[i][j] == 'D')
                *str_id = 48;
            else if (argv[i][j] == 'Z')
                *str_id = 15;
            else
                fiber_usage();
    }
    if (argc == i + 1) {
        strcpy(pdbfile, argv[i]);
        sprintf(FDIR, "%sStr%2.2ld/", FDIR, *str_id);
    } else
        fiber_usage();
}

void list_str_info(char *FDIR)
/* brief description of the 55 fiber structures */
{
    char README_FILE[BUF512], str[BUF512];
    FILE *fp;

    sprintf(README_FILE, "%sREADME", FDIR);
    fp = open_file(README_FILE, "r");
    while (fgets(str, sizeof str, fp) != NULL)
        if (!strncmp(str, "id# ", 4)) {
            fprintf(stderr, "\n%s", str);
            while (fgets(str, sizeof str, fp) != NULL)
                fprintf(stderr, "%s", str);
        }
    close_file(fp);
    contact_msg(0);
}

void get_twist_rise(long str_id, double *twist, double *rise)
{
    long idx;
    double twist_rise[55][3] = {
        {1, 32.7, 2.548},
        {2, 65.5, 5.095},
        {3, 0.0, 28.030},
        {4, 36.0, 3.375},
        {5, 72.0, 6.720},
        {6, 180.0, 16.864},
        {7, 38.6, 3.310},
        {8, 40.0, 3.312},
        {9, 120.0, 9.937},
        {10, 80.0, 6.467},
        {11, 80.0, 6.467},
        {12, 45.0, 3.013},
        {13, 90.0, 6.125},
        {14, -90.0, 18.500},
        {15, -60.0, 7.250},
        {16, -51.4, 7.571},
        {17, 0.0, 10.200},
        {18, 36.0, 3.230},
        {19, 36.0, 3.233},
        {20, 32.7, 2.812},
        {21, 30.0, 3.000},
        {22, 32.7, 2.560},
        {23, 32.0, 2.780},
        {24, 36.0, 3.130},
        {25, 32.7, 3.060},
        {26, 36.0, 3.010},
        {27, 32.7, 2.518},
        {28, 32.7, 2.596},
        {29, 32.7, 2.596},
        {30, 32.7, 3.160},
        {31, 30.0, 3.260},
        {32, 32.7, 3.040},
        {33, 30.0, 3.040},
        {34, 30.0, 3.290},
        {35, 31.3, 3.410},
        {36, 60.0, 3.155},
        {37, 36.0, 3.200},
        {38, 36.0, 3.240},
        {39, 72.0, 6.480},
        {40, 72.0, 6.460},
        {41, 144.0, 13.540},
        {42, 32.7, 3.040},
        {43, 36.0, 3.200},
        {44, 36.0, 3.233},
        {45, 36.0, 3.233},
        {46, 36.0, 3.380},
        {47, 40.0, 3.320},
        {48, 87.8, 6.020},
        {49, 60.0, 7.200},
        {50, 60.0, 7.200},
        {51, 31.6, 3.220},
        {52, 90.0, 6.060},
        {53, -38.7, 3.290},
        {54, 32.73, 2.560},
        {55, 36.0, 3.39}
    };

    idx = str_id - 1;
    *twist = twist_rise[idx][1];
    *rise = twist_rise[idx][2];

    fprintf(stderr, "Structure #%ld; ", str_id);
    fprintf(stderr, "Twist: %.1f (degrees); Rise: %.3f (Angstrom)\n", *twist, *rise);
}

void residue(double twist, double rise, char *b_rep, long s2r, char *FDIR, char *pdbfile)
/* repeating unit is each residue: str# 1, 4, 7, 8 & 12 + (44-45) + (46-47) + (53-55) */
{
    char fname[BUF512];
    static char *Wbase = "ATGC", *Cbase = "TACG";
    char *sChainID, *tChainID;
    char **sAtomName, **sResName;
    char **tAtomName, **tAtomName2, **tResName;
    char **bp_seq;

    double dphi, dz, phi;
    double **sxyz, **txyz, **txyz2;

    long num_atoms = 0, n_rep;
    long i, ik, j, num1, num2, tnum = 0, tnum2 = 0;
    long idx[4];
    long *sResSeq, *tResSeq;
    long **s2idx;  /* strand II index */

    /* get the base-pair sequence */
    bp_seq = cmatrix(1, MAXBASE, 1, 2);
    if (!strcmp(b_rep, ""))  /* 1, 4 & 7 + (46-47) + (53-55) */
        bp_sequence(Wbase, Cbase, &n_rep, bp_seq);
    else  /* 8 & 12 + (44, 45) */
        bp_repeat(Wbase, Cbase, b_rep, &n_rep, bp_seq);

    /* get the total number of atoms in the final structure */
    for (i = 0; i < 4; i++)
        idx[i] = 0;
    for (i = 1; i <= n_rep; i++)
        for (j = 1; j <= 2; j++)
            ++idx[strchr(Wbase, bp_seq[i][j]) - Wbase];
    for (i = 0; i < 4; i++) {
        if (idx[i]) {
            sprintf(fname, "%s%c.rpt", FDIR, Wbase[i]);
            num_atoms += idx[i] * number_of_atoms(fname, 1, "*");
        }
    }

    tAtomName = cmatrix(1, num_atoms, 0, 4);
    tResName = cmatrix(1, num_atoms, 0, 3);
    tChainID = cvector(1, num_atoms);
    tResSeq = lvector(1, num_atoms);
    txyz = dmatrix(1, num_atoms, 1, 3);
    tAtomName2 = cmatrix(1, num_atoms, 0, 4);
    txyz2 = dmatrix(1, num_atoms, 1, 3);

    sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
    sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
    sChainID = cvector(1, NUM_RESIDUE_ATOMS);
    sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
    sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);

    s2idx = lmatrix(1, n_rep, 1, 2);

    for (i = 1; i <= n_rep; i++) {
        dphi = -(i - 1) * twist;  /* negative */
        dz = -(i - 1) * rise;  /* negative */

        if (s2r == -1) {  /* (46, 47): positive z-axis */
            dphi = -dphi;
            dz = -dz;
        }

        /* Strand I */
        sprintf(fname, "%s%c.rpt", FDIR, bp_seq[i][1]);
        num1 = read_pdb(fname, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
                        NULL, 1, "*");
        for (j = 1; j <= num1; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], sAtomName[j]);
            sprintf(tResName[ik], "  %c", bp_seq[i][1]);
            tChainID[ik] = 'A';
            tResSeq[ik] = i;

            phi = deg2rad(dphi + sxyz[j][3]);
            txyz[ik][1] = sxyz[j][2] * cos(phi);
            txyz[ik][2] = sxyz[j][2] * sin(phi);
            txyz[ik][3] = dz + sxyz[j][1];
        }

        /* Strand II */
        sprintf(fname, "%s%c.rpt", FDIR, bp_seq[i][2]);
        num2 = read_pdb(fname, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
                        NULL, 1, "*");

        for (j = 1; j <= num2; j++) {
            if (s2r) {  /* 1, 4, 7, 8 & 12 */
                sxyz[j][1] = -sxyz[j][1];  /* -z */
                sxyz[j][3] = -sxyz[j][3];  /* -phi */
            }

            ik = tnum2 + j;
            strcpy(tAtomName2[ik], sAtomName[j]);

            phi = deg2rad(dphi + sxyz[j][3]);
            txyz2[ik][1] = sxyz[j][2] * cos(phi);
            txyz2[ik][2] = sxyz[j][2] * sin(phi);
            txyz2[ik][3] = dz + sxyz[j][1];
        }

        s2idx[i][1] = tnum2;
        s2idx[i][2] = num2;

        tnum += num1;
        tnum2 += num2;
    }

    /* "reverse" strand II and combined with I */
    reverse_stnd2(n_rep, bp_seq, s2idx, &tnum, tAtomName, tResName,
                  tChainID, tResSeq, txyz, tAtomName2, txyz2, 0);

    write_pdb(tnum, tAtomName, tResName, tChainID, tResSeq, txyz, NULL, pdbfile);

    free_cmatrix(bp_seq, 1, MAXBASE, 1, 2);
    free_pdb(num_atoms, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
    free_cmatrix(tAtomName2, 1, num_atoms, 0, 4);
    free_dmatrix(txyz2, 1, num_atoms, 1, 3);
    free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL);
    free_lmatrix(s2idx, 1, n_rep, 1, 2);
}

void multi_bp(long str_id, double twist, double rise, char *FDIR, char *pdbfile)
{
    char base_rpt[BUF512], fname[BUF512];
    char *ChainID1, *ChainID2, *tChainID;
    char **AtomName1, **AtomName2, **ResName1, **ResName2, **Miscs1;
    char **tAtomName, **tAtomName2, **tResName, **tResName2;

    double dphi, dz, phi;
    double **xyz1, **xyz2, **txyz, **txyz2;

    long i, ik, j, n_rep, nr, num1, num2, tnum, tnum2;
    long *ResSeq1, *ResSeq2, *tResSeq, *tResSeq2;
    long **s2idx;  /* strand II index */

    switch (str_id) {
    case 2:
        strcpy(base_rpt, "ABr5U:ABr5U");
        break;
    case 3:
        strcpy(base_rpt, "ATCGGAATGGT:ACCATTCCGAT");
        break;
    case 5:
        strcpy(base_rpt, "CG:CG");
        break;
    case 6:
        strcpy(base_rpt, "CCCCC:GGGGG");
        break;
    case 9:
        strcpy(base_rpt, "GGT:ACC");
        break;
    case 10:
        strcpy(base_rpt, "AG:CT");
        break;
    case 11:
        strcpy(base_rpt, "AG:CT");
        break;
    case 13:
        strcpy(base_rpt, "CI:CI");
        break;
    case 14:
        strcpy(base_rpt, "ATATAT:ATATAT");
        break;
    case 15:
        strcpy(base_rpt, "GC:GC");
        break;
    case 16:
        strcpy(base_rpt, "As4T:As4T");
        break;
    case 17:
        strcpy(base_rpt, "GC:GC");
        break;
    case 18:
        strcpy(base_rpt, "A:T");
        break;
    case 19:
        strcpy(base_rpt, "A:T");
        break;
    case 20:
        strcpy(base_rpt, "A:U");
        break;
    case 21:
        strcpy(base_rpt, "I:C");
        break;
    case 22:
        strcpy(base_rpt, "A:T");
        break;
    case 23:
        strcpy(base_rpt, "G:C");
        break;
    case 24:
        strcpy(base_rpt, "I:C");
        break;
    case 25:
        strcpy(base_rpt, "A:U");
        break;
    case 26:
        strcpy(base_rpt, "X:X");
        break;
    case 27:
        strcpy(base_rpt, "X:X");
        break;
    case 28:
        strcpy(base_rpt, "s2U:s2U");
        break;
    case 29:
        strcpy(base_rpt, "s2U:s2U");
        break;
    case 37:
        strcpy(base_rpt, "A:U");
        break;
    case 38:
        strcpy(base_rpt, "A:T");
        break;
    case 39:
        strcpy(base_rpt, "AI:CT");
        break;
    case 40:
        strcpy(base_rpt, "AI:CT");
        break;
    case 41:
        strcpy(base_rpt, "AATT:AATT");
        break;
    case 43:
        strcpy(base_rpt, "A:U");
        break;
    case 48:
        strcpy(base_rpt, "AT:AT");
        break;
    case 49:
        strcpy(base_rpt, "CG:CG");
        break;
    case 50:
        strcpy(base_rpt, "GC:GC");
        break;
    case 51:
        strcpy(base_rpt, "A:T");
        break;
    case 52:
        strcpy(base_rpt, "AT:AT");
        break;
    default:
        fatal("wrong structure id# for this function\n");
    }
    fprintf(stderr, "Repeating unit: %s\n", base_rpt);

    n_rep = repeat_num();

    sprintf(fname, "%sS1.rpt", FDIR);
    num1 = number_of_atoms(fname, 1, "*");
    AtomName1 = cmatrix(1, num1, 0, 4);
    ResName1 = cmatrix(1, num1, 0, 3);
    ChainID1 = cvector(1, num1);
    ResSeq1 = lvector(1, num1);
    xyz1 = dmatrix(1, num1, 1, 3);
    Miscs1 = cmatrix(1, num1, 0, NMISC);
    read_pdb(fname, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, Miscs1, 1, "*");
    residue_idx(num1, ResSeq1, Miscs1, ChainID1, ResName1, &nr);

    sprintf(fname, "%sS2.rpt", FDIR);
    num2 = number_of_atoms(fname, 1, "*");
    AtomName2 = cmatrix(1, num2, 0, 4);
    ResName2 = cmatrix(1, num2, 0, 3);
    ChainID2 = cvector(1, num2);
    ResSeq2 = lvector(1, num2);
    xyz2 = dmatrix(1, num2, 1, 3);
    read_pdb(fname, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz2, NULL, 1, "*");
    tnum2 = n_rep * num2;

    tnum = n_rep * num1 + tnum2;

    tChainID = cvector(1, tnum);
    tAtomName = cmatrix(1, tnum, 0, 4);
    tResName = cmatrix(1, tnum, 0, 3);
    tResSeq = lvector(1, tnum);
    txyz = dmatrix(1, tnum, 1, 3);

    tAtomName2 = cmatrix(1, tnum2, 0, 4);
    tResName2 = cmatrix(1, tnum2, 0, 3);
    tResSeq2 = lvector(1, tnum2);
    txyz2 = dmatrix(1, tnum2, 1, 3);

    s2idx = lmatrix(1, n_rep, 1, 2);

    tnum = 0;
    tnum2 = 0;
    for (i = 1; i <= n_rep; i++) {
        if (lval_in_range(str_id, 38, 41) || lval_in_range(str_id, 48, 52)) {  /* 38-41/48-52 */
            dphi = (i - 1) * twist;
            dz = (i - 1) * rise;
        } else {
            dphi = -(i - 1) * twist;  /* negative */
            dz = -(i - 1) * rise;  /* negative */
        }

        /* Strand I */
        for (j = 1; j <= num1; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], ResName1[j]);
            tChainID[ik] = 'A';
            tResSeq[ik] = ResSeq1[j] + (i - 1) * nr;

            phi = deg2rad(dphi + xyz1[j][3]);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        /* Strand II */
        for (j = 1; j <= num2; j++) {
            ik = tnum2 + j;
            strcpy(tAtomName2[ik], AtomName2[j]);
            strcpy(tResName2[ik], ResName2[j]);
            tResSeq2[ik] = ResSeq2[j] + nr * n_rep + (n_rep - i) * nr;

            phi = deg2rad(dphi + xyz2[j][3]);
            txyz2[ik][1] = xyz2[j][2] * cos(phi);
            txyz2[ik][2] = xyz2[j][2] * sin(phi);
            txyz2[ik][3] = dz + xyz2[j][1];
        }

        s2idx[i][1] = tnum2;
        s2idx[i][2] = num2;

        tnum += num1;
        tnum2 += num2;
    }

    /* "reverse" strand II and combined with I */
    for (i = n_rep; i >= 1; i--) {
        for (j = s2idx[i][1] + 1; j <= s2idx[i][1] + s2idx[i][2]; j++) {
            ik = tnum + j - s2idx[i][1];
            strcpy(tAtomName[ik], tAtomName2[j]);
            strcpy(tResName[ik], tResName2[j]);
            tChainID[ik] = 'B';
            tResSeq[ik] = tResSeq2[j];
            cpxyz(txyz2[j], txyz[ik]);
        }
        tnum += s2idx[i][2];
    }

    write_pdb(tnum, tAtomName, tResName, tChainID, tResSeq, txyz, NULL, pdbfile);

    free_pdb(num1, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, Miscs1);
    free_pdb(num2, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz2, NULL);
    free_pdb(tnum, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
    free_pdb(tnum2, NULL, tAtomName2, tResName2, NULL, tResSeq2, txyz2, NULL);
    free_lmatrix(s2idx, 1, n_rep, 1, 2);
}

void triplex(long str_id, double twist, double rise, char *FDIR, char *pdbfile)
/* structures 30-34 & 42 */
{
    char triplet[BUF512], fname[BUF512];
    char *ChainID1, *ChainID2, *tChainID;
    char **AtomName1, **AtomName2, **ResName1, **ResName2;
    char **tAtomName, **tAtomName2, **tResName, **tResName2;
    char **bp_seq;  /* for calling reverse_stnd2 only */

    double dphi, dz, phi;
    double **xyz1, **xyz2, **xyz3, **txyz, **txyz2;

    long i, ik, j, n_rep, num1, num2, tnum, tnum2, offset3;
    long *ResSeq1, *ResSeq2, *tResSeq, *tResSeq2;
    long **s2idx;  /* strand II index */

    switch (str_id) {
    case 30:
        strcpy(triplet, "CIC");
        break;
    case 31:
        strcpy(triplet, "TAT");
        break;
    case 32:
        strcpy(triplet, "UAU");
        break;
    case 33:
        strcpy(triplet, "UAU");
        break;
    case 34:
        strcpy(triplet, "IAI");
        break;
    case 42:
        strcpy(triplet, "UAU");
        break;
    default:
        fatal("wrong structure id# for function <triplex>\n");
        break;
    }  /* strands 1 & 3 have the same residue */
    fprintf(stderr, "Triplex repeating unit: %s\n", triplet);

    n_rep = repeat_num();

    /* I: middle strand, R residue, 5'-->3' */
    sprintf(fname, "%s%c2.rpt", FDIR, triplet[1]);
    num1 = number_of_atoms(fname, 1, "*");
    AtomName1 = cmatrix(1, num1, 0, 4);
    ResName1 = cmatrix(1, num1, 0, 3);
    ChainID1 = cvector(1, num1);
    ResSeq1 = lvector(1, num1);
    xyz1 = dmatrix(1, num1, 1, 3);
    read_pdb(fname, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL, 1, "*");

    /* II: 3'-->5', anti-parallel to the middle strand */
    sprintf(fname, "%s%c1.rpt", FDIR, triplet[0]);
    num2 = number_of_atoms(fname, 1, "*");
    AtomName2 = cmatrix(1, num2, 0, 4);
    ResName2 = cmatrix(1, num2, 0, 3);
    ChainID2 = cvector(1, num2);
    ResSeq2 = lvector(1, num2);
    xyz2 = dmatrix(1, num2, 1, 3);
    read_pdb(fname, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz2, NULL, 1, "*");
    tnum2 = n_rep * num2;

    /* III: 5'-->3' parallel to the middle strand */
    sprintf(fname, "%s%c3.rpt", FDIR, triplet[2]);
    xyz3 = dmatrix(1, num2, 1, 3);
    read_pdb(fname, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz3, NULL, 1, "*");

    offset3 = n_rep * num1 + tnum2;
    tnum = offset3 + tnum2;

    tChainID = cvector(1, tnum);
    tAtomName = cmatrix(1, tnum, 0, 4);
    tResName = cmatrix(1, tnum, 0, 3);
    tResSeq = lvector(1, tnum);
    txyz = dmatrix(1, tnum, 1, 3);

    tAtomName2 = cmatrix(1, tnum2, 0, 4);
    tResName2 = cmatrix(1, tnum2, 0, 3);
    tResSeq2 = lvector(1, tnum2);
    txyz2 = dmatrix(1, tnum2, 1, 3);

    s2idx = lmatrix(1, n_rep, 1, 2);

    tnum = 0;
    tnum2 = 0;
    for (i = 1; i <= n_rep; i++) {
        dphi = -(i - 1) * twist;  /* negative */
        dz = -(i - 1) * rise;  /* negative */

        /* Strand I */
        for (j = 1; j <= num1; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            sprintf(tResName[ik], "  %c", triplet[1]);
            tChainID[ik] = 'A';
            tResSeq[ik] = i;

            phi = deg2rad(dphi + xyz1[j][3]);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        /* Strand II */
        for (j = 1; j <= num2; j++) {
            ik = tnum2 + j;
            strcpy(tAtomName2[ik], AtomName2[j]);

            phi = deg2rad(dphi + xyz2[j][3]);
            txyz2[ik][1] = xyz2[j][2] * cos(phi);
            txyz2[ik][2] = xyz2[j][2] * sin(phi);
            txyz2[ik][3] = dz + xyz2[j][1];
        }

        /* Strand III */
        for (j = 1; j <= num2; j++) {
            ik = offset3 + j;
            strcpy(tAtomName[ik], AtomName2[j]);
            sprintf(tResName[ik], "  %c", triplet[2]);
            tChainID[ik] = 'C';
            tResSeq[ik] = 2 * n_rep + i;

            phi = deg2rad(dphi + xyz3[j][3]);
            txyz[ik][1] = xyz3[j][2] * cos(phi);
            txyz[ik][2] = xyz3[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz3[j][1];
        }

        s2idx[i][1] = tnum2;
        s2idx[i][2] = num2;

        tnum += num1;
        tnum2 += num2;
        offset3 += num2;
    }

    /* "reverse" strand II and combined with I */
    bp_seq = cmatrix(1, n_rep, 1, 2);
    for (i = 1; i <= n_rep; i++)
        bp_seq[i][2] = triplet[0];
    reverse_stnd2(n_rep, bp_seq, s2idx, &tnum, tAtomName, tResName,
                  tChainID, tResSeq, txyz, tAtomName2, txyz2, 0);
    tnum += tnum2;

    write_pdb(tnum, tAtomName, tResName, tChainID, tResSeq, txyz, NULL, pdbfile);

    free_pdb(num1, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL);
    free_pdb(num2, NULL, AtomName2, ResName2, ChainID2, ResSeq2, xyz2, NULL);
    free_dmatrix(xyz3, 1, num2, 1, 3);
    free_pdb(tnum, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
    free_pdb(tnum2, NULL, tAtomName2, tResName2, NULL, tResSeq2, txyz2, NULL);
    free_lmatrix(s2idx, 1, n_rep, 1, 2);
    free_cmatrix(bp_seq, 1, n_rep, 1, 2);
}

void quadruplex(double twist, double rise, char *FDIR, char *pdbfile)
/* structure 35 Poly(I).Poly(I).Poly(I).Poly(I) */
{
    char fname[BUF512];
    char *ChainID1, *tChainID;
    char **AtomName1, **ResName1, **tAtomName, **tResName;
    double dphi, dz, phi;
    double **xyz1, **txyz;
    long i, ik, j, n_rep, num1, tnum;
    long offset1, offset2, offset3, offset4;
    long *ResSeq1, *tResSeq;

    /* all 4 strands are parallel */
    fprintf(stderr, "Quadruplex repeating unit: I.I.I.I\n");

    n_rep = repeat_num();

    sprintf(fname, "%sI.rpt", FDIR);
    num1 = number_of_atoms(fname, 1, "*");
    AtomName1 = cmatrix(1, num1, 0, 4);
    ResName1 = cmatrix(1, num1, 0, 3);
    ChainID1 = cvector(1, num1);
    ResSeq1 = lvector(1, num1);
    xyz1 = dmatrix(1, num1, 1, 3);
    read_pdb(fname, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL, 1, "*");

    offset1 = 0;
    offset2 = n_rep * num1;
    offset3 = 2 * n_rep * num1;
    offset4 = 3 * n_rep * num1;
    tnum = 4 * n_rep * num1;

    tChainID = cvector(1, tnum);
    tAtomName = cmatrix(1, tnum, 0, 4);
    tResName = cmatrix(1, tnum, 0, 3);
    tResSeq = lvector(1, tnum);
    txyz = dmatrix(1, tnum, 1, 3);

    for (i = 1; i <= n_rep; i++) {
        dphi = -(i - 1) * twist;  /* negative */
        dz = -(i - 1) * rise;  /* negative */

        /* Strand I */
        for (j = 1; j <= num1; j++) {
            ik = offset1 + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], "  I");
            tChainID[ik] = 'A';
            tResSeq[ik] = i;

            phi = deg2rad(dphi + xyz1[j][3]);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        /* Strand II */
        for (j = 1; j <= num1; j++) {
            ik = offset2 + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], "  I");
            tChainID[ik] = 'B';
            tResSeq[ik] = n_rep + i;

            phi = deg2rad(dphi + xyz1[j][3] + 90.0);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        /* Strand III */
        for (j = 1; j <= num1; j++) {
            ik = offset3 + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], "  I");
            tChainID[ik] = 'C';
            tResSeq[ik] = 2 * n_rep + i;

            phi = deg2rad(dphi + xyz1[j][3] + 180.0);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        /* Strand IV */
        for (j = 1; j <= num1; j++) {
            ik = offset4 + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], "  I");
            tChainID[ik] = 'D';
            tResSeq[ik] = 3 * n_rep + i;

            phi = deg2rad(dphi + xyz1[j][3] - 90.0);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        offset1 += num1;
        offset2 += num1;
        offset3 += num1;
        offset4 += num1;
    }

    write_pdb(tnum, tAtomName, tResName, tChainID, tResSeq, txyz, NULL, pdbfile);

    free_pdb(num1, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL);
    free_pdb(tnum, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
}

void shelix(double twist, double rise, char *FDIR, char *pdbfile)
/* structure 36 Poly(eC) */
{
    char fname[BUF512];
    char *ChainID1, *tChainID;
    char **AtomName1, **ResName1, **tAtomName, **tResName;
    double dphi, dz, phi;
    double **xyz1, **txyz;
    long i, ik, j, n_rep, num1, tnum;
    long *ResSeq1, *tResSeq;

    fprintf(stderr, "Single helix repeating unit: eC\n");

    n_rep = repeat_num();

    sprintf(fname, "%sC.rpt", FDIR);
    num1 = number_of_atoms(fname, 1, "*");
    AtomName1 = cmatrix(1, num1, 0, 4);
    ResName1 = cmatrix(1, num1, 0, 3);
    ChainID1 = cvector(1, num1);
    ResSeq1 = lvector(1, num1);
    xyz1 = dmatrix(1, num1, 1, 3);
    read_pdb(fname, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL, 1, "*");

    tnum = n_rep * num1;

    tChainID = cvector(1, tnum);
    tAtomName = cmatrix(1, tnum, 0, 4);
    tResName = cmatrix(1, tnum, 0, 3);
    tResSeq = lvector(1, tnum);
    txyz = dmatrix(1, tnum, 1, 3);

    tnum = 0;
    for (i = 1; i <= n_rep; i++) {
        dphi = -(i - 1) * twist;  /* negative */
        dz = -(i - 1) * rise;  /* negative */

        for (j = 1; j <= num1; j++) {
            ik = tnum + j;
            strcpy(tAtomName[ik], AtomName1[j]);
            strcpy(tResName[ik], "  C");
            tChainID[ik] = 'A';
            tResSeq[ik] = i;

            phi = deg2rad(dphi + xyz1[j][3]);
            txyz[ik][1] = xyz1[j][2] * cos(phi);
            txyz[ik][2] = xyz1[j][2] * sin(phi);
            txyz[ik][3] = dz + xyz1[j][1];
        }

        tnum += num1;
    }

    write_pdb(tnum, tAtomName, tResName, tChainID, tResSeq, txyz, NULL, pdbfile);

    free_pdb(num1, NULL, AtomName1, ResName1, ChainID1, ResSeq1, xyz1, NULL);
    free_pdb(tnum, NULL, tAtomName, tResName, tChainID, tResSeq, txyz, NULL);
}

void bp_repeat(char *Wbase, char *Cbase, char *crepeat, long *n_rep, char **bp_seq)
/* repeat unit "crepeat", adopted from function: bp_sequence */
{
    char str[BUF512];
    long i, j, k, num_bp, n = 0;

    num_bp = strlen(crepeat);  /* no checking for W/C bases */
    fprintf(stderr, "Repeating unit: %s\n", crepeat);

    fprintf(stderr, "Number of repeats (Dft: 10): ");
    if (fgets(str, sizeof str, stdin) != NULL) {
        j = sscanf(str, "%ld", &n);
        if (!j || j == EOF)
            n = 10;
    } else
        fatal("error in reading number of repeats\n");
    if (n < 0) {
        fprintf(stderr, "number of repeats < 0. reset to positive\n");
        n = -n;
    }
    if (n * num_bp > MAXBASE)
        fatal("number of bases out of limit\n");
    for (i = 1; i <= n; i++) {
        k = (i - 1) * num_bp + 1;
        for (j = 0; j < num_bp; j++)
            bp_seq[k + j][1] = crepeat[j];
    }
    num_bp *= n;

    for (i = 1; i <= num_bp; i++) {
        j = strchr(Wbase, bp_seq[i][1]) - Wbase;
        bp_seq[i][2] = *(Cbase + j);
    }

    *n_rep = num_bp;
}

long repeat_num(void)
{
    char str[BUF512];
    long i, n = 0;

    fprintf(stderr, "Number of repeats (Dft: 10): ");
    if (fgets(str, sizeof str, stdin) != NULL) {
        i = sscanf(str, "%ld", &n);
        if (!i || i == EOF)
            n = 10;
    } else
        fatal("error in reading number of repeats\n");
    if (n < 0) {
        fprintf(stderr, "number of repeats < 0. reset to positive\n");
        n = -n;
    }
    return n;
}
