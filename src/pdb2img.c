#include "x3dna.h"

int main(int argc, char *argv[])
/* PDB format to XFIG or PS schematic presentation */
{
    char alcfile[BUF512], inpfile[BUF512], outfile[BUF512];
    double scale_factor = 0.0;
    long opts[BUF512] = { 0 };  /* all initialized to zeros */
    /* 0: 0 for PS, 1 (F) for XFIG, and 2 (R) for Raster3D
     * 1: (B) frame_box
     * 2: (C) is_color
     * 3: oline: NOT used, just for compatibility with ALC2IMG
     * 4: (P) peptide
     * 5: helix_axis: NOT used
     * 6: ref_axis: NOT used
     * 7: (I) same_faces
     * 8: (N) no_header
     * 9: (U) updown
     * 10: (M) isomM */

    remove_file(BBLKALC_FILE);
    remove_file(PBLKALC_FILE);

    pdb2img_cmdline(argc, argv, inpfile, outfile, &scale_factor, opts);

    (opts[4]) ? strcpy(alcfile, PBLKALC_FILE) : strcpy(alcfile, BBLKALC_FILE);
    pdb2alc(inpfile, alcfile, opts[4]);
    process_alc(alcfile, outfile, scale_factor, opts);

    return 0;
}

void pdb2img_usage(void)
{
    help3dna("pdb2img");
    contact_msg(0);
}

void pdb2img_cmdline(int argc, char *argv[], char *inpfile, char *outfile,
                     double *scale_factor, long *opts)
{
    long i, j;

    if (argc < 3)
        pdb2img_usage();
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'F')  /* XFIG */
                opts[0] = 1;
            else if (argv[i][j] == 'R')  /* Raster3D */
                opts[0] = 2;
            else if (argv[i][j] == 'B')  /* frame box */
                opts[1] = 1;
            else if (argv[i][j] == 'C')  /* color image */
                opts[2] = 1;
            else if (argv[i][j] == 'P')  /* peptide */
                opts[4] = 1;
            else if (argv[i][j] == 'I')  /* iso-faces for block models */
                opts[7] = 1;
            else if (argv[i][j] == 'N')  /* no-header section for Raster3D output */
                opts[8] = 1;
            else if (argv[i][j] == 'U')  /* up-down faces */
                opts[9] = 1;
            else if (argv[i][j] == 'M')  /* 5 + minor face */
                opts[10] = 1;
            else if (argv[i][j] == 'S') {
                if (sscanf(argv[i], "-S=%lf", scale_factor) != 1) {
                    fprintf(stderr, "wrong format for setting scale factor\n");
                    pdb2img_usage();
                } else
                    break;  /* -s=factor not combined with others */
            } else
                pdb2img_usage();
    }
    if (argc == i + 2) {
        strcpy(inpfile, argv[i]);
        strcpy(outfile, argv[i + 1]);
    } else
        pdb2img_usage();
    if (opts[10])  /* is 4 + minor/major faces, set -i option */
        opts[7] = 1;
}

void pdb2alc(char *inpfile, char *alcfile, long peptide)
{
    char BDIR[BUF512];
    char *ChainID, *bseq, **AtomName, **ResName, **Miscs;
    double **orien, **org, **xyz;
    long np, num, num_residue;
    long *cidx, *mchain, *resid, *ResSeq, *RY, **seidx;

    /* read in the PDB file */
    num = number_of_atoms(inpfile, 1, "*");
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");

    /* get the numbering information of each residue */
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

    /* get the reference frame for each base */
    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 4);  /* column 4 for RMS of peptide bond */

    if (peptide) {  /* for peptide bonds */
        get_BDIR(BDIR, "trans_pep.pdb", 1);

        cidx = lvector(1, num_residue);
        resid = lvector(1, num_residue);
        np = 4 * num_residue;
        mchain = lvector(1, np);
        peptide_info(num_residue, seidx, AtomName, ResName, ChainID,
                     ResSeq, Miscs, xyz, resid, cidx, mchain);
        peptide_frame(num_residue, BDIR, resid, mchain, xyz, orien, org);
        peptide_blks(num_residue, BDIR, cidx, orien, org, alcfile);

        free_lvector(cidx, 1, num_residue);
        free_lvector(resid, 1, num_residue);
        free_lvector(mchain, 1, np);
    } else {  /* for base residues */
        /* get base sequence, RY identification & serial number */
        bseq = cvector(1, num_residue);
        RY = lvector(1, num_residue);
        get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz,
                bseq, RY);

        get_BDIR(BDIR, "Atomic_A.pdb", 1);
        base_frame(num_residue, bseq, seidx, RY, AtomName, ResName,
                   ChainID, ResSeq, Miscs, xyz, BDIR, orien, org);

        /* get the block geometry in ALCHEMY format */
        get_BDIR(BDIR, "Block_R.alc", 1);
        base_blks(num_residue, RY, orien, org, bseq, BDIR, alcfile);

        free_cvector(bseq, 1, num_residue);
        free_lvector(RY, 1, num_residue);
    }

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 4);
}

void base_blks(long num_residue, long *RY, double **orien, double **org, char *bseq,
               char *BDIR, char *alcfile)
{
    char fname[BUF512];
    char **AtomName, **tAtomName;
    double temp[4];
    double **txyz, **xyz;
    long bidx, i, j, k, m, nbond, num, tnbond, tnum;
    long ia = 0, ib = 0, nb = 0;
    long *ibase, *tibase;
    long **linkage, **tlinkage;

    for (i = 1; i <= num_residue; i++)
        if (RY[i] >= 0)
            nb++;  /* number of bases */
    if (!nb)
        fatal("No base residues in this structure\n");

    strcpy(fname, BDIR);
    strcat(fname, "Block_R.alc");
    get_alc_nums(fname, &num, &nbond);
    tnum = num * nb;
    tnbond = nbond * nb;

    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 1, 2);
    tAtomName = cmatrix(1, tnum, 0, 2);
    txyz = dmatrix(1, tnum, 1, 3);
    tibase = lvector(1, tnum);
    tlinkage = lmatrix(1, tnbond, 1, 2);

    for (i = 1; i <= num_residue; i++) {
        if (RY[i] == 1)
            sprintf(fname, "%sBlock_R.alc", BDIR);
        else if (RY[i] == 0)
            sprintf(fname, "%sBlock_Y.alc", BDIR);
        else
            continue;
        read_alc(fname, &num, &nbond, AtomName, xyz, ibase, linkage);

        base_idx(i, &bseq[i], &bidx, 1);
        for (j = 1; j <= num; j++) {
            strcpy(tAtomName[ia + j], AtomName[j]);
            tibase[ia + j] = bidx;
            for (k = 1; k <= 3; k++) {
                for (m = 1; m <= 3; m++)
                    temp[m] = orien[i][k + (m - 1) * 3];
                txyz[ia + j][k] = dot(xyz[j], temp) + org[i][k];
            }
        }

        for (j = 1; j <= nbond; j++)
            for (k = 1; k <= 2; k++)
                tlinkage[ib + j][k] = ia + linkage[j][k];

        ia += num;
        ib += nbond;
    }
    write_alc(tnum, tnbond, tAtomName, txyz, tibase, tlinkage, alcfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_alc(tnum, tnbond, tAtomName, txyz, tibase, 1, tlinkage);
}

    /* ----- for block representation of peptide bond in protein ----- */
void peptide_info(long num_residue, long **seidx, char **AtomName, char **ResName,
                  char *ChainID, long *ResSeq, char **Miscs, double **xyz, long *resid,
                  long *cidx, long *mchain)
/* get residue index, color index and main chain atoms for the peptide bond */
{
    static char *aa3_code[20] = { "ALA", "VAL", "PHE", "PRO", "MET", "ILE", "LEU",
        /* hydrophobic */
        "ASP", "GLU", "LYS", "ARG",  /* charged */
        "SER", "THR", "TYR", "HIS", "CYS", "ASN", "GLN", "TRP",
        /* polar */
        "GLY"  /* smallest */
    };  /* cf. function <residue_wtype> in GET_PART */
    char *col_code = "GGGGGGGAAAATTTTTTTTCX";
    char *cmn_base = CX_LIST, idmsg[BUF512];
    long i, ib, ie, ioffset, j, naa = 20;

    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];
        resid[i] = residue_ident(AtomName, xyz, Miscs, ib, ie);
        if (resid[i] != -1)  /* non-amino-acid */
            continue;
        sprintf(idmsg, ": residue name %s, chain %c, number %4ld%c",
                ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2]);
        ioffset = (i - 1) * 4;
        mchain[ioffset + 1] = find_1st_atom(" N  ", AtomName, ib, ie, idmsg);
        mchain[ioffset + 2] = find_1st_atom(" CA ", AtomName, ib, ie, idmsg);
        mchain[ioffset + 3] = find_1st_atom(" C  ", AtomName, ib, ie, idmsg);
        mchain[ioffset + 4] = find_1st_atom(" O  ", AtomName, ib, ie, idmsg);
        for (j = 0; j < naa; j++)
            if (!strcmp(ResName[ib], aa3_code[j]))
                break;
        cidx[i] = strchr(cmn_base, col_code[j]) - cmn_base;
    }
}

void peptide_frame(long num_residue, char *BDIR, long *resid, long *mchain, double **xyz,
                   double **orien, double **org)
/* get the local reference frame for each peptide unit defined by 6 atoms:
   CA(i), C(i), O(i), N(i+1), H(i+1) & CA(i+1) */
{
    char fname[BUF512], *pfile = "trans_pep.pdb";
    char *sChainID, **sAtomName, **sResName;
    double rms_fit, orgi[4];
    double **ePxyz, **fitted_xyz, **sPxyz, **sxyz, **R;
    long i, ioffset, j, nmatch, pnum = 6;
    long idx[7], *sResSeq;

    strcpy(fname, BDIR);
    strcat(fname, pfile);
    fprintf(stderr, " ...... reading file: %s ...... \n", pfile);

    sAtomName = cmatrix(1, pnum, 0, 4);
    sResName = cmatrix(1, pnum, 0, 3);
    sChainID = cvector(1, pnum);
    sResSeq = lvector(1, pnum);
    sxyz = dmatrix(1, pnum, 1, 3);
    read_pdb(fname, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL, 1, "*");

    ePxyz = dmatrix(1, pnum, 1, 3);
    sPxyz = dmatrix(1, pnum, 1, 3);
    fitted_xyz = dmatrix(1, pnum, 1, 3);
    R = dmatrix(1, 3, 1, 3);
    for (i = 1; i < num_residue; i++) {
        org[i][4] = -1.0;  /* as a default indicator */
        if (resid[i] != -1 || resid[i + 1] != -1)  /* non-amino-acid */
            continue;
        ioffset = (i - 1) * 4;
        idx[1] = mchain[ioffset + 2];  /* CA (i) */
        idx[2] = mchain[ioffset + 3];  /* C (i) */
        idx[3] = mchain[ioffset + 4];  /* O (i) */
        ioffset += 4;  /* for residue i + 1 */
        idx[4] = mchain[ioffset + 1];  /* N (i + 1) */
        idx[5] = 0;  /* H (i + 1): not counted */
        idx[6] = mchain[ioffset + 2];  /* CA (i + 1) */

        nmatch = 0;
        for (j = 1; j <= pnum; j++) {
            if (idx[j]) {
                ++nmatch;
                cpxyz(xyz[idx[j]], ePxyz[nmatch]);
                cpxyz(sxyz[j], sPxyz[nmatch]);
            }
        }

        rms_fit = ls_fitting(sPxyz, ePxyz, nmatch, fitted_xyz, R, orgi);

        org[i][4] = rms_fit;  /* to be checked if too large */
        cpxyz(orgi, org[i]);
        for (j = 1; j <= 3; j++) {
            orien[i][j] = R[j][1];  /* x-axis */
            orien[i][j + 3] = R[j][2];  /* y-axis */
            orien[i][j + 6] = R[j][3];  /* z=axis */
        }
    }

    free_pdb(pnum, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz, NULL);
    free_dmatrix(ePxyz, 1, pnum, 1, 3);
    free_dmatrix(sPxyz, 1, pnum, 1, 3);
    free_dmatrix(fitted_xyz, 1, pnum, 1, 3);
    free_dmatrix(R, 1, 3, 1, 3);
}

void peptide_blks(long num_residue, char *BDIR, long *cidx, double **orien, double **org,
                  char *alcfile)
{
    char fname[BUF512], *pfile = "trans_pep.alc";
    char **AtomName, **tAtomName;
    double tprms = 0.5, temp[4], **txyz, **xyz;
    long i, j, k, m, nbond, num, tnbond, tnum;
    long ia = 0, ib = 0, num_pep = 0;
    long *ibase, *tibase, **linkage, **tlinkage;

    strcpy(fname, BDIR);
    strcat(fname, pfile);
    fprintf(stderr, " ...... reading file: %s ...... \n", pfile);

    get_alc_nums(fname, &num, &nbond);
    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);  /* of no use here */
    linkage = lmatrix(1, nbond, 1, 2);
    read_alc(fname, &num, &nbond, AtomName, xyz, ibase, linkage);

    for (i = 1; i < num_residue; i++)
        if (dval_in_range(org[i][4], 0.0, tprms))
            num_pep++;
    if (!num_pep)
        fatal("No peptide bonds in this structure\n");
    tnum = num * num_pep;
    tnbond = nbond * num_pep;

    tAtomName = cmatrix(1, tnum, 0, 2);
    txyz = dmatrix(1, tnum, 1, 3);
    tibase = lvector(1, tnum);
    tlinkage = lmatrix(1, tnbond, 1, 2);

    for (i = 1; i < num_residue; i++) {
        if (!dval_in_range(org[i][4], 0.0, tprms))
            continue;  /* not a trans-form peptide bond */
        for (j = 1; j <= num; j++) {
            strcpy(tAtomName[ia + j], AtomName[j]);
            tibase[ia + j] = cidx[i];
            for (k = 1; k <= 3; k++) {
                for (m = 1; m <= 3; m++)
                    temp[m] = orien[i][k + (m - 1) * 3];
                txyz[ia + j][k] = dot(xyz[j], temp) + org[i][k];
            }
        }
        for (j = 1; j <= nbond; j++)
            for (k = 1; k <= 2; k++)
                tlinkage[ib + j][k] = ia + linkage[j][k];
        ia += num;
        ib += nbond;
    }
    write_alc(tnum, tnbond, tAtomName, txyz, tibase, tlinkage, alcfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_alc(tnum, tnbond, tAtomName, txyz, tibase, 1, tlinkage);
}
