#include "x3dna.h"

int main(int argc, char *argv[])
/* A 3DNA generated structure is always oriented w.r.t. the first
   base-pair (or base). A data file containing the reference frame
   (default name "ref_frames.dat") for each base-pair (or base) is
   also available, which is used here to set the new reference frame
   for the structure. The new frame could that of any base-pair (or
   base) or the middle frame or middle helical frame define by any two
   base-pairs (or bases) */
{
    char parfile[BUF512], inpfile[BUF512], outfile[BUF512], **bp_seq;
    double morg[4], **org, **orien, **mst;
    long side_view = 0, helix_axis = 0, ref_axis = 0;
    long bp_str = 0, is_helical = 0, n1, n2 = 0, frame = 0;
    long ds, num, num_bp;

    frame_mol_cmdline(argc, argv, parfile, inpfile, outfile, &side_view, &helix_axis,
                      &ref_axis, &is_helical, &n1, &n2, &bp_str, &frame);

    num_bp = get_num_bp(parfile);

    org = dmatrix(1, num_bp, 1, 3);
    orien = dmatrix(1, num_bp, 1, 9);

    bp_seq = cmatrix(1, num_bp, 1, 2);
    read_frames(parfile, &num_bp, bp_seq, org, orien, &ds);

    if (bp_str && ds == 2) {
        ref_duplex_block(num_bp, bp_seq, org, orien, outfile);
        return 0;
    } else if (frame) {
        pdbalc_frame(ref_axis, num_bp, org, orien, inpfile, outfile);
        return 0;
    }

    if (!lval_in_range(n1, 1, num_bp))
        fatal("base-pair number out of range\n");
    if (n2 && (!lval_in_range(n2, 1, num_bp) || n2 <= n1))
        fatal("base-pair number out of range\n");

    /* define the new reference frame  */
    mst = dmatrix(1, 3, 1, 3);
    get_mstorg(n1, n2, is_helical, org, orien, morg, mst);

    num = number_of_atoms(inpfile, 1, "*");
    (num) ? frame_pdbfile(num, side_view, inpfile, outfile, morg, mst) :
        frame_alcfile(num_bp, side_view, helix_axis, ref_axis, inpfile,
                      outfile, org, orien, morg, mst);

    free_cmatrix(bp_seq, 1, num_bp, 1, 2);
    free_dmatrix(org, 1, num_bp, 1, 3);
    free_dmatrix(orien, 1, num_bp, 1, 9);
    free_dmatrix(mst, 1, 3, 1, 3);

    return 0;
}

void frame_mol_usage(void)
{
    help3dna("frame_mol");
    contact_msg(0);
}

void frame_mol_cmdline(int argc, char *argv[], char *parfile, char *inpfile,
                       char *outfile, long *side_view, long *helix_axis, long *ref_axis,
                       long *is_helical, long *n1, long *n2, long *bp_str, long *frame)
{
    long i, j, n = 0;

    if (argc < 4)
        frame_mol_usage();
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'M')
                *side_view = 1;
            else if (argv[i][j] == 'A')
                *helix_axis = 1;
            else if (argv[i][j] == 'G')
                *ref_axis = 1;
            else if (argv[i][j] == 'S')
                *bp_str = 1;
            else if (argv[i][j] == 'F')
                *frame = 1;
            else if (argv[i][j] == 'X')
                *is_helical = 1;
            else if (isdigit((int) argv[i][j])) {
                n = sscanf(argv[i], "-%ld,%ld", n1, n2);
                if (n == 1) {
                    *n2 = 0;  /* just to make sure */
                    if (*is_helical) {
                        fprintf(stderr, "two numbers -n1,n2 required to"
                                " define a the helical frame\n");
                        frame_mol_usage();
                    }
                } else if (n == 2 || *frame == 1);  /* okay */
                else
                    frame_mol_usage();
                break;
            } else  /* wrong options or -h */
                frame_mol_usage();
    }
    if (argc == i + 3 && (n || *frame)) {
        strcpy(parfile, argv[i]);
        strcpy(inpfile, argv[i + 1]);
        strcpy(outfile, argv[i + 2]);
    } else if (argc == i + 2 && *bp_str) {
        strcpy(parfile, argv[i]);
        strcpy(outfile, argv[i + 1]);
    } else
        frame_mol_usage();
}

void read_frames(char *parfile, long *num_bp, char **bp_seq, double **org,
                 double **orien, long *ds)
/* read in reference frames: origin and orientation */
{
    char *format = "%lf %lf %lf", stmp[BUF512];
    char str[BUF512];
    long i, num1 = 0, num2 = 0;
    FILE *fp;

    fp = open_file(parfile, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", num_bp) != 1)
        fatal("error reading ref-frame file: number of base-pairs\n");

    for (i = 1; i <= *num_bp; i++) {
        if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%*s %*d %s", stmp) != 1)
            fatal("error reading title line\n");
        else {
            bp_seq[i][1] = stmp[0];
            if (strlen(stmp) == 3) {
                num2++;
                bp_seq[i][2] = stmp[2];
            } else
                num1++;
        }
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &org[i][1], &org[i][2], &org[i][3]) != 3)
            fatal("error reading origin line\n");
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &orien[i][1], &orien[i][2], &orien[i][3]) != 3)
            fatal("error reading x direction cosines\n");
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &orien[i][4], &orien[i][5], &orien[i][6]) != 3)
            fatal("error reading y direction cosines\n");
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &orien[i][7], &orien[i][8], &orien[i][9]) != 3)
            fatal("error reading z direction cosines\n");
    }
    if (num1 && num2)
        fatal("duplex and single helix mixed\n");
    if ((num1 && num1 != *num_bp) || (num2 && num2 != *num_bp))
        fatal("reference frame data file not consistent\n");
    *ds = (num2) ? 2 : 1;

    close_file(fp);
}

void attach_helix(long nbpm1, double **orien, double **org, long *num, char **AtomName,
                  long *ibase, double **xyz, long *nbond, long **linkage)
/* attach helix axes to ALCHEMY file */
{
    double half_rise;
    double morg[4], pars[7];
    double **mst, **rot1, **rot2;
    long i, j, k;

    rot1 = dmatrix(1, 3, 1, 3);
    rot2 = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= nbpm1; i++) {
        orien2mst(orien[i], 0, rot1);
        orien2mst(orien[i + 1], 0, rot2);

        helical_par(rot1, org[i], rot2, org[i + 1], pars, mst, morg);
        k = *num + 2 * i - 1;
        half_rise = 0.5 * pars[3];
        for (j = 1; j <= 3; j++) {
            xyz[k][j] = morg[j] - half_rise * mst[j][3];
            xyz[k + 1][j] = morg[j] + half_rise * mst[j][3];
        }
        for (j = k; j <= k + 1; j++) {
            strcpy(AtomName[j], "OH");
            ibase[j] = NON_WC_IDX;  /* non-common base atoms */
        }
        linkage[*nbond + i][1] = k;
        linkage[*nbond + i][2] = k + 1;
    }
    *num += 2 * nbpm1;
    *nbond += nbpm1;

    free_dmatrix(rot1, 1, 3, 1, 3);
    free_dmatrix(rot2, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
}

void ref_duplex_block(long num_bp, char **bp_seq, double **org, double **orien,
                      char *outfile)
/* generate base-pair block representation from reference frame */
{
    char bpi[3], BDIR[BUF512], **AtomName, **tAtomName;
    double **rmtx, **xyz, **txyz;
    long i, ik, ia = 0, ib = 0, j, k, bidx, nbond, num, tnbond, tnum;
    long *ibase, *tibase, **linkage, **tlinkage;

    get_BDIR(BDIR, "Block_BP.alc", 1);
    strcat(BDIR, "Block_BP.alc");
    get_alc_nums(BDIR, &num, &nbond);

    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 1, 2);

    tnum = (num + 1) * num_bp;  /* including bp origin */
    tnbond = (nbond + 1) * num_bp - 1;  /* bp origin linkage */
    tAtomName = cmatrix(1, tnum, 0, 2);
    txyz = dmatrix(1, tnum, 1, 3);
    tibase = lvector(1, tnum);
    tlinkage = lmatrix(1, tnbond, 1, 2);

    rmtx = dmatrix(1, 3, 1, 3);

    read_alc(BDIR, &num, &nbond, AtomName, xyz, ibase, linkage);

    for (i = 1; i <= num_bp; i++) {
        orien2mst(orien[i], 0, rmtx);  /* get rotation matrix */

        sprintf(bpi, "%c%c", toupper(bp_seq[i][1]), toupper(bp_seq[i][2]));
        bidx = basepair_idx(bpi);

        for (j = 1; j <= num; j++) {
            ik = ia + j;
            strcpy(tAtomName[ik], AtomName[j]);
            tibase[ik] = bidx;
            for (k = 1; k <= 3; k++)
                txyz[ik][k] = dot(xyz[j], rmtx[k]) + org[i][k];
        }

        for (j = 1; j <= nbond; j++) {
            ik = ib + j;
            for (k = 1; k <= 2; k++)
                tlinkage[ik][k] = ia + linkage[j][k];
        }

        ia += num;
        ib += nbond;

    }

    cnct_org(num_bp, ia, ib, tAtomName, txyz, tibase, tlinkage, org);

    write_alc(tnum, tnbond, tAtomName, txyz, tibase, tlinkage, outfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_alc(tnum, tnbond, tAtomName, txyz, tibase, 1, tlinkage);
    free_dmatrix(rmtx, 1, 3, 1, 3);
}

void pdbalc_frame(long ref_axis, long num_bp, double **org, double **orien,
                  char *inpfile, char *outfile)
/* Uncommented: generate reference frame axes for illustration purpose */
{
    static char *ANPDB[4] = { " OO ", " OX ", " OY ", " OZ " };
    static char ANALC[3] = "OH";
    static char *ref_symbol[4] = { "OO", "OX", "OY", "OZ" };
    double d[4][3];
    static double ref_xyz[4][3] = {
        {0.0, 0.0, 0.0},  /* origin */
        {1.0, 0.0, 0.0},  /* x-axis */
        {0.0, 1.0, 0.0},  /* y-axis */
        {0.0, 0.0, 1.0}  /* z-axis */
    };
    long i, j, k, vnum0, vnum, nbond;
    long ipp = 0, ipa = 0, ipn = 0, ir = 0;
    FILE *fpdb, *falc;

    fpdb = open_file(inpfile, "w");
    falc = open_file(outfile, "w");

    vnum = vnum0 = 6 * num_bp;
    nbond = 3 * num_bp;
    if (ref_axis) {  /* with global reference frame */
        vnum += 4;
        nbond += 3;
    }

    /* Check later on with write_pdb/alc */
    fprintf(falc, "%5ld ATOMS, %5ld BONDS\n", vnum, nbond);

    for (i = 1; i <= num_bp; i++) {
        for (j = 1; j <= 3; j++) {
            d[0][j - 1] = org[i][j];
            for (k = 1; k <= 3; k++)
                d[j][k - 1] = BOND_UPPER_LIMIT * orien[i][(j - 1) * 3 + k] + org[i][k];
        }
        ir++;
        for (j = 0; j <= 3; j++)
            fprintf(fpdb, "ATOM  %5ld %4s %3s %c%4ld    %8.3f%8.3f%8.3f\n", ++ipp,
                    ANPDB[j], "AXE", 'Z', ir, d[j][0], d[j][1], d[j][2]);
        for (j = 1; j <= 3; j++) {
            fprintf(falc, "%5ld %2s   %9.4f%9.4f%9.4f\n", ++ipa, ANALC, d[0][0], d[0][1], d[0][2]);  /* origin */
            fprintf(falc, "%5ld %2s   %9.4f%9.4f%9.4f\n", ++ipa, ANALC, d[j][0], d[j][1], d[j][2]);  /* x, y, z tips */
        }
    }

    if (ref_axis) {  /* add global reference frame */
        ir++;
        for (j = 0; j <= 3; j++) {
            fprintf(fpdb, "ATOM  %5ld %4s %3s %c%4ld    ", ++ipp,
                    ANPDB[j], "AXE", 'Z', ir);
            fprintf(falc, "%5ld %2s   ", ++ipa, ref_symbol[j]);
            for (k = 0; k <= 2; k++) {
                fprintf(fpdb, "%8.3f", AXIS_LENGTH * ref_xyz[j][k]);
                fprintf(falc, "%9.4f", AXIS_LENGTH * ref_xyz[j][k]);
            }
            fprintf(fpdb, "\n");
            fprintf(falc, "\n");
        }
    }

    for (i = 1; i <= ir; i++) {
        k = (i - 1) * 4;
        fprintf(fpdb, "CONECT%5ld%5ld%5ld%5ld\n", k + 1, k + 2, k + 3, k + 4);
        fprintf(fpdb, "CONECT%5ld%5ld\n", k + 2, k + 1);
        fprintf(fpdb, "CONECT%5ld%5ld\n", k + 3, k + 1);
        fprintf(fpdb, "CONECT%5ld%5ld\n", k + 4, k + 1);
    }

    ipn = 0;
    for (i = 1; i <= vnum0; i += 2)
        fprintf(falc, "%5ld %5ld %5ld\n", ++ipn, i, i + 1);
    if (ref_axis)
        for (j = 1; j <= 3; j++)
            fprintf(falc, "%5ld %5ld %5ld\n", ++ipn, vnum0 + 1, vnum0 + j + 1);

    close_file(fpdb);
    close_file(falc);
}

void get_mstorg(long n1, long n2, long is_helical, double **org, double **orien,
                double *morg, double **mst)
/* Get mst and org */
{
    double pars[7];
    double **rot1, **rot2;
    long i, bz = 1, bz_junction = 0, z_step = 0;

    /* define the new reference frame  */
    rot1 = dmatrix(1, 3, 1, 3);
    rot2 = dmatrix(1, 3, 1, 3);
    if (!n2) {  /* single base-pair */
        cpxyz(org[n1], morg);
        orien2mst(orien[n1], 0, mst);
    } else {  /* two base-pairs */
        orien2mst(orien[n1], 0, rot1);
        orien2mst(orien[n2], 0, rot2);

        bz_check(rot1, org[n1], rot2, org[n2], bz, &bz_junction, &z_step);
        if (is_helical)
            helical_par(rot1, org[n1], rot2, org[n2], pars, mst, morg);
        else
            bpstep_par(rot1, org[n1], rot2, org[n2], pars, mst, morg);
        fprintf(stderr, "the six parameters are:\n");
        for (i = 1; i <= 6; i++)
            fprintf(stderr, "%8.2f", pars[i]);
        fprintf(stderr, "\n");
    }

    free_dmatrix(rot1, 1, 3, 1, 3);
    free_dmatrix(rot2, 1, 3, 1, 3);
}

void frame_pdbfile(long num, long side_view, char *inpfile, char *outfile, double *morg,
                   double **mst)
{
    char *ChainID, **AtomName, **ResName, **Miscs;
    double **xyz;
    long *ResSeq;

    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);

    read_pdb(inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");
    change_xyz(side_view, morg, mst, num, xyz);
    write_pdb(num, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, outfile);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
}

void frame_alcfile(long num_bp, long side_view, long helix_axis, long ref_axis,
                   char *inpfile, char *outfile, double **org, double **orien,
                   double *morg, double **mst)
{
    char **AtomName;
    double **xyz;
    long nbond, nbpm1, num, tnbond, tnum;
    long *ibase, **linkage;

    nbpm1 = num_bp - 1;

    get_alc_nums(inpfile, &num, &nbond);
    tnum = num + 2 * nbpm1 + 4;
    tnbond = nbond + nbpm1 + 3;

    AtomName = cmatrix(1, tnum, 0, 2);
    xyz = dmatrix(1, tnum, 1, 3);
    ibase = lvector(1, tnum);
    linkage = lmatrix(1, tnbond, 1, 2);

    read_alc(inpfile, &num, &nbond, AtomName, xyz, ibase, linkage);

    if (helix_axis)
        attach_helix(nbpm1, orien, org, &num, AtomName, ibase, xyz, &nbond, linkage);

    change_xyz(side_view, morg, mst, num, xyz);

    if (ref_axis)
        add_3axes(&num, AtomName, ibase, xyz, &nbond, linkage, side_view);

    write_alc(num, nbond, AtomName, xyz, ibase, linkage, outfile);

    free_alc(tnum, tnbond, AtomName, xyz, ibase, 1, linkage);
}
