#include "x3dna.h"

int main(int argc, char *argv[])
/* set a standard base geometry using three parameters: does NOT work for P/p
 *
 * Revised on 03-27-2006: using Atomic_?.pdb as the reference to reset
 * an experimental residue.
*/
{
    char bname = '\0';
    char inpfile[BUF512], outfile[BUF512];
    char *ChainID;
    char **AtomName, **ResName, **Miscs;
    double c1_n_org, dist_n2org, n_org_x;
    double n_c[4], n_c1[4], org[4], vec_n_org[4], x[4], y[4], z[4];
    double **nxyz, **xyz;
    long i, idxC1, idxN = 0, idxC = 0, ry, num, fit_std = 0;
    long *ResSeq;

    std_base_cmdline(argc, argv, inpfile, outfile, &bname, &c1_n_org, &dist_n2org,
                     &n_org_x, &fit_std);

    /* read in the PDB file */
    num = number_of_atoms(inpfile, 1, "*");
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");

    i = num_strmatch(" C1'", AtomName, 1, num);
    if (i != 1)
        fatal("one residue with C1' atom required\n");

    if (bname && !fit_std)  /* should we keep the original? */
        for (i = 1; i <= num; i++)
            sprintf(ResName[i], "  %c", bname);
    for (i = 1; i <= num; i++) {
        ChainID[i] = 'A';
        ResSeq[i] = 1;
    }

    idxC1 = find_1st_atom(" C1'", AtomName, 1, num, inpfile);

    ry = residue_ident(AtomName, xyz, Miscs, 1, num);
    if (ry < 0)
        fatal("this residue is NOT a legal base\n");
    else if (ry == 1) {  /* R base */
        idxN = find_1st_atom(" N9 ", AtomName, 1, num, inpfile);
        idxC = find_1st_atom(" C8 ", AtomName, 1, num, inpfile);
    } else {  /* Y base */
        idxN = find_1st_atom(" N1 ", AtomName, 1, num, inpfile);
        idxC = find_1st_atom(" C6 ", AtomName, 1, num, inpfile);
    }

    nxyz = dmatrix(1, num, 1, 3);
    if (fit_std) {
        fit2standard_base(num, AtomName, xyz, bname, ry, nxyz);
        goto CLN_UP;
    }

    /* project all atoms onto the least-squares base plane */
    prj2plane(num, (ry == 1) ? 9 : 6, AtomName, xyz, 0, nxyz);

    ddxyz(nxyz[idxN], nxyz[idxC], n_c);  /* RN9/YN1--->RC8/YC6 */
    ddxyz(nxyz[idxN], nxyz[idxC1], n_c1);  /* RN9/YN1--->C1' */
    cross(n_c, n_c1, z);
    vec_norm(z);

    /* get the origin of the standard base */
    get_vector(n_c1, z, c1_n_org, vec_n_org);
    for (i = 1; i <= 3; i++)
        org[i] = nxyz[idxN][i] + dist_n2org * vec_n_org[i];

    /* get x- and y-axes */
    get_vector(vec_n_org, z, n_org_x, x);
    cross(z, x, y);

    /* reset the base to its reference frame */
    move_position(xyz, num, 3, org);
    for (i = 1; i <= num; i++) {
        nxyz[i][1] = dot(xyz[i], x);
        nxyz[i][2] = dot(xyz[i], y);
        nxyz[i][3] = dot(xyz[i], z);
    }

  CLN_UP:
    /* write the standard base geometry */
    write_pdb(num, AtomName, ResName, ChainID, ResSeq, nxyz, Miscs, outfile);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_dmatrix(nxyz, 1, num, 1, 3);

    return 0;
}

void fit2standard_base(long num, char **AtomName, double **xyz, char bname, long ry,
                       double **nxyz)
{
    static char *RingAtom[] = { RA_LIST };
    char BDIR[BUF512], spdb[BUF512];
    char *sChainID, **sAtomName, **sResName, **sMiscs;

    double org[4];
    double **eRing_xyz, **fitted_xyz, **sRing_xyz, **sxyz, **R;

    long i, RingAtom_num, exp_katom, std_katom, nmatch, snum;
    long *sResSeq;

    sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
    sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
    sChainID = cvector(1, NUM_RESIDUE_ATOMS);
    sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
    sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
    sMiscs = cmatrix(1, NUM_RESIDUE_ATOMS, 0, NMISC);

    sRing_xyz = dmatrix(1, 9, 1, 3);
    eRing_xyz = dmatrix(1, 9, 1, 3);
    fitted_xyz = dmatrix(1, 9, 1, 3);
    R = dmatrix(1, 3, 1, 3);

    RingAtom_num = (ry == 1) ? 9 : 6;

    get_BDIR(BDIR, "Atomic_A.pdb", 1);
    sprintf(spdb, "%sAtomic_%c.pdb", BDIR, bname);
    snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
                    sMiscs, 1, "*");

    nmatch = 0;
    for (i = 0; i < RingAtom_num; i++) {
        exp_katom = find_1st_atom(RingAtom[i], AtomName, 1, num, "");
        std_katom = find_1st_atom(RingAtom[i], sAtomName, 1, snum, "");
        if (exp_katom && std_katom) {
            ++nmatch;
            cpxyz(xyz[exp_katom], eRing_xyz[nmatch]);
            cpxyz(sxyz[std_katom], sRing_xyz[nmatch]);
        }
    }

    ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, org);

    copy_matrix(xyz, num, 3, nxyz);
    change_xyz(0, org, R, num, nxyz);

    free_dmatrix(eRing_xyz, 1, 9, 1, 3);
    free_dmatrix(sRing_xyz, 1, 9, 1, 3);
    free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
             sMiscs);
    free_dmatrix(fitted_xyz, 1, 9, 1, 3);
    free_dmatrix(R, 1, 3, 1, 3);
}

void std_base_usage(void)
{
    help3dna("std_base");
    contact_msg(0);
}

void std_base_cmdline(int argc, char *argv[], char *inpfile, char *outfile, char *bname,
                      double *c1_n_org, double *dist_n2org, double *n_org_x,
                      long *fit_std)
{
    static char options[5] = "ATUGC";
    /* three parameters used by 3DNA */
    double pars_used[6][3] = {
        {141.51, 4.68, 73.99},  /* A */
        {141.43, 4.68, 74.07},  /* T */
        {141.43, 4.68, 74.07},  /* U */
        {141.31, 4.73, 74.19},  /* G */
        {141.60, 4.72, 74.20},  /* C */
        {141.46, 4.70, 74.10}  /* overall mean */
    };
    long i, j, nbase = 5, idx_base = 5;

    if (argc < 3)
        std_base_usage();

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        upperstr(argv[i]);
        if (!strcmp(argv[i], "-H"))
            std_base_usage();
        for (j = 0; j < nbase; j++)
            if (argv[i][1] == options[j]) {
                *bname = options[j];
                idx_base = j;
                break;
            }
        if (j >= nbase) {
            if (!strncmp(argv[i], "-FIT", 4)) {
                *fit_std = 1;
                continue;
            } else
                std_base_usage();
        }
    }

    if (*fit_std && *bname == '\0')
        fatal("you must set an one-letter option (e.g., -A) with -fit\n");

    j = argc - i;  /* command-line options left over */
    if (j == 2) {
        strcpy(inpfile, argv[i]);
        strcpy(outfile, argv[i + 1]);
    } else if (j == 5) {
        if (sscanf(argv[i], "%lf", c1_n_org) != 1)
            fatal("error reading C1'--RN9/YN1--Origin angle\n");
        if (sscanf(argv[i + 1], "%lf", dist_n2org) != 1)
            fatal("error reading RN9/YN1--Origin distance\n");
        if (sscanf(argv[i + 2], "%lf", n_org_x) != 1)
            fatal("error reading RN9/YN1--Origin--X-Axis angle\n");
        idx_base = -1;
        strcpy(inpfile, argv[i + 3]);
        strcpy(outfile, argv[i + 4]);
    } else
        std_base_usage();

    if (*fit_std)
        return;

    if (idx_base >= 0) {
        *c1_n_org = pars_used[idx_base][0];
        *dist_n2org = pars_used[idx_base][1];
        *n_org_x = pars_used[idx_base][2];
    }

    fprintf(stderr, "three parameters used: %8.2f%8.2f%8.2f\n",
            *c1_n_org, *dist_n2org, *n_org_x);
}
