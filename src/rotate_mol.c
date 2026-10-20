#include "x3dna.h"

int main(int argc, char *argv[])
/* rotate a structure [PDB or ALCHEMY] to adjust its orientation
   ---------------------------------------------------------------
   rotating or translating a PDB file in RasMol and then write it
   in PDB format does NOT change its coordinates. One way around
   is to "write molscript MOL-FILE", and then extract the rotation
   information from it, which looks like this:

   by rotation x 180.0
   by rotation z 105.669
   by rotation y -9.26690
   by rotation x 92.6453;

   the rotation matrix used here will be (note the negative values):
   rotx(-180)*rotz(-105.669)......

   to have a specific control over the orientation of a structure,
   either in PDB or ALCHEMY format, user may also need to create
   such a file by hand (with "by rotation" prefix).

   fnum = 0 means neither transformation matrix or rotation file
        = 1 means rotation file as from RasMol
        = 2 means transformation matrix
        = 3 means to read in the reference frame for aligning a structure TO
   ---------------------------------------------------------------
 */
{
    char rotfile[BUF512] = "", inpfile[BUF512], outfile[BUF512];
    double **rotmat;
    long num, fnum = 0, to_center = 0, pmi = 2;  /* pmi of base residues */
    long to_frame = 0;  /* set for aligning the structure TO the reference frame */

    rotate_mol_cmdline(argc, argv, rotfile, inpfile, outfile, &fnum, &to_center, &pmi,
                       &to_frame);
    if (fnum) {
        to_center = 0;
        pmi = 0;
        if (strcmp(rotfile, ROTMAT_FILE))
            remove_file(ROTMAT_FILE);
    } else
        remove_file(ROTMAT_FILE);
    remove_file(VIEW1_FILE);
    remove_file(VIEW2_FILE);
    remove_file(VIEW3_FILE);

    rotmat = dmatrix(1, 4, 1, 3);
    identity_matrix(rotmat, 3);
    if (fnum == 1)
        molscr_rotmat(rotfile, rotmat);
    else if (fnum == 2)
        read_rotmat(rotfile, rotmat);
    else if (fnum == 3)
        read_bp1frame(rotfile, rotmat);

    num = number_of_atoms(inpfile, 1, "*");
    (num) ? rot_pdbfile(num, inpfile, outfile, to_center, rotmat, pmi, to_frame) :
        rot_alcfile(inpfile, outfile, to_center, rotmat, pmi, to_frame);

    free_dmatrix(rotmat, 1, 4, 1, 3);

    return 0;
}

void rotmol_usage(void)
{
    help3dna("rotate_mol");
    contact_msg(0);
}

void rotate_mol_cmdline(long argc, char *argv[], char *rotfile, char *inpfile,
                        char *outfile, long *fnum, long *to_center, long *pmi,
                        long *to_frame)
{
    char str[BUF512];  /* keep a copy of the original */
    long i, j;

    if (argc < 3)
        rotmol_usage();
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        strcpy(str, argv[i]);
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++) {  /* skip - */
            if (argv[i][j] == 'C')
                *to_center = 1;
            else if (argv[i][j] == 'A')
                *pmi = 1;  /* all atoms */
            else if (argv[i][j] == 'B')
                *pmi = 2;  /* base */
            else if (argv[i][j] == 'P')
                *pmi = 3;  /* protein */
            else if (argv[i][j] == 'L')
                *to_frame = 1;  /* to align along the reference frame */
            else if (argv[i][j] == 'R') {
                if (sscanf(str, "-R=%s", rotfile) != 1 &&
                    sscanf(str, "-r=%s", rotfile) != 1) {
                    fprintf(stderr, "wrong format for reading rotation angles\n");
                    rotmol_usage();
                } else {
                    *fnum = 1;
                    break;  /* -r=transMtx file not combined with others */
                }
            } else if (argv[i][j] == 'T') {
                if (sscanf(str, "-T=%s", rotfile) != 1 &&
                    sscanf(str, "-t=%s", rotfile) != 1) {
                    fprintf(stderr, "wrong format for reading transformation matrix\n");
                    rotmol_usage();
                } else {
                    *fnum = 2;
                    break;  /* -t=rotfile not combined with others */
                }
            } else if (argv[i][j] == 'O') {  /* to covert rebuilt str to the original frame */
                if (sscanf(str, "-O=%s", rotfile) != 1 &&
                    sscanf(str, "-o=%s", rotfile) != 1) {
                    fprintf(stderr, "wrong format for reading reference frame file\n");
                    rotmol_usage();
                } else {
                    *fnum = 3;
                    *to_frame = 1;
                    break;  /* -o=bp1_reference_frame not combined with others */
                }
            } else
                rotmol_usage();
        }
    }
    if (argc == i + 2) {
        strcpy(inpfile, argv[i++]);
        strcpy(outfile, argv[i]);
    } else
        rotmol_usage();
    if (*to_frame && *fnum < 2)
        fprintf(stderr, "ignoring -l option: it must be used with -t=TRNSMAT\n");
}

void molscr_rotmat(char *rotfile, double **rotmat)
/* get rotation matrix based on MolScript input (as from RasMol) */
{
    char axis, str[BUF512], str_by[BUF512], str_rotation[BUF512];
    double rot_ang, **roti, **temp;
    long naxis = 0;
    FILE *fp;

    fp = open_file(rotfile, "r");
    roti = dmatrix(1, 3, 1, 3);
    temp = dmatrix(1, 3, 1, 3);
    while (fgets(str, sizeof str, fp) != NULL) {
        upperstr(str);  /* change to UPPER case */
        if (sscanf(str, "%s %s", str_by, str_rotation) == 2 &&
            !strcmp(str_by, "BY") && !strcmp(str_rotation, "ROTATION")) {
            if (sscanf(str, "%*s %*s %c %lf", &axis, &rot_ang) != 2)
                fatal("wrong <by rotation> format\n");
            if (axis == 'X' || axis == 'Y' || axis == 'Z')
                naxis++;
            else
                fatal("none of x-, y- or z-axis\n");
            if (axis == 'X')
                rotx(-rot_ang, roti);
            else if (axis == 'Y')
                roty(-rot_ang, roti);
            else  /* axis == 'Z' */
                rotz(-rot_ang, roti);
            multi_matrix(rotmat, 3, 3, roti, 3, 3, temp);
            copy_matrix(temp, 3, 3, rotmat);
        }
    }
    close_file(fp);
    if (!naxis)
        fprintf(stderr, "NO rotation applied to the structure\n");

    free_dmatrix(roti, 1, 3, 1, 3);
    free_dmatrix(temp, 1, 3, 1, 3);
}

void read_rotmat(char *rotfile, double **rotmat)
/* read in a rotation matrix: x-axis can be any vector;
   y-axis orthogonal to it; z-axis forms a right-handed system */
{
    char str[BUF512], *format = "%lf %lf %lf";
    double tmp[4], **tmprot;
    long i, row_wise;
    FILE *fp;

    fp = open_file(rotfile, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", &row_wise) != 1)
        fatal("error reading rotmat file: row/column-wise indicator\n");
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, format, &rotmat[4][1], &rotmat[4][2], &rotmat[4][3]) != 3)
        fatal("error reading origin line\n");

    tmprot = dmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++) {
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &tmp[1], &tmp[2], &tmp[3]) != 3)
            fatal("error reading x/y/z directions\n");
        vec_norm(tmp);  /* normalize it */
        cpxyz(tmp, tmprot[i]);
    }
    close_file(fp);

    if (!row_wise) {  /* make tmprot row-vector for xyz axes */
        transpose_matrix(tmprot, 3, 3, rotmat);
        copy_matrix(rotmat, 3, 3, tmprot);
    }
    if (fabs(dot(tmprot[1], tmprot[2])) >= 0.005)
        fprintf(stderr, "Warning: x- and y-axes are not orthogonal\n");
    vec_orth(tmprot[2], tmprot[1]);  /* orthogonal component to x-axis as y-axis */
    cross(tmprot[1], tmprot[2], tmp);
    if (fabs(dot(tmp, tmprot[3])) <= 0.9995)
        fprintf(stderr, "Warning: z-axis is not perpendicular to xy-plane\n");
    cpxyz(tmp, tmprot[3]);
    transpose_matrix(tmprot, 3, 3, rotmat);  /* make rotmat column-wise */

    free_dmatrix(tmprot, 1, 3, 1, 3);
}

void read_bp1frame(char *rotfile, double **rotmat)
/* read in the reference frame of the 1st base-pair. this is used to
 * re-orient 3DNA rebuilt structure to the original coordinate frame.
 * the reference frame file can be obtained by either "find_pair" or
 * "analyze". It needs to be renamed to a difference name (e.g.,
 * bp1_frame.dat) from the system default "ref_frames.dat". see also
 * function "read_frames" in file frame_mol.c */
{
    char *format = "%lf %lf %lf", str[BUF512];
    long num_bp, inum = 0;
    FILE *fp;

    fp = open_file(rotfile, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", &num_bp) != 1)
        fatal("error reading ref-frame file: number of base-pairs\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%*s %ld", &inum) != 1)
        fatal("error reading title line\n");

    if (inum != 1)
        fprintf(stderr, "Warning: <%ld> is not the 1st base-pair!\n", inum);

    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, format, &rotmat[4][1], &rotmat[4][2], &rotmat[4][3]) != 3)
        fatal("error reading origin line\n");

    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, format, &rotmat[1][1], &rotmat[2][1], &rotmat[3][1]) != 3)
        fatal("error reading x direction cosines\n");
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, format, &rotmat[1][2], &rotmat[2][2], &rotmat[3][2]) != 3)
        fatal("error reading y direction cosines\n");
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, format, &rotmat[1][3], &rotmat[2][3], &rotmat[3][3]) != 3)
        fatal("error reading z direction cosines\n");

    close_file(fp);
}

void rot_pdbfile(long num, char *inpfile, char *outfile, long to_center, double **rotmat,
                 long pmi, long to_frame)
{
    char *ChainID, **AtomName, **ResName, **Miscs;
    double **nxyz, **xyz, **refmat, **tmprot;
    long nnum, *ResSeq;

    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    nxyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");

    if (to_frame) {
        align_str(rotmat, num, xyz, nxyz);
        goto CLEAN_UP_PDB;
    }

    refmat = dmatrix(1, 4, 1, 3);
    identity_matrix(refmat, 3);  /* initialize it to identity matrix */
    if (pmi)
        pdb_pmi(num, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, pmi, nxyz,
                refmat, rotmat, &nnum);

    if (to_center)
        (pmi >= 2 && nnum) ?  /* nucleic acid base or protein */
            ave_dmatrix(nxyz, nnum, 3, rotmat[4]) : ave_dmatrix(xyz, num, 3, rotmat[4]);  /* all atoms */
    move_position(xyz, num, 3, rotmat[4]);

    if (pmi) {  /* three views */
        tmprot = dmatrix(1, 3, 1, 3);
        multi_matrix(xyz, num, 3, rotmat, 3, 3, nxyz);  /* default: most extended */
        write_pdb(num, AtomName, ResName, ChainID, ResSeq, nxyz, Miscs, VIEW1_FILE);

        pmi_view(rotmat, tmprot, 3, 1, 2);  /* top view */
        multi_matrix(xyz, num, 3, tmprot, 3, 3, nxyz);
        write_pdb(num, AtomName, ResName, ChainID, ResSeq, nxyz, Miscs, VIEW2_FILE);

        pmi_view(rotmat, tmprot, 2, 3, 1);  /* landscape view */
        multi_matrix(xyz, num, 3, tmprot, 3, 3, nxyz);
        write_pdb(num, AtomName, ResName, ChainID, ResSeq, nxyz, Miscs, VIEW3_FILE);
        free_dmatrix(tmprot, 1, 3, 1, 3);
    }

    multi_matrix(xyz, num, 3, rotmat, 3, 3, nxyz);
    cmb_transformation(refmat, rotmat);
    free_dmatrix(refmat, 1, 4, 1, 3);

  CLEAN_UP_PDB:
    write_pdb(num, AtomName, ResName, ChainID, ResSeq, nxyz, Miscs, outfile);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_dmatrix(nxyz, 1, num, 1, 3);

}

void pdb_pmi(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
             double **xyz, char **Miscs, long pmi, double **nxyz, double **refmat,
             double **rotmat, long *nnum)
/* get rotation matrix based on PMI of a PDB file */
{
    char BDIR[BUF512], *bseq;
    double **orien, **org;
    long i, j, num_residue, nnum2 = 0, *RY, **seidx;

    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);
    bseq = cvector(1, num_residue);
    RY = lvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);

    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    get_BDIR(BDIR, "Atomic_A.pdb", 1);
    base_frame(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq, Miscs,
               xyz, BDIR, orien, org);

    for (i = 1; i <= num_residue; i++)
        if (RY[i] >= 0) {
            reorient_xyz(num, xyz, orien[i], org[i], refmat);
            break;
        }
    if (i > num_residue)  /* no base residue */
        reorient_xyz(num, xyz, NULL, NULL, refmat);
    for (i = 1; i <= num_residue; i++)
        if ((pmi == 2 && RY[i] >= 0) ||  /* base residue */
            (pmi == 3 && RY[i] == -1))  /* protein */
            for (j = seidx[i][1]; j <= seidx[i][2]; j++) {
                nnum2++;
                cpxyz(xyz[j], nxyz[nnum2]);
            }
    if (pmi >= 2 && nnum2) {
        pmi_rotmat(nnum2, nxyz, rotmat);
        multi_matrix(xyz, num, 3, rotmat, 3, 3, nxyz);
        if (pmi == 2) {  /* nucleic acid bases */
            base_frame(num_residue, bseq, seidx, RY, AtomName, ResName,
                       ChainID, ResSeq, Miscs, nxyz, BDIR, orien, org);
            base1_check(num_residue, RY, org, rotmat);
        }
    } else
        pmi_rotmat(num, xyz, rotmat);  /* all atoms used */
    *nnum = nnum2;

    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lvector(RY, 1, num_residue);
    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);
}

void rot_alcfile(char *inpfile, char *outfile, long to_center, double **rotmat, long pmi,
                 long to_frame)
{
    char **AtomName;
    double **nxyz, **xyz, **refmat, **tmprot;
    long nbond, num, *ibase, **linkage;

    get_alc_nums(inpfile, &num, &nbond);
    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 1, 2);
    nxyz = dmatrix(1, num, 1, 3);

    read_alc(inpfile, &num, &nbond, AtomName, xyz, ibase, linkage);

    if (to_frame) {
        align_str(rotmat, num, xyz, nxyz);
        goto CLEAN_UP_ALC;
    }

    if (to_center)
        ave_dmatrix(xyz, num, 3, rotmat[4]);
    move_position(xyz, num, 3, rotmat[4]);

    refmat = dmatrix(1, 4, 1, 3);
    identity_matrix(refmat, 3);  /* initialize it to identity matrix */

    tmprot = dmatrix(1, 3, 1, 3);
    if (pmi) {  /* use all xyz coordinates */
        reorient_xyz(num, xyz, NULL, NULL, refmat);

        pmi_rotmat(num, xyz, rotmat);  /* default: most extended */
        multi_matrix(xyz, num, 3, rotmat, 3, 3, nxyz);
        write_alc(num, nbond, AtomName, nxyz, ibase, linkage, VIEW1_FILE);

        pmi_view(rotmat, tmprot, 3, 1, 2);  /* top view: rotx(-90) * rotz(-90) */
        multi_matrix(xyz, num, 3, tmprot, 3, 3, nxyz);
        write_alc(num, nbond, AtomName, nxyz, ibase, linkage, VIEW3_FILE);

        pmi_view(rotmat, tmprot, 2, 3, 1);  /* landscape view */
        multi_matrix(xyz, num, 3, tmprot, 3, 3, nxyz);
        write_alc(num, nbond, AtomName, nxyz, ibase, linkage, VIEW2_FILE);
    }
    free_dmatrix(tmprot, 1, 3, 1, 3);

    cmb_transformation(refmat, rotmat);
    free_dmatrix(refmat, 1, 4, 1, 3);

    multi_matrix(xyz, num, 3, rotmat, 3, 3, nxyz);

  CLEAN_UP_ALC:
    write_alc(num, nbond, AtomName, nxyz, ibase, linkage, outfile);

    free_alc(num, nbond, AtomName, xyz, ibase, 1, linkage);
    free_dmatrix(nxyz, 1, num, 1, 3);
}

void cmb_transformation(double **refmat, double **rotmat)
/* combine two transformation matrices: applied to PMI view
 * [(xyz - o1) * r1 - o2] * r2
 *   ====> [(xyz - o1) * r1 - o2] * (inv(r1) * r1 * r2)
 *   ====> [xyz - (o1 + o2 * inv(r1))] * (r1 * r2)
 *  O = o1 + o2 * inv(r1); R = r1 * r2
 * inv(r1) = transpose of r1
 * Here r1: refmat; r2: rotmat */
{
    double **cmbrot;

    cmbrot = dmatrix(1, 4, 1, 3);

    multi_vec_Tmatrix(rotmat[4], 3, refmat, 3, 3, cmbrot[4]);
    sumxyz(refmat[4], cmbrot[4], cmbrot[4]);
    multi_matrix(refmat, 3, 3, rotmat, 3, 3, cmbrot);
    write_rotmat(cmbrot);

    free_dmatrix(cmbrot, 1, 3, 1, 3);
}

void write_rotmat(double **rotmat)
/* Write the position and orientation of the transformation to ROTMAT_FILE */
{
    char *format = "%12.4f%12.4f%12.4f\n";
    long i;
    FILE *fp;

    fp = open_file(ROTMAT_FILE, "w");
    fprintf(fp, "    1  # x-, y-, z-axes row-rise\n");
    fprintf(fp, format, rotmat[4][1], rotmat[4][2], rotmat[4][3]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, rotmat[1][i], rotmat[2][i], rotmat[3][i]);

    close_file(fp);
}

void base1_check(long num_residue, long *RY, double **org, double **rotmat)
/* Double check the positioning of base 1 origin: it should be in
 * lower-left compared to the rest. further checking: 1DGC */
{
    double b1[4], bmax[4], bmin[4], **tmpmtx, **rot180;
    long i, j, nb = 0;

    for (i = 1; i <= 3; i++) {
        b1[i] = 0.0;
        bmax[i] = -XBIG;
        bmin[i] = XBIG;
    }
    for (i = 1; i <= num_residue; i++)
        if (RY[i] >= 0) {
            nb++;
            for (j = 1; j <= 3; j++) {
                if (nb == 1)
                    b1[j] = org[i][j];
                bmax[j] = dval_max(bmax[j], org[i][j]);
                bmin[j] = dval_min(bmin[j], org[i][j]);
            }
        }

    tmpmtx = dmatrix(1, 3, 1, 3);
    copy_matrix(rotmat, 3, 3, tmpmtx);

    rot180 = dmatrix(1, 3, 1, 3);
    if ((b1[1] - bmin[1]) > 0.5 * (bmax[1] - bmin[1])) {  /* x-axis reversed */
        roty(180.0, rot180);
        multi_matrix(tmpmtx, 3, 3, rot180, 3, 3, rotmat);
        copy_matrix(rotmat, 3, 3, tmpmtx);  /* copy again! */
    }
    if ((b1[2] - bmin[2]) > 0.5 * (bmax[2] - bmin[2])) {  /* y-axis reversed */
        rotx(180.0, rot180);
        multi_matrix(tmpmtx, 3, 3, rot180, 3, 3, rotmat);
    }

    free_dmatrix(tmpmtx, 1, 3, 1, 3);
    free_dmatrix(rot180, 1, 3, 1, 3);
}

void pmi_rotmat(long num, double **xyz, double **rotmat)
/* Get the principle moment of inertia. By default, x-axis points
   along the normal to the best plane (least variance), z-axis along
   the best line (most variance), and y-axis in between. Reorder it as
   [2 3 1] to get the most-extended view and make sure a(?)-axis has
   positive components: 1ber */
{
    double tx[4], ty[4], tz[4], **cov_mtx, **rot180;
    long i;

    if (!num)
        return;
    cov_mtx = dmatrix(1, 3, 1, 3);
    cov_matrix(xyz, xyz, num, 3, cov_mtx);
    jacobi(cov_mtx, 3, tx, rotmat);

    /* Make sure (:, 3) == cross((:, 1),  (:, 2)) */
    for (i = 1; i <= 3; i++) {
        tx[i] = rotmat[i][1];
        ty[i] = rotmat[i][2];
    }
    cross(tx, ty, tz);
    for (i = 1; i <= 3; i++)
        rotmat[i][3] = tz[i];

    pmi_view(rotmat, cov_mtx, 2, 3, 1);  /* 2-3-1 view */
    copy_matrix(cov_mtx, 3, 3, rotmat);

    rot180 = dmatrix(1, 3, 1, 3);
    if (cov_mtx[1][1] < 0.0) {  /* x-axis reversed */
        roty(180.0, rot180);
        multi_matrix(cov_mtx, 3, 3, rot180, 3, 3, rotmat);
        copy_matrix(rotmat, 3, 3, cov_mtx);  /* copy back */
    }
    if (cov_mtx[2][2] < 0.0) {  /* y-axis reversed */
        rotx(180.0, rot180);
        multi_matrix(cov_mtx, 3, 3, rot180, 3, 3, rotmat);
    }

    free_dmatrix(cov_mtx, 1, 3, 1, 3);
    free_dmatrix(rot180, 1, 3, 1, 3);
}

void pmi_view(double **rotmat, double **tmprot, long c1, long c2, long c3)
/* get the new view of PMI */
{
    long i;

    for (i = 1; i <= 3; i++) {
        tmprot[i][1] = rotmat[i][c1];  /* x-axis */
        tmprot[i][2] = rotmat[i][c2];  /* y-axis */
        tmprot[i][3] = rotmat[i][c3];  /* z-axis */
    }
}

void reorient_xyz(long num, double **xyz, double *orien1, double *org1, double **refmat)
/* Reorient the structure with regard to the first base reference frame if any.
   Otherwise, use the first three atoms to define the reference frame */
{
    double d, tmpx[4], tmpy[4], tmpz[4];
    long i;

    if (num < 2)
        return;

    if (orien1 == NULL) {  /* no base reference */
        cpxyz(xyz[1], refmat[4]);  /* 1st atom as origin */
        ddxyz(xyz[1], xyz[2], tmpx);  /* 1 ---> 2 as x-axis */
        for (i = 3; i <= num; i++) {
            ddxyz(xyz[1], xyz[i], tmpy);
            d = magang(tmpx, tmpy);
            if (dval_in_range(d, 5, 175))  /* nonlinear */
                break;
        }
        if (i > num) {
            fprintf(stderr, "too few atoms to change orientation\n");
            return;
        }
        vec_norm(tmpx);
        vec_orth(tmpy, tmpx);  /* orthogonal component as y-axis */
        cross(tmpx, tmpy, tmpz);
        x_y_z_2_mtx(tmpx, tmpy, tmpz, refmat);  /* column vector for x, y & z */
    } else {  /* with base reference */
        orien2mst(orien1, 0, refmat);
        cpxyz(org1, refmat[4]);
    }

    /* ignore translations */
    init_dvector(refmat[4], 1, 3, 0.0);

    change_xyz(0, refmat[4], refmat, num, xyz);
}

void align_str(double **rotmat, long num, double **xyz, double **nxyz)
/* align the structure so that its global reference frame is with
 * reference to that defined in rotmat. note that "rotmat" is
 * column-wise for x-, y- and z-axes, and the origin is in rotmat[4] */
{
    long i, j;

    for (i = 1; i <= num; i++)
        for (j = 1; j <= 3; j++)
            nxyz[i][j] = dot(xyz[i], rotmat[j]) + rotmat[4][j];
}
