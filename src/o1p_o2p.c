#include "x3dna.h"

int main(int argc, char *argv[])
/* There is a confusion about O1P/O2P labeling in nucleic acid
   structures which affects least-fitting RMS value between two
   structures. This utility checks if O1P/O2P in a PDB file are
   properly labeled as in the NDB. A new clean PDB file is written. */
{
    char idmsg[BUF512], *ChainID, **AtomName, **ResName, **Miscs;
    double d, O1P_O2P[4], O2P_O5[4], O_P[4], temp[4], **xyz;
    long i, ib, ie, inum = 0, j, num, num_residue, O1P, O2P, O5, P;
    long *ResSeq, **seidx;
    FILE *fp;

    if (argc != 3) {
        help3dna("o1p_o2p");
        contact_msg(0);
    }

    /* read in PDB file */
    num = number_of_atoms(argv[1], 1, "*");
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(argv[1], NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");

    /* get the numbering information of each residue */
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

    /* check if O1P/O2P labeling are okay */
    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];

        if (num_strmatch(" P  ", AtomName, ib, ie)) {  /* with P atom */
            sprintf(idmsg, ": residue name %s, chain %c, number %4ld%c",
                    ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2]);
            O1P = find_1st_atom(" O1P", AtomName, ib, ie, idmsg);
            O2P = find_1st_atom(" O2P", AtomName, ib, ie, idmsg);
            O5 = find_1st_atom(" O5'", AtomName, ib, ie, idmsg);
            P = find_1st_atom(" P  ", AtomName, ib, ie, idmsg);

            if (O1P && O2P && O5 && P) {
                ddxyz(xyz[O1P], xyz[O2P], O1P_O2P);
                ddxyz(xyz[O2P], xyz[O5], O2P_O5);
                for (j = 1; j <= 3; j++) {
                    d = (xyz[O1P][j] + xyz[O2P][j] + xyz[O5][j]) / 3.0;
                    O_P[j] = xyz[P][j] - d;
                }
                cross(O1P_O2P, O2P_O5, temp);
                if (dot(temp, O_P) < 0) {
                    fprintf(stderr,
                            "*PO4 group at residue #%ld %s on chain %c"
                            " does NOT conform to the NDB convention*\n",
                            i, ResName[ib], ChainID[ib]);
                    for (j = 1; j <= 3; j++)
                        dval_swap(&xyz[O1P][j], &xyz[O2P][j]);
                }
            }
        }
    }
    fp = open_file(argv[2], "w");
    print_pdb_title(argv[1], "*", fp);
    pdb_record(1, num, &inum, 0, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, fp);
    fprintf(fp, "END\n");
    close_file(fp);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);

    return 0;
}
