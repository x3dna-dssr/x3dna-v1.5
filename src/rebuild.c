#include "x3dna.h"

int main(int argc, char *argv[])
{
    char BDIR[BUF512], inpfile[BUF512], outfile[BUF512];
    char user_input[BUF512];
    char *options[] = { "-ATOMIC", "-BASE_P", "-BLOCK1", "-BLOCK2" };
    char **bp_seq;
    double **bp_par, **step_par;
    long neg_xdir = 0, is_single = 0, is_helical = 0, xdir = 0;
    long parallel = 0, num_bp, num_atoms;
    long *pidx;
    time_t time0;

    time(&time0);

    remove_file(REF_FILE);

    strcpy(user_input, options[2]);  /* block1 */
    rebuild_cmdline(argc, argv, options, user_input, inpfile, outfile, &neg_xdir);

    num_bp = get_num_bp(inpfile);

    bp_seq = cmatrix(1, num_bp, 0, 2);
    bp_par = dmatrix(1, num_bp, 1, 6);
    step_par = dmatrix(1, num_bp, 1, 6);
    pidx = lvector(1, num_bp);

    read_par(inpfile, &num_bp, &is_single, &is_helical, &xdir, &parallel, bp_seq,
             bp_par, step_par, pidx);
    if (xdir && !neg_xdir)
        fprintf(stderr, "Twist < 0.0, taken as Z-DNA and reverse x- and z-axes\n");
    if (neg_xdir)
        xdir = !xdir;

    if (!strcmp(user_input, options[0])) {
        get_BDIR(BDIR, "Atomic_A.pdb", 1);
        num_atoms = num_PDB_atoms(num_bp, is_single, bp_seq, BDIR);
        if (is_single)
            atomic_pdb1(num_bp, num_atoms, is_helical, xdir, bp_seq, step_par, BDIR,
                        outfile);
        else
            atomic_pdb2(parallel, num_bp, num_atoms, is_helical, xdir, bp_seq,
                        bp_par, step_par, BDIR, outfile);
    } else if (!strcmp(user_input, options[1])) {
        get_BDIR(BDIR, "Atomic_A.pdb", 1);
        if (is_single)
            fatal("<base_p> option only valid for a duplex\n");
        else {
            num_atoms = num_PDB_atoms(num_bp, is_single, bp_seq, BDIR);
            atomic_base_p(parallel, num_bp, num_atoms, is_helical, xdir, bp_seq,
                          bp_par, step_par, BDIR, pidx, outfile);
        }
    } else if (!strcmp(user_input, options[2])) {
        get_BDIR(BDIR, "Block_R.alc", 1);
        block_alc1(num_bp, is_single, is_helical, xdir, bp_seq, step_par, BDIR, outfile);
    } else {
        if (is_single)
            fatal("only one block per base allowed in a single helix\n");
        else {
            get_BDIR(BDIR, "Block_R.alc", 1);
            block_alc2(num_bp, is_helical, xdir, bp_seq, bp_par, step_par, BDIR, outfile);
        }
    }

    free_cmatrix(bp_seq, 1, num_bp, 0, 2);
    free_dmatrix(bp_par, 1, num_bp, 1, 6);
    free_dmatrix(step_par, 1, num_bp, 1, 6);
    free_lvector(pidx, 1, num_bp);

    print_used_time(time0);

    return 0;
}
