#include "x3dna.h"

int main(int argc, char *argv[])
{
    long i, j, k, istart = 1, istep = 1, icnt = 0, waters = 0;
    long bz = 1;  /* for B-Z junction */
    time_t time0;

    time(&time0);

    /* clean up files with fixed name */
    remove_file(AUX_FILE);
    remove_file(BPSTEP_FILE);
    remove_file(HLXSTEP_FILE);
    remove_file(STACK_FILE);
    remove_file(HSTACK_FILE);
    remove_file(REF_FILE);
    remove_file(POC_FILE);
    remove_file(SEVEN_FILE);

    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'C')
                icnt = 1;
            else if (argv[i][j] == 'W')
                waters = 1;
            else if (argv[i][j] == 'N')
                bz = 0;
            else if (argv[i][j] == 'S') {
                k = sscanf(argv[i], "-S=%ld,%ld", &istep, &istart);
                if (k == 2) {
                    if (istart < 0)
                        istart = -istart;  /* make it positive */
                } else if (k == 1)
                    istart = 1;  /* redundant: already initialized */
                else {
                    fprintf(stderr, "wrong format for setting step\n");
                    analyze_usage();
                }
                fprintf(stderr, "***start at %ld, with step size: %ld***\n", istart,
                        istep);
                break;  /* -s=istep,istart not combined with others */
            } else
                analyze_usage();
    }
    if (argc == i)
        process_str("stdin", icnt, waters, istart, istep, bz);
    else
        for (j = i; j < argc; j++) {
            if (strcmp(argv[j], "tmpfile"))  /* manalyze */
                fprintf(stderr, "\n......Processing structure #%ld: <%s>......\n",
                        j - i + 1, argv[j]);
            process_str(argv[j], icnt, waters, istart, istep, bz);
        }

    print_used_time(time0);

    return 0;
}

void analyze_usage(void)
{
    help3dna("analyze");
    contact_msg(0);
}

void process_str(char *inpfile, long icnt, long waters, long istart, long istep, long bz)
/* process a structure */
{
    char pdbfile[BUF512], outfile[BUF512], HB_ATOM[BUF512], ALT_LIST[BUF512];
    char *ChainID, *bseq, **AtomName, **bp_seq, **ResName, **Miscs;

    double mtwist, HB_UPPER[2];
    double *mst_org, *mst_orgH, *mst_orien, *mst_orienH;
    double **org, **orien, **twist_rise, **xyz;
    double **nt_bb_torsion;

    long i, ip, bbexist = 0, ds, hetatm = 0, parallel = 0, nbpm1, num;
    long num_bp, num_residue, str_type = 0;  /* for duplex only */
    long *idx, *WC_info, *ResSeq, *bphlx, *RY;
    long **c6_c8, **chi, **pair_num, **phos;
    long **seidx, **sugar, **o3p_brk, **htm_water;

    FILE *fp;

    /* get PDB file and pairing information */
    pair_num = read_input(inpfile, pdbfile, outfile, &ds, &num_bp, &ip, &hetatm);

    bphlx = lvector(1, num_bp);  /* helix break marker */
    for (i = 1; i <= num_bp; i++)
        bphlx[i] = icnt ? 0 : pair_num[ds + 1][i];

    fp = open_file(outfile, "w");

    hb_crt_alt(HB_UPPER, HB_ATOM, ALT_LIST);

    /* read in PDB file */
    num = number_of_atoms(pdbfile, hetatm, ALT_LIST);
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(pdbfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
             hetatm, ALT_LIST);

    /* get the numbering information of each residue */
    seidx = residue_idx(num, ResSeq, Miscs, ChainID, ResName, &num_residue);

    /* base-pairing residue number checking */
    pair_checking(ip, ds, num_residue, pdbfile, &num_bp, pair_num);

    /* check strand direction: O3'-P linkage, parallel etc */
    o3p_brk = lmatrix(1, ds, 1, num_bp);
    drct_checking(ds, num_bp, pair_num, seidx, AtomName, xyz, &parallel,
                  &bbexist, o3p_brk, fp);

    /* get base or base-pair sequence */
    bp_seq = cmatrix(0, ds, 1, num_bp);  /* need to be kept */
    RY = lvector(1, num_residue);  /* simplified */
    get_bpseq(ds, num_bp, pair_num, seidx, AtomName, ResName, ChainID,
              ResSeq, Miscs, xyz, bp_seq, RY);

    /* for base overlap */
    bseq = cvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);  /* RY re-calculated */

    print_header(num_bp, num, pdbfile, fp);

    /* get atom list for each residue */
    phos = lmatrix(1, ds + 4, 1, num_bp);  /* also include O1P & O2P atoms */
    c6_c8 = lmatrix(1, ds, 1, num_bp);
    sugar = lmatrix(1, ds, 1, num_bp * 5);
    chi = lmatrix(1, ds, 1, num_bp * 4);
    atom_list(ds, num_bp, pair_num, seidx, RY, bp_seq, AtomName, ResName,
              ChainID, ResSeq, Miscs, phos, c6_c8, sugar, chi);

    nt_bb_torsion = dmatrix(1, num_residue, 1, 6);
    get_nt_bb_torsion(nt_bb_torsion, num_residue, seidx, RY, AtomName, ResName,
                      ChainID, ResSeq, Miscs, xyz);

    /* get the reference frame for each base */
    orien = dmatrix(1, ds, 1, num_bp * 9);
    org = dmatrix(1, ds, 1, num_bp * 3);
    WC_info = lvector(1, num_bp);
    ref_frames(ds, num_bp, pair_num, bp_seq, seidx, RY, AtomName, ResName, ChainID,
               ResSeq, Miscs, xyz, fp, orien, org, WC_info, &str_type, 0, o3p_brk);

    if (ds == 2)  /* H-Bond information */
        hb_information(num_bp, pair_num, bp_seq, seidx, AtomName, xyz, WC_info, fp);

    base_overlap(ds, num_bp, num, num_residue, pair_num, RY, bp_seq, seidx, AtomName, xyz, orien, org, fp);  /* base overlap area in A^2 */

    /* get and print out the parameters */
    nbpm1 = num_bp - 1;
    mst_orien = dvector(1, nbpm1 * 9);
    mst_org = dvector(1, nbpm1 * 3);
    mst_orienH = dvector(1, nbpm1 * 9);
    mst_orgH = dvector(1, nbpm1 * 3);
    twist_rise = dmatrix(1, nbpm1, 1, 2);
    get_parameters(ds, num_bp, bp_seq, orien, org, fp, twist_rise, mst_orien,
                   mst_org, mst_orienH, mst_orgH, bphlx, istart, istep, bz, &str_type);

    if (ds == 2) {
        mtwist = get_mtwist(nbpm1, bphlx, WC_info, twist_rise);

        str_classify(mtwist, str_type, parallel, num_bp, fp);
        lambda_d3(num_bp, bp_seq, chi, c6_c8, xyz, fp);

        idx = lvector(1, num);
        atom_idx(num, AtomName, NULL, idx);

        htm_water = lmatrix(1, 4, 0, num);  /* HETATM and water index */
        init_htm_water(waters, num, num_residue, idx, htm_water);
        identify_htw(num_residue, seidx, RY, AtomName, ResName, ChainID, ResSeq,
                     Miscs, xyz, htm_water);

        write_mst(num_bp, pair_num, bp_seq, mst_orien, mst_org, seidx, AtomName,
                  ResName, ChainID, ResSeq, xyz, Miscs, htm_water, STACK_FILE);
        write_mst(num_bp, pair_num, bp_seq, mst_orienH, mst_orgH, seidx, AtomName,
                  ResName, ChainID, ResSeq, xyz, Miscs, htm_water, HSTACK_FILE);

        free_lvector(idx, 1, num);
        free_lmatrix(htm_water, 1, 4, 0, num);

        if (bbexist) {
            print_PP(parallel, mtwist, twist_rise, num_bp, bp_seq, phos, mst_orien,
                     mst_org, mst_orienH, mst_orgH, xyz, WC_info, bphlx, fp);
            groove_width(parallel, num_bp, bp_seq, phos, xyz, bphlx, fp);
        }
        other_pars(num_bp, bp_seq, bphlx, orien, org);
    }

    /* global analysis */
    global_analysis(ds, num_bp, num, bp_seq, chi, phos, xyz, fp);

    if (bbexist) {
        backbone_torsion(ds, num_bp, pair_num, bp_seq, sugar, chi, xyz,
                         nt_bb_torsion, fp);
        p_c1_dist(ds, num_bp, bp_seq, phos, chi, xyz, bphlx, fp);
        helix_radius(ds, num_bp, bp_seq, orien, org, phos, chi, xyz, bphlx, fp);
    }
    get_helix_axis(ds, num_bp, bp_seq, orien, org, bphlx, fp);

    /* free allocated vectors & matrices */
    free_lmatrix(pair_num, 1, ds + 1, 1, num_bp);
    free_lvector(bphlx, 1, num_bp);
    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_lmatrix(o3p_brk, 1, ds, 1, num_bp);
    free_cmatrix(bp_seq, 0, ds, 1, num_bp);
    free_lvector(RY, 1, num_residue);
    free_cvector(bseq, 1, num_residue);
    free_lmatrix(phos, 1, ds + 4, 1, num_bp);
    free_lmatrix(c6_c8, 1, ds, 1, num_bp);
    free_dmatrix(nt_bb_torsion, 1, num_residue, 1, 6);
    free_lmatrix(sugar, 1, ds, 1, num_bp * 5);
    free_lmatrix(chi, 1, ds, 1, num_bp * 4);
    free_dmatrix(orien, 1, ds, 1, num_bp * 9);
    free_dmatrix(org, 1, ds, 1, num_bp * 3);
    free_lvector(WC_info, 1, num_bp);
    free_dvector(mst_orien, 1, nbpm1 * 9);
    free_dvector(mst_org, 1, nbpm1 * 3);
    free_dvector(mst_orienH, 1, nbpm1 * 9);
    free_dvector(mst_orgH, 1, nbpm1 * 3);
    free_dmatrix(twist_rise, 1, nbpm1, 1, 2);

    close_file(fp);
}
