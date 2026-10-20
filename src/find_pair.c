#include "x3dna.h"

int main(int argc, char *argv[])
{
    char inpfile[BUF512], outfile[BUF512];
    long ds = 2, curves = 0, divide = 0, hetatm = 1, pairs = 0;
    long detailed = 0, waters = 0;
    time_t time0;

    time(&time0);

    /* clean up files with fixed name */
    remove_file(MUL_FILE);
    remove_file(ALLP_FILE);
    remove_file(BPORDER_FILE);
    remove_file(BESTP_FILE);
    remove_file(REF_FILE);
    remove_file(MREF_FILE);
    remove_file(HLXREG_FILE);
    remove_file(COLCHN_FILE);
    remove_file(COLHLX_FILE);
    remove_file(MULBP_FILE);
    remove_file(TMP_FILE);

    fp_cmdline(argc, argv, inpfile, outfile, &ds, &curves, &divide, &hetatm, &pairs,
               &detailed, &waters);
    handle_str(inpfile, outfile, ds, curves, divide, hetatm, pairs, detailed, waters);

    print_used_time(time0);

    return 0;
}

void fp_usage(void)
{
    help3dna("find_pair");
    contact_msg(0);
}

void fp_cmdline(int argc, char *argv[], char *inpfile, char *outfile, long *ds,
                long *curves, long *divide, long *hetatm, long *pairs,
                long *detailed, long *waters)
{
    long i, j;

    if (argc < 3)
        fp_usage();
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'S' || argv[i][j] == '1')
                *ds = 1;
            else if (argv[i][j] == 'C')
                *curves = 1;
            else if (argv[i][j] == 'D')
                *divide = 1;
            else if (argv[i][j] == 'P')
                *pairs = 1;
            else if (argv[i][j] == 'T')
                *hetatm = 1;
            else if (argv[i][j] == 'A')
                *hetatm = 0;
            else if (argv[i][j] == 'Z')
                *detailed = 1;
            else if (argv[i][j] == 'W')
                *waters = 1;
            else
                fp_usage();
    }
    if (argc == i + 2) {
        strcpy(inpfile, argv[i]);
        strcpy(outfile, argv[i + 1]);
    } else
        fp_usage();
    if (*pairs) {
        if ((*ds == 1) + *curves + *divide)
            fprintf(stderr, "for -p, ignore other options except for -t\n");
        return;
    }
    if (*ds == 1) {
        if (*curves) {
            fprintf(stderr, "no input to Curves for single strand: -c ignored\n");
            *curves = 0;
        }
        if (*divide) {
            fprintf(stderr, "no dividing necessary for single strand: -d ignored\n");
            *divide = 0;
        }
    }
    if (*waters == 1)
        *hetatm = 1;
}

void handle_str(char *pdbfile, char *outfile, long ds, long curves, long divide,
                long hetatm, long pairs, long detailed, long waters)
{
    char b1[BUF512], parfile[BUF512], HB_ATOM[BUF512], ALT_LIST[BUF512];
    char *ChainID, *bseq;
    char **AtomName, **ResName, **Miscs;
    double HB_UPPER[2], **xyz;
    long i, ir, num, num_bp = 0, num_residue;
    long *ResSeq, *RY, **seidx;
    FILE *fp;

    del_extension(pdbfile, parfile);
    hb_crt_alt(HB_UPPER, HB_ATOM, ALT_LIST);

    /* read in the PDB file */
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

    /* get base sequence, RY identification */
    bseq = cvector(1, num_residue);
    RY = lvector(1, num_residue);
    get_seq(num_residue, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, xyz, bseq, RY);

    if (ds == 1) {  /* single helix */
        fp = open_file(outfile, "w");
        for (i = 1; i <= num_residue; i++)
            if (RY[i] >= 0)
                num_bp++;
        fprintf(fp, "%s\n", pdbfile);
        fprintf(fp, "%s.outs\n", parfile);
        fprintf(fp, "    1      # single helix\n");
        fprintf(fp, "%5ld      # number of bases\n", num_bp);
        fprintf(fp, "    1 %4ld # explicit bp numbering/hetero atoms\n", hetatm);
        for (i = 1; i <= num_residue; i++) {
            ir = seidx[i][1];
            base_str(ChainID[ir], ResSeq[ir], Miscs[ir], ResName[ir], bseq[i], 1, b1);
            if (RY[i] >= 0)
                fprintf(fp, "%5ld      # %s\n", i, b1);
        }
        close_file(fp);
    } else {  /* duplex */
        if (pairs)
            multi_bps(pdbfile, parfile);
        duplex(num, num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
               ResSeq, Miscs, xyz, curves, divide, hetatm, pairs, detailed,
               waters, pdbfile, outfile, parfile);
    }

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lmatrix(seidx, 1, num_residue, 1, 2);
    free_cvector(bseq, 1, num_residue);
    free_lvector(RY, 1, num_residue);
}

void multi_bps(char *pdbfile, char *parfile)
/* output a temporary file for multiplets */
{
    FILE *mulbp;

    mulbp = open_file(MULBP_FILE, "w");
    fprintf(mulbp, "%s\n", pdbfile);
    fprintf(mulbp, "%s.outm\n", parfile);
    close_file(mulbp);
}

void multiplets(long max_ple, long num_residue, long **pair_info, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                double **orien, double **org, long **seidx, char *bseq, long hetatm,
                long **htm_water, FILE * fp)
/* print out multiplets information and write out the coordinates */
{
    char b1[BUF512], pairstr[BUF512], tmp[BUF512];
    double morg[4], morien[10];
    long i, j, jr, k, m, num_match, num_partial;
    long **match_list, **partial_list;
    FILE *mfp, *mulbp, *rframe;

    mfp = open_file(MUL_FILE, "w");
    mulbp = open_file(MULBP_FILE, "a");

    fprintf(fp, "\nSummary of triplets and higher multiplets\n");
    match_list = lmatrix(1, num_residue, 0, NP);
    partial_list = lmatrix(1, num_residue, 0, NP);

    base_compatibility(num_residue, pair_info, &num_match, match_list,
                       &num_partial, partial_list, fp);

    fprintf(mulbp, "%5ld         # number of bases per layer\n", max_ple);
    fprintf(mulbp, "%5ld         # number of layers\n", num_match);
    fprintf(mulbp, "    1 %5ld   # explicit bp numbering/hetero atoms\n", hetatm);

    /* perfect-match set of base associations */
    fprintf(fp, "===================== perfect match =====================\n");
    fprintf(mfp, "REMARK    =========================== perfect match"
            " ===========================\n");

    rframe = open_file(MREF_FILE, "w");
    fprintf(rframe, "%5ld base multiplets\n", num_match + num_partial);
    for (i = 1; i <= num_match; i++) {
        fprintf(fp, "%5ld: #%ld ", i, match_list[i][0]);
        pairstr[0] = '\0';
        fprintf(rframe, "... %5ld ", i);
        for (j = 1; j <= match_list[i][0]; j++) {
            jr = match_list[i][j];
            k = seidx[jr][1];
            base_str(ChainID[k], ResSeq[k], Miscs[k], ResName[k], bseq[jr], 1, b1);
            sprintf(tmp, "[%ld]%s%s", jr, b1, (j == match_list[i][0]) ? "" : " + ");
            strcat(pairstr, tmp);
            fprintf(mulbp, "%5ld", jr);
            fprintf(rframe, "%c", bseq[jr]);
        }
        fprintf(fp, "%s\n", pairstr);
        fprintf(mulbp, "\n");
        fprintf(rframe, " ...\n");

        fprintf(mfp, "REMARK    Section #%4.4ld %ld\n", i, match_list[i][0]);
        fprintf(mfp, "REMARK    %s\n", pairstr);
        fprintf(mfp, "REMARK    %s\n", X3DNA_VER);
        pair2mst(match_list[i][0], match_list[i], AtomName, ResName, ChainID,
                 ResSeq, Miscs, xyz, orien, org, seidx, morien, morg, htm_water, mfp);
        fprintf(mfp, "END\n");

        write_fpmst(morg, morien, rframe);
    }
    close_file(mulbp);

    /* partial-match set of base associations */
    fprintf(fp, "\n===================== partial match =====================\n");
    fprintf(mfp, "\nREMARK    =========================== partial match"
            " ===========================\n");
    for (i = 1; i <= num_partial; i++) {
        m = i + num_match;
        fprintf(fp, "%5ld: #%ld ", m, partial_list[i][0]);
        pairstr[0] = '\0';
        fprintf(rframe, "... %5ld ", m);
        for (j = 1; j <= partial_list[i][0]; j++) {
            jr = partial_list[i][j];
            k = seidx[jr][1];
            base_str(ChainID[k], ResSeq[k], Miscs[k], ResName[k], bseq[jr], 1, b1);
            sprintf(tmp, "[%ld]%s%s", jr, b1, (j == partial_list[i][0]) ? "" : " + ");
            strcat(pairstr, tmp);
            fprintf(rframe, "%c", bseq[jr]);
        }
        fprintf(fp, "%s\n", pairstr);
        fprintf(rframe, " ...\n");

        fprintf(mfp, "REMARK    Section #%4.4ld %ld\n", m, partial_list[i][0]);
        fprintf(mfp, "REMARK    %s\n", pairstr);
        fprintf(mfp, "REMARK    %s\n", X3DNA_VER);
        pair2mst(partial_list[i][0], partial_list[i], AtomName, ResName, ChainID,
                 ResSeq, Miscs, xyz, orien, org, seidx, morien, morg, htm_water, mfp);
        fprintf(mfp, "END\n");

        write_fpmst(morg, morien, rframe);
    }
    close_file(mfp);
    close_file(rframe);

    free_lmatrix(match_list, 1, num_residue, 0, NP);
    free_lmatrix(partial_list, 1, num_residue, 0, NP);
}

void base_compatibility(long num_residue, long **pair_info, long *num_match,
                        long **match_list, long *num_partial, long **partial_list,
                        FILE * fp)
/* check for compatibility of all bases in the network */
{
    long i, j, k, numok, numb;
    long nmatch = 0, npartial = 0;
    long *idx, **full_list;

    idx = lvector(1, NP);
    full_list = lmatrix(1, num_residue, 0, NP);

    /* make a copy of pair info & sort into order */
    for (i = 1; i <= num_residue; i++) {
        if (pair_info[i][NP] > 1) {
            numb = pair_info[i][NP] + 1;
            full_list[i][1] = i;
            for (j = 1; j <= pair_info[i][NP]; j++)
                full_list[i][j + 1] = pair_info[i][j];
            lsort(numb, full_list[i], idx);  /* sort into order */
            full_list[i][0] = numb;
        }
    }
    print_list(num_residue, full_list, "initial lists", fp);

    /* full match list */
    for (i = 1; i <= num_residue; i++) {
        numb = full_list[i][0];
        if (numb <= 0)
            continue;
        numok = 0;
        for (j = 1; j <= numb; j++) {
            k = full_list[i][j];
            if (isequal_list(full_list[i], full_list[k]))
                numok++;
        }
        if (numok == numb) {  /* all match */
            for (j = 1; j <= numb; j++) {  /* mark the rest out */
                k = full_list[i][j];
                full_list[k][0] = (k == i) ? -numb : 0;
            }
        }
    }
    for (i = 1; i <= num_residue; i++)  /* make a clean copy */
        if (full_list[i][0] < 0) {
            nmatch++;
            match_list[nmatch][0] = -full_list[i][0];
            for (j = 1; j <= match_list[nmatch][0]; j++)
                match_list[nmatch][j] = full_list[i][j];
        }
    print_list(num_residue, match_list, "perfect match", fp);

    /* partial match list */
    for (i = 1; i <= num_residue; i++)
        if (full_list[i][0] > 0) {  /* only partial base list */
            for (j = 1; j <= npartial; j++)
                if (isequal_list(full_list[i], partial_list[j]))
                    break;
            if (j > npartial) {
                partial_list[++npartial][0] = full_list[i][0];
                for (j = 1; j <= full_list[i][0]; j++)
                    partial_list[npartial][j] = full_list[i][j];
            }
        }
    print_list(num_residue, partial_list, "partial match", fp);

    *num_match = nmatch;
    *num_partial = npartial;

    free_lvector(idx, 1, NP);
    free_lmatrix(full_list, 1, num_residue, 0, NP);
}

void print_list(long num_residue, long **mylist, char *my_str, FILE * fp)
/* print out selected list for checking: temporary */
{
    long i, j, numb;

    fprintf(fp, "====================== %s ======================\n", my_str);
    for (i = 1; i <= num_residue; i++) {
        numb = mylist[i][0];
        if (!numb)
            continue;
        fprintf(fp, "%4ld[%+2ld]: ", i, numb);
        for (j = 1; j <= labs(numb); j++)
            fprintf(fp, "%5ld", mylist[i][j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

long isequal_list(long *listA, long *listB)
/* check if two lists are identical: same # & residue order */
{
    long i, is_equal = 0;

    if (listA[0] == listB[0]) {
        for (i = 1; i <= listA[0]; i++)
            if (listA[i] != listB[i])
                break;
        if (i > listA[0])
            is_equal = 1;
    }
    return is_equal;
}

void bp_network(long num_residue, long *RY, long **seidx, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, char **Miscs, long *idx,
                double **xyz, long **ring_atom, char *bseq, long **pair_info,
                double **NC1xyz, double **orien, double **org, double *BPRS,
                long hetatm, long **htm_water, FILE * fp)
/* get the base-pair networking system: triple, etc */
{
    char b1[BUF512];
    long i, ir, j, inum_base;
    long tnum_base = 0, max_ple = -1;
    long *ivec;

    for (i = 1; i <= num_residue; i++)
        if (RY[i] >= 0)
            tnum_base++;  /* total number of bases */
    ivec = lvector(1, tnum_base);

    fprintf(fp, "\nDetailed pairing information for each base\n");
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0)
            continue;

        ir = seidx[i][1];
        base_str(ChainID[ir], ResSeq[ir], Miscs[ir], ResName[ir], bseq[i], 1, b1);

        /* list of direct pairing */
        fprintf(fp, "%5ld %s: [%2ld]", i, b1, pair_info[i][NP]);
        for (j = 1; j <= pair_info[i][NP]; j++)
            fprintf(fp, "%5ld", pair_info[i][j]);
        fprintf(fp, "\n");

        inum_base = allbase_cncts(i, tnum_base, ivec, pair_info, fp);
        bases_elimination(i, inum_base, ivec, bseq, seidx, ring_atom, xyz, NC1xyz,
                          orien, org, AtomName, BPRS, b1, idx, &max_ple, pair_info, fp);
    }

    if (max_ple > 1)
        multiplets(max_ple, num_residue, pair_info, AtomName, ResName, ChainID, ResSeq,
                   Miscs, xyz, orien, org, seidx, bseq, hetatm, htm_water, fp);

    free_lvector(ivec, 1, tnum_base);
}

long allbase_cncts(long i, long tnum_base, long *ivec, long **pair_info, FILE * fp)
/* get a list of all base connections */
{
    long ir, j, m, inum_base;

    ivec[1] = i;
    init_lvector(ivec, 2, tnum_base, 0);
    inum_base = 1;

    m = 1;
    while (ivec[m] && m <= tnum_base) {
        ir = ivec[m++];
        for (j = 1; j <= pair_info[ir][NP]; j++)
            if (!lval_in_set(pair_info[ir][j], 1, inum_base, ivec))
                ivec[++inum_base] = pair_info[ir][j];
    }

    /* list of networked pairing */
    fprintf(fp, "                      [%2ld]", inum_base - 1);
    for (j = 2; j <= inum_base; j++)
        fprintf(fp, "%5ld", ivec[j]);
    fprintf(fp, "\n");

    return inum_base;
}

void bases_elimination(long i, long inum_base, long *ivec, char *bseq, long **seidx,
                       long **ring_atom, double **xyz, double **NC1xyz, double **orien,
                       double **org, char **AtomName, double *BPRS, char *b1, long *idx,
                       long *max_ple, long **pair_info, FILE * fp)
/* eliminate incompatible bases: keep only the ones that have {dv, angle,
   and dNN} in range. check for isolated bases */
{
    double rtn_val[RTNNUM];
    long j, k, m, bpid, num_kept = 0, num_final = 0;
    long *idx1, *idx2;

    idx1 = lvector(1, inum_base);
    idx2 = lvector(1, inum_base);

    for (j = 1; j <= inum_base - 1; j++) {
        if (ivec[j] < 0)
            break;
        m = 0;
        for (k = j + 1; k <= inum_base; k++) {
            if (ivec[k] < 0)  /* already ruled out */
                break;
            m++;
            check_pair(ivec[j], ivec[k], bseq, seidx, xyz, NC1xyz, orien, org, idx,
                       AtomName, BPRS, rtn_val, &bpid, ring_atom, 1);
            if (!bpid) {  /* ivec[k] not in a network with ivec[j] */
                idx1[m] = lround(MFACTOR * 12.0);  /* ruled out by making it big */
                ivec[k] = -ivec[k];
            } else
                idx1[m] = lround(MFACTOR * rtn_val[2]);  /* vertical distance */
        }
        if (m > 1) {
            lsort(m, idx1, idx2);  /* sort according to dv */
            for (k = 1; k <= m; k++)
                idx1[k] = ivec[j + idx2[k]];
            for (k = 1; k <= m; k++)
                ivec[k + j] = idx1[k];  /* keep the matched bases */
        }
    }

    for (j = 1; j <= inum_base; j++) {
        if (ivec[j] < 0)
            break;
        num_kept++;
        idx1[num_kept] = 0;  /* clear it up */
    }
    fprintf(fp, "                      [%2ld]", num_kept - 1);
    for (j = 2; j <= num_kept; j++)
        fprintf(fp, "%5ld", ivec[j]);
    fprintf(fp, "\n");

    /* check for connections: delete isolated bases or fragments not
       containing base residue i (i.e, ivec[1]) */
    idx1[1] = 1;  /* start with base i */
    while (1) {
        for (j = 1; j <= num_kept; j++)  /* find the residue to start with */
            if (idx1[j] > 0)
                break;
        if (j > num_kept)
            break;  /* all done! */
        idx1[j] = -1;  /* already checked */
        for (k = 1; k <= num_kept; k++) {
            if (idx1[k])
                continue;
            check_pair(ivec[j], ivec[k], bseq, seidx, xyz, NC1xyz, orien, org, idx,
                       AtomName, BPRS, rtn_val, &bpid, ring_atom, 0);
            if (bpid)
                idx1[k] = 1;
        }
    }

    for (j = 1; j <= num_kept; j++)
        if (idx1[j])
            num_final++;
    fprintf(fp, "                     %s[%2ld]",
            (num_final != num_kept) ? "*" : " ", num_final - 1);
    for (j = 2; j <= num_kept; j++)
        if (idx1[j])
            fprintf(fp, "%5ld", ivec[j]);
    fprintf(fp, "\n");

    /* re-set pair_info for later processing */
    k = 0;
    for (j = 2; j <= num_kept; j++) {
        if (!idx1[j])
            continue;  /* skip isolated base */
        if (++k >= NP) {
            fprintf(stderr, "residue %s has over %ld pairs\n", b1, NP - 1);
            --k;
            break;
        } else
            pair_info[i][k] = ivec[j];
    }
    pair_info[i][NP] = (idx1[1]) ? k++ : 0;  /* could < direct pair_info[i][NP] */

    if (k > 1 && k > *max_ple)
        *max_ple = k;

    free_lvector(idx1, 1, inum_base);
    free_lvector(idx2, 1, inum_base);
}

void all_pairs(long num_residue, long *RY, double **NC1xyz, double **orien,
               double **org, double *BPRS, long **seidx, double **xyz, long *idx,
               long **ring_atom, char **AtomName, char **ResName, char *ChainID,
               long *ResSeq, char **Miscs, char *bseq, long hetatm, long **htm_water,
               char *pdbfile, char *outfile)
/* find all possible base-pairs */
{
    char wc[4], HB_ATOM[BUF512];
    char b1[BUF512], b2[BUF512], idmsg[BUF512];
    char *doc_str[] = {
        "Six-line information for each base-pair as follows:\n",
        "   #1: overall serial number, local serial number, paired residue numbers,\n",
        "       detailed pairing residue information.\n",
        "   #2: one-letter base-pair followed by six base-pair parameters (shear,\n",
        "       stretch, stagger, buckle, propeller, opening). The parameters are\n",
        "       with respect to the Watson-Crick base reference frame.  There are\n",
        "       two types of base-pair orientation: A-B means the two bases have\n",
        "       opposite orientations as in Watson-Crick base-pair; A+B means the\n",
        "       two bases have the same local orientations as in Hoogsteen base-\n",
        "       pair. All possible base pairing patterns can then be classified\n",
        "       based on the six parameters, among which shear, stretch and opening\n",
        "       are most discriminative.\n",
        "   #3: H-bonding information (atom pair followed by their distance).\n",
        "   #4: Overall classification of the base-pair (anti-parallel vs parallel\n",
        "       based on relative z-axis of the two bases, cis vs trans based on\n",
        "       x-axis and C1-RN9/YN1 directions).\n",
        "   #5: Relative directions of the three axes and their numerical values.\n",
        "       The last 3 numbers are the angles between the glycosidic bonds, and\n",
        "       the two chi torsion angles.\n",
        "   #6: The actual parameters used to locate the base-pair in question.\n\n"
    };
    double morg[4], morien[10];
    double rtn_val[RTNNUM], HB_UPPER[2], *chi;
    long bpid, i, inum, ir, j, jr, num_bp = 0;
    long inum_base = 2, ivec[3];  /* for middle frame orientation */
    long **pair_info;
    FILE *fp, *mfp, *rtmp, *rframe;

    hb_crt_alt(HB_UPPER, HB_ATOM, b1);
    pair_info = lmatrix(1, num_residue, 1, NP);  /* detailed base-pair network */
    mfp = open_file(ALLP_FILE, "w");

    fp = open_file(outfile, "w");
    fprintf(fp, "PDB data file name: %s\n", pdbfile);
    print_bp_crit(BPRS, HB_ATOM, fp);

    ir = sizeof doc_str / sizeof doc_str[0];
    for (i = 0; i < ir; i++)
        fprintf(fp, "%s", doc_str[i]);

    chi = dvector(1, num_residue);
    get_chi_angle(num_residue, RY, bseq, seidx, xyz, AtomName, ResName, ChainID,
                  ResSeq, Miscs, chi, NULL);

    rtmp = open_file(TMP_FILE, "w");
    for (i = 1; i < num_residue; i++) {
        if (RY[i] < 0)
            continue;
        inum = 0;
        for (j = i + 1; j <= num_residue; j++) {
            if (RY[j] < 0)
                continue;
            check_pair(i, j, bseq, seidx, xyz, NC1xyz, orien, org, idx, AtomName,
                       BPRS, rtn_val, &bpid, ring_atom, 0);
            if (bpid) {
                num_bp++;
                inum++;
                sprintf(wc, "%c%c%c", (bpid == 2) ? '-' : '*', (bpid > 0) ? '-' : '*',
                        (rtn_val[35] < 0.0) ? '-' : '+');
                ir = seidx[i][1];
                jr = seidx[j][1];
                base_str(ChainID[ir], ResSeq[ir], Miscs[ir], ResName[ir], bseq[i], 1, b1);
                base_str(ChainID[jr], ResSeq[jr], Miscs[jr], ResName[jr], bseq[j], 2, b2);
                fprintf(fp, "%5ld%5ld%5ld%5ld %s-%s-%s\n", num_bp, inum, i, j, b1, wc,
                        b2);
                print_pairinfo(i, j, bseq[i], bseq[j], rtn_val, chi, HB_ATOM, HB_UPPER,
                               seidx, AtomName, xyz, bseq, 1, fp);

                /* multiple structure file containing all base-pairs */
                sprintf(idmsg, "%s-%s-%s", b1, wc, b2);
                fprintf(mfp, "REMARK    Section #%4.4ld %s\n", num_bp, idmsg);
                fprintf(mfp, "REMARK    %s\n", X3DNA_VER);
                ivec[1] = i;
                ivec[2] = j;
                pair2mst(inum_base, ivec, AtomName, ResName, ChainID, ResSeq, Miscs,
                         xyz, orien, org, seidx, morien, morg, htm_water, mfp);
                fprintf(mfp, "END\n");

                fprintf(rtmp, "... %5ld %c%c%c ...\n", num_bp, bseq[i], wc[2], bseq[j]);
                write_fpmst(morg, morien, rtmp);

                if (++pair_info[i][NP] >= NP) {
                    fprintf(stderr, "residue %s has over %ld pairs\n", b1, NP - 1);
                    --pair_info[i][NP];
                    break;
                } else
                    pair_info[i][pair_info[i][NP]] = j;
                if (++pair_info[j][NP] >= NP) {
                    fprintf(stderr, "residue %s has over %ld pairs\n", b2, NP - 1);
                    --pair_info[j][NP];
                    break;
                } else
                    pair_info[j][pair_info[j][NP]] = i;
            }
        }
    }
    close_file(rtmp);

    rtmp = open_file(TMP_FILE, "r");
    rframe = open_file(REF_FILE, "w");
    fprintf(rframe, "%5ld base-pairs\n", num_bp);  /* get total # of base-pairs */
    while (fgets(b1, sizeof b1, rtmp) != NULL)
        fputs(b1, rframe);
    close_file(rtmp);
    close_file(rframe);

    bp_network(num_residue, RY, seidx, AtomName, ResName, ChainID, ResSeq, Miscs, idx,
               xyz, ring_atom, bseq, pair_info, NC1xyz, orien, org, BPRS, hetatm,
               htm_water, fp);

    free_lmatrix(pair_info, 1, num_residue, 1, NP);
    free_dvector(chi, 1, num_residue);
    close_file(fp);
    close_file(mfp);
}

void best_pair(long i, long num_residue, long *RY, long **seidx, double **xyz, long *idx,
               double **NC1xyz, long *matched_idx, double **orien, double **org,
               long **ring_atom, char **AtomName, char *bseq, double *BPRS,
               long *pair_stat)
/* find the best-paired residue id#
 * pair_stat[PSTNUM]:
 *     j, bpid, d, dv, angle, dNN, dsum, bp-org,  x1,   y1,   z1,   x2,   y2,   z2
 *col# 1   2    3   4    5     6     7     8-10  11-13 14-16 17-19 20-22 23-25 26-28
 */
{
    double ddmin = XBIG, rtn_val[RTNNUM];
    long bpid, j, k, nout = PSTNUM - 1;

    init_lvector(pair_stat, 1, nout, 0);

    for (j = 1; j <= num_residue; j++) {
        if (j == i || RY[j] < 0 || matched_idx[j])
            continue;
        check_pair(i, j, bseq, seidx, xyz, NC1xyz, orien, org, idx, AtomName, BPRS,
                   rtn_val, &bpid, ring_atom, 0);
        if (bpid && rtn_val[5] < ddmin) {
            ddmin = rtn_val[5];
            pair_stat[1] = j;
            pair_stat[2] = bpid;
            for (k = 1; k <= nout - 2; k++)
                pair_stat[2 + k] = lround(MFACTOR * rtn_val[k]);
        }
    }
}

long bp_coplanar(long i, double d, double d2, double *txyz, double *txyz2, long n,
                 long *ddidx, double **bp_xyz)
/* check if 2 bps co-planar based on simple geometry criterion */
{
    double *dp = NULL;  /* use pointer for efficiency */
    long j = 0, co_planar = 0;

    if (fabs(d) < OLCRT) {
        j = ddidx[1];
        dp = txyz;
    } else if (fabs(d2) < OLCRT) {
        j = ddidx[n];
        dp = txyz2;
    }
    if (j && fabs(dot(&bp_xyz[i][9], dp)) < OLCRT &&
        fabs(dot(&bp_xyz[i][18], dp)) < OLCRT &&
        fabs(dot(&bp_xyz[j][9], dp)) < OLCRT && fabs(dot(&bp_xyz[j][18], dp)) < OLCRT)
        co_planar = 1;
    return co_planar;
}

void bp_context(long num_bp, double HELIX_CHG, double **bp_xyz, long **bp_order,
                long **end_list, long *num_ends, FILE * tfp)
/* find base-pair neighbors using simple geometric criterion for re_ordering */
{
    double d = EMPTY_NUMBER, d2, d3, ddmin[9], txyz[4], txyz2[4], txyz3[4], zave[4];
    long i, j, k, m, n, overlap = 0, quadruple = 0, cnum = 8, ddidx[9];

    fprintf(tfp, "\nBase-pair context information\n");

    for (i = 1; i <= num_bp; i++) {
        init_dvector(ddmin, 1, cnum, XBIG);
        init_lvector(ddidx, 1, cnum, 0);
        d = dot(&bp_xyz[i][9], &bp_xyz[i][18]);  /* between base normals */
        (d <= 0.0) ? ddxyz(bp_xyz[i] + 18, bp_xyz[i] + 9, zave) :
            sumxyz(bp_xyz[i] + 18, bp_xyz[i] + 9, zave);
        vec_norm(zave);
        for (j = 1; j <= num_bp; j++) {
            if (j == i)
                continue;
            ddxyz(bp_xyz[j], bp_xyz[i], txyz);
            d = veclen(txyz);
            for (k = 1; k <= cnum; k++)
                if (d < ddmin[k]) {
                    for (m = cnum; m > k; m--)
                        if (ddidx[n = m - 1]) {
                            ddmin[m] = ddmin[n];
                            ddidx[m] = ddidx[n];
                        }
                    ddmin[k] = d;
                    ddidx[k] = j;
                    break;
                }
        }
        if (ddidx[1] && ddidx[2]) {  /* at least 2 nearest neighbors */
            if (ddmin[1] > HELIX_CHG) {  /* isolated bp */
                end_list[++*num_ends][1] = i;  /* [i 0 0] */
                n = 2;
            } else {
                if (!overlap && ddmin[1] < OLCRT)
                    overlap = 1;
                ddxyz(bp_xyz[ddidx[1]], bp_xyz[i], txyz);  /* i's nearest neighbor */
                d = dot(zave, txyz);

                /* check for swapped bps in deformed pdr026/ur0009 */
                if (ddidx[3] && ddmin[2] <= HELIX_CHG && ddmin[3] <= HELIX_CHG) {
                    ddxyz(bp_xyz[ddidx[2]], bp_xyz[i], txyz2);
                    ddxyz(bp_xyz[ddidx[3]], bp_xyz[i], txyz3);
                    d2 = dot(zave, txyz2);
                    d3 = dot(zave, txyz3);
                    if (d * d2 < 0.0 && d * d3 < 0.0 && fabs(d2) > fabs(d3)) {
                        lval_swap(&ddidx[2], &ddidx[3]);
                        dval_swap(&ddmin[2], &ddmin[3]);
                        fprintf(stderr, "[swapping 2nd & 3rd] %4ld %8.2f %8.2f "
                                "%8.2f %8.2f\n", i, ddmin[2], ddmin[3], d2, d3);
                    }
                }
                n = 0;
                for (j = 2; j <= cnum && ddidx[j]; j++) {
                    if (ddmin[j] > HELIX_CHG)
                        break;
                    ddxyz(bp_xyz[ddidx[j]], bp_xyz[i], txyz2);
                    d2 = dot(zave, txyz2);
                    if (d * d2 < 0.0) {  /* vertical direction */
                        n = j;
                        bp_order[i][1] = -1;  /* middle base-pair */
                        bp_order[i][2] = ddidx[1];
                        bp_order[i][3] = ddidx[j];
                        break;
                    }
                }
                if (!n) {  /* terminal bp */
                    n = 2;
                    end_list[++*num_ends][1] = i;
                    end_list[*num_ends][2] = ddidx[1];
                    bp_order[i][2] = ddidx[1];
                    ddxyz(bp_xyz[ddidx[1]], bp_xyz[ddidx[n]], txyz2);
                    d2 = dot(zave, txyz2);  /* check for normal terminal */
                    if (d * d2 < 0.0 && veclen(txyz2) <= HELIX_CHG) {
                        end_list[*num_ends][3] = ddidx[n];
                        bp_order[i][3] = ddidx[n];
                    }
                }
            }
            fprintf(tfp, "%4ld: %4ld %4ld %4ld %8.2f %8.2f%c", i, bp_order[i][1],
                    bp_order[i][2], bp_order[i][3], ddmin[1], ddmin[n],
                    (ddmin[n] > HELIX_CHG) ? '*' : ' ');
            if (!bp_order[i][2])
                fprintf(tfp, "  isolated base-pairs\n");
            else {
                (!bp_order[i][3]) ? fprintf(tfp, " (%4ld)", ddidx[n]) :
                    fprintf(tfp, "       ");
                ddxyz(bp_xyz[ddidx[n]], bp_xyz[i], txyz2);
                d2 = dot(zave, txyz2);
                if (!overlap && !quadruple)
                    quadruple = bp_coplanar(i, d, d2, txyz, txyz2, n, ddidx, bp_xyz);
                fprintf(tfp, " ==> %8.2f %8.2f %c", d, d2, (d * d2 > 0.0) ? '*' : ' ');
                d = magang(txyz, txyz2);  /* magang normalize txyz & txyz2 */
                fprintf(tfp, " (%8.2f%c)\n", d, (d <= 90.0) ? '*' : ' ');
            }
        }
    }

    if (!*num_ends) {  /* num_bp == 1 || 2 */
        end_list[++*num_ends][1] = 1;
        if (num_bp == 2) {
            if (d <= HELIX_CHG) {
                end_list[*num_ends][2] = 2;  /* 1 2 0 && 2 1 0 */
                end_list[++*num_ends][1] = 2;
                end_list[*num_ends][2] = 1;
            } else
                end_list[++*num_ends][1] = 2;  /* 1 0 0 && 2 0 0 */
        }
    }
    fprintf(tfp, "\nEnd base-pair list\n");
    for (i = 1; i <= *num_ends; i++)
        fprintf(tfp, "%4ld: %4ld %4ld %4ld\n", i, end_list[i][1], end_list[i][2],
                end_list[i][3]);
    if (overlap)
        fprintf(stderr, "***Warning: structure with overlapped base-pairs***\n");
    else if (quadruple)
        fprintf(stderr, "***Warning: structure with 2 neighbor bps co-planar ***\n");
}

void locate_helix(long num_bp, long **helix_idx, long num_ends, long *num_helix,
                  long **end_list, long **bp_order, long *bp_idx, long *helix_marker)
/* locate all possible helical regions, including isolated base-pairs */
{
    long i, ip = 0, j, k, k0, k2, k3, m;
    long *matched_idx;

    helix_idx[*num_helix][1] = 1;

    matched_idx = lvector(1, num_bp);  /* indicator for used bps */

    for (i = 1; i <= num_ends && ip < num_bp; i++) {
        k = 0;
        k0 = 0;
        for (j = 1; j <= 3; j++)
            if (end_list[i][j]) {
                k += matched_idx[end_list[i][j]];
                k0++;
            }
        if (k == k0)
            continue;  /* end point of a processed helix */
        for (j = 1; j <= 3 && ip < num_bp; j++) {
            k = end_list[i][j];
            if (k && !matched_idx[k]) {
                bp_idx[++ip] = k;
                matched_idx[k] = 1;
            }
        }
        for (j = 1; j <= num_bp; j++) {
            k = bp_idx[ip];
            k2 = bp_order[k][2];
            k3 = bp_order[k][3];
            if (!bp_order[k][1]) {  /* end-point */
                if (k2 && !matched_idx[k2] && !k3) {  /* e.g. 0 10 0 */
                    bp_idx[++ip] = k2;
                    matched_idx[k2] = 1;
                }
                break;  /* normal case */
            }
            m = matched_idx[k2] + matched_idx[k3];
            if (m == 2 || m == 0)
                break;  /* chain terminates */
            if (k2 == bp_idx[ip - 1]) {
                bp_idx[++ip] = k3;
                matched_idx[k3] = 1;
            } else if (k3 == bp_idx[ip - 1]) {
                bp_idx[++ip] = k2;
                matched_idx[k2] = 1;
            } else
                break;  /* no direct connection */
        }
        helix_idx[*num_helix][2] = ip;
        helix_marker[ip] = 1;  /* helix_marker & helix_idx are parallel */
        if (ip < num_bp)
            helix_idx[++*num_helix][1] = ip + 1;
    }

    if (ip < num_bp) {  /* all un-classified bps */
        fprintf(stderr, "[%ld %ld]: complicated structure, left over"
                " base-pairs put into the last region [%ld]\n", ip, num_bp, *num_helix);
        helix_idx[*num_helix][2] = num_bp;
        helix_marker[num_bp] = 1;
        for (j = 1; j <= num_bp; j++)
            if (!matched_idx[j])
                bp_idx[++ip] = j;
    }
    free_lvector(matched_idx, 1, num_bp);
}

void first_step(long i, long **helix_idx, long *bp_idx, long *swapped, long **base_pairs,
                double **o3_p)
/* make strand I of the first step in 5'--->3' direction */
{
    long i1, i2, j, j1, j2, k, m, n;

    if (helix_idx[i][3] == 1)
        return;
    j = helix_idx[i][1];
    m = bp_idx[j];
    n = bp_idx[j + 1];
    get_ij(m, swapped, base_pairs, &i1, &j1);
    get_ij(n, swapped, base_pairs, &i2, &j2);
    k = is_linked(i1, i2, o3_p);
    if (k == -1)
        swapped[m] = !swapped[m];
    else if (!k) {  /* no linkage: reverse the strand direction */
        lreverse(helix_idx[i][1], helix_idx[i][3], bp_idx);
        j = helix_idx[i][1];
        m = bp_idx[j];
        n = bp_idx[j + 1];
        get_ij(m, swapped, base_pairs, &i1, &j1);
        get_ij(n, swapped, base_pairs, &i2, &j2);
        k = is_linked(i1, i2, o3_p);
        if (k == -1)
            swapped[m] = !swapped[m];
        else if (!k)  /* still no connection */
            lreverse(helix_idx[i][1], helix_idx[i][3], bp_idx);  /* reverse back */
    }
}

long chain1dir(long m, long n, long *swapped, long **base_pairs, double **o3_p)
/* make sure strand I is in the 5'--->3' direction */
{
    long i1, i2, j1, j2, k, irev0 = 0, irev1 = 1;

    get_ij(m, swapped, base_pairs, &i1, &j1);
    get_ij(n, swapped, base_pairs, &i2, &j2);
    k = is_linked(i1, i2, o3_p);
    return (k == -1) ? irev1 : irev0;
}

double distance_ab(double **o3_p, long ia, long ib, long ipa, long ipb)
/* calculate distance between ia & ib: ipa & ipb mark if ia/ib exist */
{
    return (o3_p[ia][ipa] > 0.0 && o3_p[ib][ipb] > 0.0) ?
        p1p2_dist(o3_p[ia] + ipa - 4, o3_p[ib] + ipb - 4) : -1.0;
}

void get_ij(long m, long *swapped, long **base_pairs, long *n1, long *n2)
/* get the two base indices [*n1 and *n2] of base-pair m */
{
    if (swapped[m]) {
        *n1 = base_pairs[m][2];
        *n2 = base_pairs[m][1];
    } else {
        *n1 = base_pairs[m][1];
        *n2 = base_pairs[m][2];
    }
}

void get_bidx(long m, long *swapped, long *idx1, long *idx2)
/* Get base normal position index in bp_xyz */
{
    if (swapped[m]) {
        *idx1 = 18;
        *idx2 = 9;
    } else {
        *idx1 = 9;
        *idx2 = 18;
    }
}

double wcbp_zdir(long m, long n, long idxm1, long idxm2, long idxn1, long idxn2,
                 double **bp_xyz)
/* get the dot product between base-pair m & n normals in WC geometry */
{
    double dm[4], dn[4];

    ddxyz(bp_xyz[m] + idxm2, bp_xyz[m] + idxm1, dm);  /* opposite direction */
    ddxyz(bp_xyz[n] + idxn2, bp_xyz[n] + idxn1, dn);
    vec_norm(dm);
    vec_norm(dn);

    return dot(dm, dn);
}

double wcbp_xang(long m, long n, double **bp_xyz)
/* get the dot product between base-pair m & n x-axes in WC geometry */
{
    double dm[4], dn[4];

    sumxyz(bp_xyz[m] + 3, bp_xyz[m] + 12, dm);  /* same direction */
    sumxyz(bp_xyz[n] + 3, bp_xyz[n] + 12, dn);

    return magang(dm, dn);
}

long wc_bporien(long m, long n, long *swapped, long **base_pairs, double **bp_xyz,
                double **o3_p)
/* get the relative base-pair orientations of WC-geometry step */
{
    long i1, i2, j1, j2, idxm1, idxm2, idxn1, idxn2, irev = 0;

    if (base_pairs[m][3] > 0 && base_pairs[n][3] > 0) {  /* WC geometry */
        get_ij(m, swapped, base_pairs, &i1, &j1);
        get_ij(n, swapped, base_pairs, &i2, &j2);
        if (wcbp_xang(m, n, bp_xyz) > END_STACK_XANG ||  /* bdl070: 105 degrees */
            is_linked(i1, i2, o3_p) || is_linked(j1, j2, o3_p))
            return irev;  /* apposite bp orientation */
        get_bidx(m, swapped, &idxm1, &idxm2);
        get_bidx(n, swapped, &idxn1, &idxn2);
        if (wcbp_zdir(m, n, idxm1, idxm2, idxn1, idxn2, bp_xyz) < 0.0 &&
            wcbp_zdir(m, n, idxm1, idxm2, idxn2, idxn1, bp_xyz) > 0.0)
            irev = 1;
    }
    return irev;
}

long check_o3dist(long m, long n, long *swapped, long **base_pairs, double **o3_p)
/* Use O3* distance criterion to decide if the top bp should be flipped:
 *    I       II
 *    i2 ---- j2
 *    |       |
 *    i1 ---- j1 */
{
    double di1_i2, di1_j2, dj1_i2, dj1_j2;
    long i1, i2, j1, j2, irev = 0;

    get_ij(m, swapped, base_pairs, &i1, &j1);
    get_ij(n, swapped, base_pairs, &i2, &j2);
    di1_i2 = distance_ab(o3_p, i1, i2, 4, 4);
    di1_j2 = distance_ab(o3_p, i1, j2, 4, 4);
    dj1_i2 = distance_ab(o3_p, j1, i2, 4, 4);
    dj1_j2 = distance_ab(o3_p, j1, j2, 4, 4);
    if ((di1_i2 > 0.0 && di1_j2 > 0.0 && di1_i2 > di1_j2) &&
        (dj1_i2 > 0.0 && dj1_j2 > 0.0 && dj1_j2 > dj1_i2))
        irev = 1;
    return irev;
}

long is_linked(long i, long j, double **o3_p)
/* Use O3' and P distance to check if residues i & j are directly connected */
{
    double d;

    d = distance_ab(o3_p, i, j, 4, 8);  /* O3'[i]-P[j] */
    if (dval_in_range(d, 0.0, O3P_UPPER))
        return 1L;  /* 5'--->3' linkage */
    d = distance_ab(o3_p, j, i, 4, 8);  /* O3'[j]-P[i] */
    if (dval_in_range(d, 0.0, O3P_UPPER))
        return -1L;  /* 3'--->5' linkage */
    return 0L;  /* no direct linkage */
}

long check_schain(long m, long n, long *swapped, long **base_pairs, double **o3_p)
/* Use continuous single chain (if any) to control if to flip the top bp */
{
    long i1, i2, j1, j2, irev = 0;

    get_ij(m, swapped, base_pairs, &i1, &j1);
    get_ij(n, swapped, base_pairs, &i2, &j2);
    if (!is_linked(i1, i2, o3_p) && !is_linked(j1, j2, o3_p) &&
        (is_linked(i1, j2, o3_p) || is_linked(j1, i2, o3_p)))
        irev = 1;
    return irev;
}

long check_others(long m, long n, long *swapped, long **base_pairs, double **o3_p,
                  double **bp_xyz)
/* Check for situations where wc_bporien, check_o3dist & check_schain do not apply */
{
    double d[5], a1[4], a2[4], r1[4], r2[4];
    long i, j, i1, i2, j1, j2, idxm1, idxm2, idxn1, idxn2;
    long irev0 = 0, irev1 = 1;

    /* check if any two residues are directly connected */
    get_ij(m, swapped, base_pairs, &i1, &j1);
    get_ij(n, swapped, base_pairs, &i2, &j2);
    if (is_linked(i1, i2, o3_p) || is_linked(j1, j2, o3_p) ||
        is_linked(i1, j2, o3_p) || is_linked(j1, i2, o3_p))
        return irev0;

    get_bidx(m, swapped, &idxm1, &idxm2);
    get_bidx(n, swapped, &idxn1, &idxn2);

    /* original orientation */
    for (i = 1; i <= 3; i++) {
        j = (3 - i) * 3;
        a1[i] = dot(&bp_xyz[m][idxm1 - j], &bp_xyz[n][idxn1 - j]);
        a2[i] = dot(&bp_xyz[m][idxm2 - j], &bp_xyz[n][idxn2 - j]);
    }
    i1 = (a1[1] > 0.0 && a1[2] > 0.0 && a1[3] > 0.0);
    i2 = (a2[1] > 0.0 && a2[2] > 0.0 && a2[3] > 0.0);
    if (i1 && i2)  /* as is */
        return irev0;

    /* swap the two bases in the top base-pair */
    for (i = 1; i <= 3; i++) {
        j = (3 - i) * 3;
        r1[i] = dot(&bp_xyz[m][idxm1 - j], &bp_xyz[n][idxn2 - j]);
        r2[i] = dot(&bp_xyz[m][idxm2 - j], &bp_xyz[n][idxn1 - j]);
    }
    j1 = (r1[1] > 0.0 && r1[2] > 0.0 && r1[3] > 0.0);
    j2 = (r2[1] > 0.0 && r2[2] > 0.0 && r2[3] > 0.0);
    if (!i1 && !i2) {
        if (j1 || j2)  /* reverse */
            return irev1;
        else if (!j1 && !j2)  /* keep as is */
            return irev0;
    }
    d[1] = dot2ang(a1[1]) + dot2ang(a1[2]) + dot2ang(a1[3]);  /* i1 */
    d[2] = dot2ang(a2[1]) + dot2ang(a2[2]) + dot2ang(a2[3]);  /* i2 */
    d[3] = dot2ang(r1[1]) + dot2ang(r1[2]) + dot2ang(r1[3]);  /* j1 */
    d[4] = dot2ang(r2[1]) + dot2ang(r2[2]) + dot2ang(r2[3]);  /* j2 */

    if (i1 && j1)
        return (d[1] > d[3]) ? irev1 : irev0;
    if (i1 && j2)
        return (d[1] > d[4]) ? irev1 : irev0;
    if (i2 && j1)
        return (d[2] > d[3]) ? irev1 : irev0;
    if (i2 && j2)
        return (d[2] > d[4]) ? irev1 : irev0;

    return irev0;
}

void check_direction(long i, long **helix_idx, long *bp_idx, long *swapped,
                     long **base_pairs, double **o3_p, long *direction)
/* check directions of two chains, mark and reverse as necessary a */
{
    long i1, i2, j, j1, j2, k, m, n;

    /* check if strand I is in 5'-->3' direction */
    init_lvector(direction, 1, 6, 0);
    for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
        m = bp_idx[j];
        n = bp_idx[j + 1];
        get_ij(m, swapped, base_pairs, &i1, &j1);
        get_ij(n, swapped, base_pairs, &i2, &j2);
        k = is_linked(i1, i2, o3_p);
        (k == 1) ? ++direction[1] : (k == -1) ? ++direction[2] : ++direction[3];
        k = is_linked(j1, j2, o3_p);
        (k == 1) ? ++direction[4] : (k == -1) ? ++direction[5] : ++direction[6];
    }
    if ((direction[1] && direction[2]) || (direction[4] && direction[5])) {
        helix_idx[i][7] = 1;
        return;
    }
    if (direction[1] + direction[2] + direction[4] + direction[5] == 0)
        return;  /* all empty: no direct connection */

    m = bp_idx[helix_idx[i][1]];  /* start bp */
    n = bp_idx[helix_idx[i][2]];  /* end bp */
    get_ij(m, swapped, base_pairs, &i1, &j1);
    get_ij(n, swapped, base_pairs, &i2, &j2);
    if (direction[3] || direction[6])
        helix_idx[i][5] = 1;  /* broken O3'[i] to P[i+1] linkage */
    if (direction[1] && !direction[2]) {  /* 5'--->3' direction */
        if (!direction[4] && direction[5]) {  /* anti-parallel */
            if (i1 > j2) {  /* swap and reverse the two strands */
                for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++)
                    swapped[bp_idx[j]] = !swapped[bp_idx[j]];
                lreverse(helix_idx[i][1], helix_idx[i][3], bp_idx);
            }
        } else if (direction[4] && !direction[5]) {  /* parallel */
            helix_idx[i][6] = 1;
            if (i1 > j1)  /* swapped the two strands */
                for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++)
                    swapped[bp_idx[j]] = !swapped[bp_idx[j]];
        }
    }
}

void check_strand2(long i, long **helix_idx, long *bp_idx, double **bp_xyz,
                   long *swapped, long **base_pairs, double **o3_p, long *direction,
                   FILE * tfp)
/* direction of strand II: either anti-parallel or parallel to strand I */
{
    long i1, i2, j, j1, j2, k, m, n, anti_p, parallel;

    if (!helix_idx[i][7]) {
        if (direction[1] + direction[2] + direction[4] + direction[5] == 0)
            return;
        init_lvector(helix_idx[i], 5, 7, 0);
        for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
            m = bp_idx[j];
            n = bp_idx[j + 1];
            if (wc_bporien(m, n, swapped, base_pairs, bp_xyz, o3_p))
                continue;  /* in good WC-geometry orientation */
            get_ij(m, swapped, base_pairs, &i1, &j1);
            get_ij(n, swapped, base_pairs, &i2, &j2);
            if (!is_linked(i1, i2, o3_p) && !is_linked(j1, j2, o3_p) &&
                ((is_linked(i1, j2, o3_p) == 1) ||
                 (is_linked(i1, j2, o3_p) && is_linked(j1, i2, o3_p)))) {
                swapped[n] = !swapped[n];  /* pr0004, ptr012, trna12 */
                fprintf(tfp, "                  000    [%ld-%ld]\n", m, n);
                fprintf(stderr, "000:    [%ld-%ld]\n", m, n);
            }
        }
    } else {
        helix_idx[i][7] = 0;
        anti_p = (direction[1] > direction[2]) && (direction[4] < direction[5]);
        parallel = (direction[1] > direction[2]) && (direction[4] > direction[5]);
        for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
            m = bp_idx[j];
            n = bp_idx[j + 1];
            get_ij(m, swapped, base_pairs, &i1, &j1);
            get_ij(n, swapped, base_pairs, &i2, &j2);
            k = is_linked(j1, j2, o3_p);
            if (!is_linked(i1, i2, o3_p) &&  /* strand I not linked */
                ((anti_p && k == 1) || (parallel && k == -1))) {
                swapped[n] = !swapped[n];
                fprintf(tfp, "                  2nd %2ld [%ld-%ld]\n", k, m, n);
            }
            /* another round of post-processing */
            get_ij(n, swapped, base_pairs, &i2, &j2);
            if (!is_linked(i1, i2, o3_p) && !is_linked(j1, j2, o3_p)) {
                if ((anti_p && is_linked(j1, i2, o3_p) == 1) ||
                    (parallel && is_linked(i1, j2, o3_p) == -1)) {
                    fprintf(tfp, "                  3rdL   [%ld-%ld]\n", m, n);
                    swapped[m] = !swapped[m];  /* lower bp */
                } else if ((anti_p && is_linked(i1, j2, o3_p) == 1) ||
                           (parallel && is_linked(j1, i2, o3_p) == -1)) {
                    swapped[n] = !swapped[n];  /* upper bp */
                    fprintf(tfp, "                  3rdU    [%ld-%ld]\n", m, n);
                }
            }
        }
    }
    check_direction(i, helix_idx, bp_idx, swapped, base_pairs, o3_p, direction);
}

void check_rise(long i, long **helix_idx, long *bp_idx, long *swapped, long **base_pairs,
                double **bp_xyz, double **o3_p)
/* check WC-geometry step with negative rise */
{
    double d1, d2, rise, dorg[4], mn[4];
    long i1, j1, i2, j2, idxm1, idxm2, idxn1, idxn2, j, k, m, n, num = 0;

    for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
        m = bp_idx[j];
        n = bp_idx[j + 1];
        if (base_pairs[m][3] > 0 && base_pairs[n][3] > 0) {  /* WC geometry */
            get_bidx(m, swapped, &idxm1, &idxm2);
            get_bidx(n, swapped, &idxn1, &idxn2);
            if (wcbp_zdir(m, n, idxm1, idxm2, idxn1, idxn2, bp_xyz) < 0.0 &&
                wcbp_xang(m, n, bp_xyz) > END_STACK_XANG)
                fprintf(stderr, "^^vv opposite bp direction: %ld(%ld)"
                        " %ld(%ld)-%ld(%ld)\n", i, helix_idx[i][3], m, j, n, j + 1);
            ddxyz(bp_xyz[m], bp_xyz[n], dorg);
            for (k = 1; k <= 3; k++)
                mn[k] = bp_xyz[m][k + idxm1] - bp_xyz[m][k + idxm2] +  /* bp m */
                    bp_xyz[n][k + idxn1] - bp_xyz[n][k + idxn2];  /* bp n */
            vec_norm(mn);
            rise = dot(dorg, mn);
            if (rise < 0.0) {
                get_ij(m, swapped, base_pairs, &i1, &j1);
                get_ij(n, swapped, base_pairs, &i2, &j2);
                d1 = distance_ab(o3_p, i2, i1, 4, 8);
                d2 = distance_ab(o3_p, j1, j2, 4, 8);
                if (dval_in_range(d1, 0.0, O3P_UPPER)
                    && dval_in_range(d2, 0.0, O3P_UPPER)) {
                    num++;
                    fprintf(stderr, "===> %ld(%ld) %ld-%ld [%ld-%ld (%ld)]:"
                            " %8.2f%8.2f%8.2f\n", i, helix_idx[i][3], m, n,
                            j, j + 1, helix_idx[i][2], rise, d1, d2);
                }
            }
        }
    }

    if (num)
        fprintf(stderr, "****** Please check base-pair ordering ******\n\n");
}

void five2three(long num_bp, long *num_helix, long **helix_idx, long *bp_idx,
                double **bp_xyz, long **base_pairs, double **o3_p, FILE * tfp)
/* make sure the leading strand is in the 5' to 3' direction
 * helix_idx[][7]: start#, ending#, total#, Z-DNA, break, parallel, weird
 *           col#     1       2       3       4      5       6        7
 */
{
    double do3_p;
    long i, j, k, m, n, rev_wc, rev_o3d, rev_csc, rev_oth, rev_s1;
    long direction[7], *swapped;

    /* check if O3'[i] is wrongly connected to P[i] */
    for (i = 1; i <= *num_helix; i++)
        for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
            k = base_pairs[bp_idx[j]][1];
            do3_p = distance_ab(o3_p, k, k, 4, 8);
            if (dval_in_range(do3_p, 0.0, O3P_UPPER)) {
                fprintf(stderr, "Wrong: O3'[i] connected to P[i]\n");
                return;  /* ignore 5'-->3' re-arrangement */
            }
        }

    fprintf(tfp, "\nHelix region information\n");
    swapped = lvector(1, num_bp);
    for (i = 1; i <= *num_helix; i++) {
        helix_idx[i][3] = helix_idx[i][2] - helix_idx[i][1] + 1;
        print_sep(tfp, '-', 84);
        fprintf(tfp, "Helix #%4.4ld\n", i);

        first_step(i, helix_idx, bp_idx, swapped, base_pairs, o3_p);

        /* The first base-pair is used as the global reference */
        for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
            m = bp_idx[j];
            n = bp_idx[j + 1];
            rev_wc = wc_bporien(m, n, swapped, base_pairs, bp_xyz, o3_p);
            rev_o3d = check_o3dist(m, n, swapped, base_pairs, o3_p);
            rev_csc = check_schain(m, n, swapped, base_pairs, o3_p);
            rev_oth = check_others(m, n, swapped, base_pairs, o3_p, bp_xyz);
            fprintf(tfp, "          %4ld: %2ld %2ld %2ld %2ld", j, rev_wc, rev_o3d,
                    rev_oth, rev_csc);
            if (rev_wc)
                swapped[n] = !swapped[n];
            else if (rev_o3d || rev_csc || rev_oth)  /* rev_o3d: udf027 */
                swapped[n] = !swapped[n];
            /* 5'--->3' direction of strand I as final proof */
            rev_s1 = chain1dir(m, n, swapped, base_pairs, o3_p);
            if (rev_s1)
                swapped[n] = !swapped[n];
            fprintf(tfp, " %2ld [%ld-%ld]\n", rev_s1, m, n);
        }

        /* one more round checking for WC geometry steps: 09-23-2002
         * works for Luger pol_end1.pdb & pol_end2.pdb, and bdl070.pdb */
        fprintf(tfp, "\n              ===> 2nd around checking or WC geometry steps\n");
        for (j = helix_idx[i][1]; j < helix_idx[i][2]; j++) {
            m = bp_idx[j];
            n = bp_idx[j + 1];
            rev_wc = wc_bporien(m, n, swapped, base_pairs, bp_xyz, o3_p);
            if (rev_wc) {
                swapped[m] = !swapped[m];  /* m here */
                fprintf(tfp, "          %4ld: [%ld-%ld]\n", j, m, n);
            }
        }
        fprintf(tfp, "\n");

        check_direction(i, helix_idx, bp_idx, swapped, base_pairs, o3_p, direction);
        check_strand2(i, helix_idx, bp_idx, bp_xyz, swapped, base_pairs, o3_p,
                      direction, tfp);
        check_rise(i, helix_idx, bp_idx, swapped, base_pairs, bp_xyz, o3_p);

        for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
            m = bp_idx[j];
            if (swapped[m]) {
                lval_swap(&base_pairs[m][1], &base_pairs[m][2]);
                for (k = 1; k <= 9; k++) {  /* swap base I and II xyz-axes */
                    dval_swap(&bp_xyz[m][k + 3], &bp_xyz[m][k + 12]);
                    lval_swap(&base_pairs[m][k + 11], &base_pairs[m][k + 20]);
                }
            }
        }

        fprintf(tfp, "\n%4ld [%c%c%c]:", helix_idx[i][3], helix_idx[i][5] ? 'b' : '-',
                helix_idx[i][6] ? 'p' : '-', helix_idx[i][7] ? '?' : '-');
        for (j = 1; j <= 6; j++)
            fprintf(tfp, "%6ld", direction[j]);
        fprintf(tfp, "\n           ");
        k = 0;
        for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
            m = bp_idx[j];
            fprintf(tfp, "%6ld", swapped[m] ? -m : m);
            if (++k % 10 == 0 && k != helix_idx[i][3])
                fprintf(tfp, "\n           ");
        }
        fprintf(tfp, "\n");
    }

    free_lvector(swapped, 1, num_bp);
}

void check_zdna(long *num_helix, long **helix_idx, long *bp_idx, double **bp_xyz,
                long **base_pairs, FILE * tfp)
/* check if a helical region is in left-handed Z-form */
{
    double txyz[4];
    long i, j, m, n, nweird = 0, nrev, mixed_rl = 0;

    fprintf(tfp, "\nZ-DNA helical region if any\n");
    for (i = 1; i <= *num_helix; i++) {
        if (helix_idx[i][5] || helix_idx[i][6] || helix_idx[i][7] || helix_idx[i][3] <= 1) {
            nweird++;
            continue;  /* break/parallel/weird/only one pair */
        }
        nrev = 0;
        for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
            m = bp_idx[j];
            if (helix_idx[i][3] == 1)
                continue;
            if (j < helix_idx[i][2]) {
                n = bp_idx[j + 1];
                ddxyz(bp_xyz[m], bp_xyz[n], txyz);
            }
            /* with base 1 normal & WC-geometry */
            if (dot(txyz, &bp_xyz[m][9]) < 0.0 && base_pairs[m][3] > 0)
                nrev++;
            else
                break;
        }
        if (nrev == helix_idx[i][3]) {
            helix_idx[i][4] = 1;
            mixed_rl++;
            fprintf(tfp, "Helix #%4.4ld (%4ld) is a Z-DNA\n", i, helix_idx[i][3]);
        }
    }

    if (!nweird && mixed_rl && mixed_rl != *num_helix)
        fprintf(stderr, "This structure has right-/left-handed helical regions\n");
}

void wc3(long *pair_k, char *wc)
/* Kth base-pair classification: WC, WC-geometry, Anti-parallel */
{
    long i, j;
    double z1[4], z2[4];

    for (i = 18; i <= 20; i++) {
        j = i - 17;
        z1[j] = pair_k[i] / MFACTOR;
        z2[j] = pair_k[i + 9] / MFACTOR;
    }
    sprintf(wc, "%c%c%c", (pair_k[3] == 2) ? '-' : '*', (pair_k[3] > 0) ? '-' : '*',
            (dot(z1, z2) < 0.0) ? '-' : '+');
}

void re_ordering(long num_bp, long **base_pairs, long *bp_idx, long *helix_marker,
                 long **helix_idx, double *BPRS, long *num_helix, double **o3_p,
                 char *bseq, long **seidx, char **ResName, char *ChainID,
                 long *ResSeq, char **Miscs)
/* re-order base-pairs into separate helical regions */
{
    char b1[BUF512], b2[BUF512], wc[4];
    double **bp_xyz;
    long i, i_order, j, j_order, num_ends = 0;
    long **bp_order, **end_list;
    FILE *tfp;

    tfp = open_file(BPORDER_FILE, "w");
    print_bp_crit(BPRS, NULL, tfp);
    fprintf(tfp, "Base-pair information BEFORE re-ordering\n");
    for (i = 1; i <= num_bp; i++) {
        wc3(base_pairs[i], wc);
        i_order = base_pairs[i][1];
        j_order = base_pairs[i][2];
        j = seidx[i_order][1];
        base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[i_order], 1, b1);
        j = seidx[j_order][1];
        base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[j_order], 2, b2);
        fprintf(tfp, "%5ld: %5ld%5ld %s-%s-%s", i, i_order, j_order, b1, wc, b2);
        for (j = 4; j <= 8; j++)  /* d, dv, angle, dNN */
            fprintf(tfp, "%6.2f", base_pairs[i][j] / MFACTOR);
        for (j = 9; j <= 11; j++)  /* bp origin */
            fprintf(tfp, "%8.2f", base_pairs[i][j] / MFACTOR);
        fprintf(tfp, "\n");
    }

    /* bp origin, base I/II xyz-axes copied from columns 9-29 of base_pairs
     *       bp_orgin,   x1,    y1,   z1,    x2,    y2,    z2
     * col#     1-3     4-6    7-9   10-12  13-15  16-18  19-21 */
    bp_xyz = dmatrix(1, num_bp, 1, 21);
    for (i = 1; i <= num_bp; i++)
        for (j = 1; j <= 21; j++)
            bp_xyz[i][j] = base_pairs[i][j + 8] / MFACTOR;

    bp_order = lmatrix(1, num_bp, 1, 3);
    end_list = lmatrix(1, num_bp, 1, 3);

    bp_context(num_bp, BPRS[6], bp_xyz, bp_order, end_list, &num_ends, tfp);

    locate_helix(num_bp, helix_idx, num_ends, num_helix, end_list, bp_order, bp_idx,
                 helix_marker);

    five2three(num_bp, num_helix, helix_idx, bp_idx, bp_xyz, base_pairs, o3_p, tfp);

    check_zdna(num_helix, helix_idx, bp_idx, bp_xyz, base_pairs, tfp);

    close_file(tfp);
    free_dmatrix(bp_xyz, 1, num_bp, 1, 21);
    free_lmatrix(bp_order, 1, num_bp, 1, 3);
    free_lmatrix(end_list, 1, num_bp, 1, 3);
}

void helix_info(long **helix_idx, long idx, FILE * fp)
/* print out a summary of helix information: Z-DNA, break, parallel, weird? */
{
    fprintf(fp, "%s%s%s%s\n", (helix_idx[idx][4]) ? "  ***Z-DNA***" : "",
            (helix_idx[idx][5]) ? "  ***broken O3'[i] to P[i+1] linkage***"
            : "", (helix_idx[idx][6]) ? "  ***parallel***" : "",
            (helix_idx[idx][7]) ? "  ***intra-chain direction reverse***" : "");
}

void write_bestpairs(long num_bp, long **base_pairs, long *bp_idx, char *bseq,
                     long **seidx, char **AtomName, char **ResName, char *ChainID,
                     long *ResSeq, char **Miscs, double **xyz, double **orien,
                     double **org, long **htm_water)
/* write out a PDB file containing all best pairs oriented with mean base-pair normal */
{
    char b1[BUF512], b2[BUF512], idmsg[BUF512], wc[4];
    double morg[4], morien[10];
    long i, ia, ib, j, k;
    long inum_base = 2, ivec[3];  /* for middle frame orientation */
    FILE *mfp, *rframe;

    mfp = open_file(BESTP_FILE, "w");
    rframe = open_file(REF_FILE, "w");
    fprintf(rframe, "%5ld base-pairs\n", num_bp);
    for (i = 1; i <= num_bp; i++) {
        k = bp_idx[i];
        ia = base_pairs[k][1];
        ib = base_pairs[k][2];
        j = seidx[ia][1];
        base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[ia], 1, b1);
        j = seidx[ib][1];
        base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[ib], 2, b2);
        wc3(base_pairs[k], wc);
        sprintf(idmsg, "%s-%s-%s", b1, wc, b2);
        fprintf(mfp, "REMARK    Section #%4.4ld %s\n", i, idmsg);
        fprintf(mfp, "REMARK    %s\n", X3DNA_VER);
        ivec[1] = ia;
        ivec[2] = ib;
        pair2mst(inum_base, ivec, AtomName, ResName, ChainID, ResSeq, Miscs, xyz,
                 orien, org, seidx, morien, morg, htm_water, mfp);
        fprintf(mfp, "END\n");

        fprintf(rframe, "... %5ld %c%c%c ...\n", i, bseq[ia], wc[2], bseq[ib]);
        write_fpmst(morg, morien, rframe);
    }
    close_file(mfp);
    close_file(rframe);
}

void write_fpmst(double *morg, double *morien, FILE * rframe)
{
    static char *format = "%10.4f%10.4f%10.4f\n";
    long i, j;

    fprintf(rframe, format, morg[1], morg[2], morg[3]);
    for (i = 1; i <= 3; i++) {
        j = (i - 1) * 3;
        fprintf(rframe, format, morien[j + 1], morien[j + 2], morien[j + 3]);
    }
}

void write_helix(long num_helix, long **helix_idx, long *bp_idx, long **seidx,
                 char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                 char **Miscs, double **xyz, long **base_pairs, long **htm_water)
/* write out helical regions in a multiple structure PDB file */
{
    long i, inum, j, k, m, n, tnum_res, num_residue;
    long *ivec, *ivect;
    FILE *mfp;

    num_residue = htm_water[2][0];
    ivec = lvector(1, num_residue);
    ivect = lvector(1, num_residue);
    mfp = open_file(HLXREG_FILE, "w");

    for (i = 1; i <= num_helix; i++) {
        inum = 0;
        fprintf(mfp, "REMARK    Section #%4.4ld %ld base-pairs", i, helix_idx[i][3]);
        helix_info(helix_idx, i, mfp);
        fprintf(mfp, "REMARK    %s\n", X3DNA_VER);
        k = 0;  /* residue index in the helical region */
        for (n = 1; n <= 2; n++) {  /* two strands */
            for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
                if (n == 2 && !helix_idx[i][6])  /* anti-parallel */
                    m = base_pairs[bp_idx[helix_idx[i][2] - j + helix_idx[i][1]]][n];
                else
                    m = base_pairs[bp_idx[j]][n];
                k++;
                ivec[k] = m;
            }
        }
        tnum_res = attached_residues(k, ivec, ivect, seidx, xyz, htm_water);
        for (j = 1; j <= tnum_res; j++) {
            m = ivect[j];
            pdb_record(seidx[m][1], seidx[m][2], &inum, 0, AtomName, ResName,
                       ChainID, ResSeq, xyz, Miscs, mfp);
        }
        fprintf(mfp, "END\n");
    }

    close_file(mfp);
    free_lvector(ivec, 1, num_residue);
    free_lvector(ivect, 1, num_residue);
}

void no_basepairs(char *pdbfile, char *outfile, char *parfile)
/* for non-base-pairs found in structure "pdbfile" */
{
    FILE *fp;

    fprintf(stderr, "no base-pairs found for this structure\n");
    fp = open_file(outfile, "w");
    fprintf(fp, "%s\n", pdbfile);
    fprintf(fp, "%s.out\n", parfile);
    fprintf(fp, "    2         # duplex\n");
    fprintf(fp, "    0         # number of base-pairs\n");
    close_file(fp);
}

void duplex(long num, long num_residue, char *bseq, long **seidx, long *RY,
            char **AtomName, char **ResName, char *ChainID, long *ResSeq, char **Miscs,
            double **xyz, long curves, long divide, long hetatm, long pairs,
            long detailed, long waters, char *pdbfile, char *outfile, char *parfile)
{
    double BPRS[7];
    double **orien, **org, **NC1xyz, **o3_p;
    long i, num_bp, num_helix = 1;
    long nout = PSTNUM - 1, nout_p1 = PSTNUM;
    long *idx, *bp_idx, *helix_marker, **htm_water;
    long **base_pairs, **helix_idx, **ring_atom;

    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    NC1xyz = dmatrix(1, num_residue, 1, 7);  /* RN9/YN1 & C1' atomic coordinates */
    o3_p = dmatrix(1, num_residue, 1, 8);  /* O3'/P atomic coordinates */

    idx = lvector(1, num);
    atom_idx(num, AtomName, NULL, idx);

    htm_water = lmatrix(1, 4, 0, num);  /* HETATM and water index */
    init_htm_water(waters, num, num_residue, idx, htm_water);
    identify_htw(num_residue, seidx, RY, AtomName, ResName, ChainID, ResSeq,
                 Miscs, xyz, htm_water);

    base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
              ResSeq, Miscs, xyz, orien, org, NC1xyz, o3_p, BPRS);

    /* 1-9 ring atom index, 10 # of ring atoms, 11-19 first level */
    ring_atom = lmatrix(1, num_residue, 1, 19);
    ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, ring_atom);

    if (pairs) {  /* all possible pairs */
        all_pairs(num_residue, RY, NC1xyz, orien, org, BPRS, seidx, xyz, idx,
                  ring_atom, AtomName, ResName, ChainID, ResSeq, Miscs, bseq,
                  hetatm, htm_water, pdbfile, outfile);
        goto ALL_PAIRS;  /* to clean up */
    }

    /* find best base-pairs */
    base_pairs = lmatrix(1, num_residue, 1, nout_p1);
    /* base_pairs[][PSTNUM]:
     *      i, j, bpid, d, dv, angle, dNN, dsum, bp-org,
     * col# 1  2    3   4  5     6     7     8    9-11
     *       x1,   y1,   z1,   x2,   y2,   z2
     *      12-14 15-17 18-20 21-23 24-26 27-29
     * i.e., add one more column for i from pair_stat [best_pair]
     */
    num_bp = find_bestpair(nout, base_pairs, num_residue, bseq, seidx, RY, AtomName,
                           xyz, idx, orien, org, NC1xyz, ring_atom, BPRS);
    if (!num_bp) {
        no_basepairs(pdbfile, outfile, parfile);
        goto NO_BASE_PAIR;  /* to clean up */
    }

    bp_idx = lvector(1, num_bp);
    helix_marker = lvector(1, num_bp);
    helix_idx = lmatrix(1, num_bp, 1, 7);
    re_ordering(num_bp, base_pairs, bp_idx, helix_marker, helix_idx, BPRS,
                &num_helix, o3_p, bseq, seidx, ResName, ChainID, ResSeq, Miscs);
    write_bestpairs(num_bp, base_pairs, bp_idx, bseq, seidx, AtomName, ResName, ChainID,
                    ResSeq, Miscs, xyz, orien, org, htm_water);
    write_helix(num_helix, helix_idx, bp_idx, seidx, AtomName, ResName, ChainID, ResSeq,
                Miscs, xyz, base_pairs, htm_water);

    if (curves) {  /* generate Curves input file */
        if (divide && num_helix > 1)
            for (i = 1; i <= num_helix; i++)
                curves_input(i, helix_idx[i][1], helix_idx[i][2], helix_idx[i][3],
                             pdbfile, outfile, parfile, bp_idx, base_pairs,
                             helix_idx[i][4], helix_idx[i][6]);
        else
            curves_input(0, 1, num_bp, num_bp, pdbfile, outfile, parfile,
                         bp_idx, base_pairs, helix_idx[1][4], helix_idx[1][6]);
    } else {  /* 3DNA/CEHS input file */
        col_helices(num_helix, helix_idx, bp_idx, base_pairs, seidx, pdbfile,
                    ChainID, ResSeq);
        if (divide && num_helix > 1)
            for (i = 1; i <= num_helix; i++)
                x3dna_input(i, helix_idx[i][1], helix_idx[i][2], helix_idx[i][3],
                            pdbfile, outfile, parfile, hetatm, bp_idx, helix_marker,
                            helix_idx, base_pairs, seidx, ResName, ChainID, ResSeq,
                            Miscs, bseq, BPRS, detailed);
        else  /* overall */
            x3dna_input(0, 1, num_bp, num_bp, pdbfile, outfile, parfile, hetatm,
                        bp_idx, helix_marker, helix_idx, base_pairs, seidx, ResName,
                        ChainID, ResSeq, Miscs, bseq, BPRS, detailed);
    }
    free_lvector(bp_idx, 1, num_bp);
    free_lvector(helix_marker, 1, num_bp);
    free_lmatrix(helix_idx, 1, num_bp, 1, 7);

  NO_BASE_PAIR:
    free_lmatrix(base_pairs, 1, num_residue, 1, nout_p1);

  ALL_PAIRS:
    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);
    free_dmatrix(NC1xyz, 1, num_residue, 1, 7);
    free_dmatrix(o3_p, 1, num_residue, 1, 8);
    free_lvector(idx, 1, num);
    free_lmatrix(htm_water, 1, 4, 0, num);
    free_lmatrix(ring_atom, 1, num_residue, 1, 19);
}

long find_bestpair(long nout, long **base_pairs, long num_residue, char *bseq,
                   long **seidx, long *RY, char **AtomName, double **xyz, long *idx,
                   double **orien, double **org, double **NC1xyz, long **ring_atom,
                   double *BPRS)
/* get the best pair listing */
{
    long i, j, num1 = 0, num2 = 1, num_bp = 0, *matched_idx;
    long pair_istat[PSTNUM], pair_jstat[PSTNUM];

    matched_idx = lvector(1, num_residue);
    while (num1 < num2) {
        num1 = num2;
        for (i = 1; i <= num_residue; i++) {
            if (RY[i] < 0 || matched_idx[i])
                continue;  /* non-base or paired base */
            best_pair(i, num_residue, RY, seidx, xyz, idx, NC1xyz, matched_idx,
                      orien, org, ring_atom, AtomName, bseq, BPRS, pair_istat);
            if (pair_istat[1]) {  /* with paired base */
                best_pair(pair_istat[1], num_residue, RY, seidx, xyz, idx,
                          NC1xyz, matched_idx, orien, org, ring_atom, AtomName,
                          bseq, BPRS, pair_jstat);
                if (i == pair_jstat[1]) {  /* best match between i && pair_istat[1] */
                    matched_idx[i] = 1;
                    matched_idx[pair_istat[1]] = 1;
                    base_pairs[++num_bp][1] = i;
                    for (j = 1; j <= nout; j++)
                        base_pairs[num_bp][j + 1] = pair_istat[j];
                }
            }
        }
        num2 = 0;
        for (i = 1; i <= num_residue; i++)
            if (matched_idx[i])
                num2++;
    }
    free_lvector(matched_idx, 1, num_residue);

    return num_bp;
}

void col_helices(long num_helix, long **helix_idx, long *bp_idx, long **base_pairs,
                 long **seidx, char *pdbfile, char *ChainID, long *ResSeq)
/* RasMol script for coloring the structure by chains or helical regions */
{
    static char *col_code[9] = { "violet", "red", "green", "blue", "yellow",
        "cyan", "magenta", "orange", "purple"
    };
    long i, ia, ib, ic, j, k;
    FILE *fpc, *fph;

    fpc = open_file(COLCHN_FILE, "w");
    fprintf(fpc, "zap\nload hel_regions.pdb\n");
    fprintf(fpc, "# load %s\n", pdbfile);
    fprintf(fpc, "# restrict not (protein or water)\n");
    fprintf(fpc, "\n");

    fph = open_file(COLHLX_FILE, "w");
    fprintf(fph, "zap\nload hel_regions.pdb\n");
    fprintf(fph, "# load %s\n", pdbfile);
    fprintf(fph, "# restrict not (protein or water)\n");
    fprintf(fph, "\n");

    for (i = 1; i <= num_helix; i++) {
        ic = i % 9;  /* color code */
        fprintf(fph, "\n#------Helix #%ld, color: %s------\n", i, col_code[ic]);
        for (j = helix_idx[i][1]; j <= helix_idx[i][2]; j++) {
            k = bp_idx[j];
            ia = seidx[base_pairs[k][1]][1];
            ib = seidx[base_pairs[k][2]][1];
            fprintf(fpc, "select %ld:%c\n", ResSeq[ia], ChainID[ia]);
            fprintf(fpc, "color %s\n", col_code[1]);  /* strand I */
            fprintf(fpc, "select %ld:%c\n", ResSeq[ib], ChainID[ib]);
            fprintf(fpc, "color %s\n", col_code[2]);  /* strand II */
            fprintf(fph, "select %ld:%c, %ld:%c\n", ResSeq[ia], ChainID[ia],
                    ResSeq[ib], ChainID[ib]);
            fprintf(fph, "color %s\n", col_code[ic]);
        }
    }

    fprintf(fpc, "\nselect all\n");
    close_file(fpc);
    fprintf(fph, "\nselect all\n");
    close_file(fph);
}

void x3dna_input(long idx, long start_num, long end_num, long nbp, char *pdbfile,
                 char *outfile, char *parfile, long hetatm, long *bp_idx,
                 long *helix_marker, long **helix_idx, long **base_pairs, long **seidx,
                 char **ResName, char *ChainID, long *ResSeq, char **Miscs, char *bseq,
                 double *BPRS, long detailed)
/* print out base-pairing information in a format suitable for analyze/cehs */
{
    char b1[BUF512], b2[BUF512], marker_c, wc[4];
    char outfile_new[BUF512], parfile_new[BUF512];
    double x[4];
    long i, i_order, j, j_order, k, marker_n;
    long num_helix = 0, num_1bp = 0, num_nwc = 0;
    FILE *fp;

    if (idx) {  /* structure taken as helical regions */
        sprintf(outfile_new, "%s_%4.4ld", outfile, idx);
        sprintf(parfile_new, "%s_%4.4ld", parfile, idx);
    } else {  /* structure taken as a whole */
        strcpy(outfile_new, outfile);
        strcpy(parfile_new, parfile);
    }

    fp = open_file(outfile_new, "w");
    fprintf(fp, "%s\n", pdbfile);
    fprintf(fp, "%s.out\n", parfile_new);
    fprintf(fp, "    2         # duplex\n");
    fprintf(fp, "%5ld         # number of base-pairs\n", nbp);
    fprintf(fp, "    1%5ld    # explicit bp numbering/hetero atoms\n", hetatm);

    for (i = start_num; i <= end_num; i++) {
        k = bp_idx[i];
        i_order = base_pairs[k][1];
        j_order = base_pairs[k][2];
        if (base_pairs[k][3] != 2)
            num_nwc++;
        wc3(base_pairs[k], wc);
        j = seidx[i_order][1];
        base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[i_order], 1, b1);
        j = seidx[j_order][1];
        base_str(ChainID[j], ResSeq[j], Miscs[j], ResName[j], bseq[j_order], 2, b2);

        marker_n = 0;
        marker_c = '|';
        if (helix_marker[i]) {  /* change of helix regions */
            ++num_helix;
            if ((!idx && helix_idx[num_helix][3] == 1) || (idx && helix_idx[idx][3] == 1)) {
                num_1bp++;
                marker_n = 1;
                marker_c = '+';
            } else if (i != end_num) {
                marker_n = 9;
                marker_c = 'x';
            }
        }
        fprintf(fp, "%5ld%5ld%3ld #%5ld %c %s-%s-%s", i_order, j_order, marker_n,
                i - start_num + 1, marker_c, b1, wc, b2);
        for (j = 4; j <= 8; j++)
            fprintf(fp, "%6.2f", base_pairs[k][j] / MFACTOR);
        if (detailed) {  /* print out origin + x-axis */
            for (j = 9; j <= 11; j++)
                fprintf(fp, "%8.2f", base_pairs[k][j] / MFACTOR);
            for (j = 1; j <= 3; j++)
                x[j] = (base_pairs[k][j + 11] + (base_pairs[k][j + 20])) / MFACTOR;
            vec_norm(x);
            for (j = 1; j <= 3; j++)
                fprintf(fp, "%8.2f", x[j]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "##### ");
    print_bp_crit(BPRS, NULL, fp);
    fprintf(fp, "##### %ld non-Watson-Crick base-pair%s", num_nwc,
            (num_nwc == 1) ? "" : "s");
    (num_helix == 1) ? strcpy(b1, "x") : strcpy(b1, "ces");
    fprintf(fp, ", and %ld heli%s", num_helix, b1);
    fprintf(fp, " (%ld isolated bp%s)\n", num_1bp, (num_1bp == 1) ? "" : "s");

    if (!idx) {  /* structure taken as a whole */
        for (i = 1; i <= num_helix; i++) {
            if (helix_idx[i][3] == 1)
                fprintf(fp, "##### Helix #%ld (%ld): %ld", i, helix_idx[i][3],
                        helix_idx[i][1]);
            else
                fprintf(fp, "##### Helix #%ld (%ld): %ld - %ld", i,
                        helix_idx[i][3], helix_idx[i][1], helix_idx[i][2]);
            helix_info(helix_idx, i, fp);
        }
    } else {  /* for each helical region */
        if (nbp == 1)
            fprintf(fp, "##### Helix #1 (%ld): %ld", nbp, nbp);
        else
            fprintf(fp, "##### Helix #1 (%ld): 1 - %ld", nbp, nbp);
        helix_info(helix_idx, idx, fp);
    }

    close_file(fp);
}

void curves_input(long idx, long start_num, long end_num, long nbp, char *pdbfile,
                  char *outfile, char *parfile, long *bp_idx, long **base_pairs,
                  long zdna, long parallel)
/* print out base-pairing information in a format suitable for Curves */
{
    char outfile_new[BUF512], parfile_new[BUF512];
    long i, j;
    FILE *fp;

    if (idx) {  /* structure taken as helical regions */
        sprintf(outfile_new, "%s_%4.4ld", outfile, idx);
        sprintf(parfile_new, "%s_%4.4ld", parfile, idx);
    } else {  /* structure taken as a whole */
        strcpy(outfile_new, outfile);
        strcpy(parfile_new, parfile);
    }

    fp = open_file(outfile_new, "w");
    fprintf(fp, "&inp file=%s, comb=.t., fit=.t., grv=.t., %s\n"
            "     lis=%s, pdb=%s_grp, &end\n", pdbfile,
            zdna ? "dinu=.t.," : "", parfile_new, parfile_new);
    fprintf(fp, "    2%5ld%5ld    0    0\n", nbp, (parallel) ? nbp : -nbp);
    for (i = 1; i <= 2; i++) {
        for (j = start_num; j <= end_num; j++)
            fprintf(fp, "%5ld", base_pairs[bp_idx[j]][i]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "   0.0   0.0   0.0%6.1f\n", zdna ? 180.0 : 0.0);
    close_file(fp);
}
