#include "x3dna.h"

int main(int argc, char *argv[])
/* extract a structure from a multiple-structure PDB file (NMR, stacking.pdb etc) */
{
    char inpfile[BUF512], outfile[BUF512];
    long inum = 1, NMR = 0, BIOU = 0, DELH = 0;

    exstr_cmdline(argc, argv, inpfile, outfile, &inum, &NMR, &BIOU, &DELH);
    remove_file(outfile);

    if (NMR)
        extract_NMR(inpfile, outfile, inum);
    else if (BIOU)  /* assuming NMR contains only one biological unit */
        extract_BIOU(inpfile, outfile, inum);
    else  /* 3DNA structure unit */
        extract_3DNA(inpfile, outfile, inum);

    if (DELH) {
        copy_file(outfile, TMP_FILE);
        delH_pdbfile(TMP_FILE, outfile);
    }

    return 0;
}

void exstr_usage(void)
{
    help3dna("ex_str");
    contact_msg(0);
}

void exstr_cmdline(int argc, char *argv[], char *inpfile, char *outfile, long *inum,
                   long *NMR, long *BIOU, long *DELH)
{
    long i;

    if (argc < 3)
        exstr_usage();
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        upperstr(argv[i]);
        if (sscanf(argv[i], "-%ld", inum) == 1) {
            if (*inum <= 0) {
                fprintf(stderr, "Structure #: %ld, reset to 1\n", *inum);
                *inum = 1;
            }
        } else if (!strcmp(argv[i], "-NMR"))
            *NMR = 1;
        else if (!strncmp(argv[i], "-NMRB", 5))
            *NMR = 2;
        else if (!strncmp(argv[i], "-BIOU", 5))
            *BIOU = 1;
        else if (!strncmp(argv[i], "-DELH", 5))
            *DELH = 1;
        else
            exstr_usage();
    }
    if (argc == i + 2) {
        strcpy(inpfile, argv[i]);
        strcpy(outfile, argv[i + 1]);
    } else
        exstr_usage();

    if (*NMR)
        check_ifNMR(inpfile);
    if (*NMR == 2)
        *inum = best_conformer(inpfile);
}

void check_ifNMR(char *pdbfile)
/* Check if "pdbfile" is an NMR file based on EXPDTA record */
{
    char str[BUF512];
    long nmr = 0;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        upperstr(str);
        if (!strncmp(str, "ATOM", 4) || !strncmp(str, "HETATM", 6) ||
            !strncmp(str, "END", 3))
            break;
        if (!strncmp(str, "EXPDTA", 6) && strstr(str, "NMR")) {
            nmr = 1;
            break;
        }
    }
    if (!nmr)
        fprintf(stderr, "Warning: ***unlikely an NMR file***\n");
}

long best_conformer(char *pdbfile)
/* Best representative NMR conformer, use #1 if not assigned in REMARK 210.
 * Typical cases:
 *    BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE : NULL [pdb28sp]
 *    BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE : 1    [pdb2bj2]
 *    BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE : 3    [pdb2lef]
 *    BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE : 5    [pdb3php]
 *
 * For structures pdb1qby, pdb1djd & pdb1bjd, the REMARK 210 says the
 *    best representative conformer is #12 or #13, but the data file
 *    actually does not has that model. Reset inum to 1 for such cases.
 */
{
    char str[BUF512], numstr[BUF512];
    long bestrec = 0, inum;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        upperstr(str);
        if (!strncmp(str, "ATOM", 4) || !strncmp(str, "HETATM", 6) ||
            !strncmp(str, "END", 3))
            break;
        if (!strncmp(str, "REMARK 210 BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE :",
                     59)) {
            bestrec = 1;  /* has BEST REP ... record */
            if (sscanf(str + 59, "%ld", &inum) != 1) {
                fprintf(stderr, "Best representative conformer NULL: set to 1\n");
                inum = 1;
            } else {  /* Check if model inum really exists */
                sprintf(numstr, "%ld", inum);
                while (fgets(str, sizeof str, fp) != NULL)
                    if (!strncmp(str, "MODEL ", 6) && strstr(str, numstr) != NULL) {
                        bestrec = 2;
                        break;
                    }
                if (bestrec != 2) {
                    fprintf(stderr, "Best representative conformer #%ld does NOT exist"
                            ": set to 1\n", inum);
                    inum = 1;
                }
            }
            break;
        }
    }
    if (!bestrec) {
        fprintf(stderr, "No best representative conformer assigned: set to 1\n");
        inum = 1;
    }
    close_file(fp);

    return inum;
}

void extract_NMR(char *inpfile, char *outfile, long inum)
/* Extract structure unit from an NMR PDB file */
{
    char str[BUF512];
    long found_str = 0, num0;
    FILE *fpi, *fpo;

    fpi = open_file(inpfile, "r");
    fpo = open_file(outfile, "w");

    while (fgets(str, sizeof str, fpi) != NULL) {
        if (!strncmp(str, "MODEL ", 6) && (sscanf(str, "%*s %ld", &num0) == 1) &&
            num0 == inum) {
            print_pdb_title(inpfile, "*", fpo);
            found_str = 1;
            fprintf(stderr, "Model selected #%ld\n", inum);
            if (fputs(str, fpo) == EOF)
                fatal("error in writing selected structure\n");
            while (1) {
                if (fgets(str, sizeof str, fpi) != NULL) {
                    fputs(str, fpo);
                    if (!strncmp(str, "ENDMDL", 6)) {
                        found_str = 2;
                        goto FINISHED;
                    }
                } else if (feof(fpi))
                    goto FINISHED;
                else
                    fatal("error reading PDB data file\n");
            }
        }
    }

  FINISHED:
    if (!found_str) {
        fprintf(stderr, "NMR structure #: %ld not found in file %s\n", inum, inpfile);
        if (inum == 1) {  /* only for the 1st unit */
            fprintf(stderr, "    Copy the whole file: %s ===> %s\n", inpfile, outfile);
            copy_file(inpfile, outfile);
        }
    } else if (found_str != 2) {
        fprintf(stderr, "NMR structure #: %ld has no ENDMDL record in file %s\n", inum,
                inpfile);
        fputs("ENDMDL\n", fpo);
        fputs("END\n", fpo);
    }

    close_file(fpi);
    close_file(fpo);
}

void extract_3DNA(char *inpfile, char *outfile, long inum)
/* Extract a structure unit from an 3DNA output */
{
    char str[BUF512], numstr[BUF512];
    long found_str = 0;
    FILE *fpi, *fpo;

    fpi = open_file(inpfile, "r");
    fpo = open_file(outfile, "w");

    sprintf(numstr, "Section #%4.4ld", inum);
    while (fgets(str, sizeof str, fpi) != NULL) {
        if (strstr(str, numstr) != NULL) {
            found_str = 1;
            fprintf(stderr, "Structure selected #%ld\n", inum);
            if (fputs(str, fpo) == EOF)
                fatal("error in writing selected structure\n");
            while (1) {
                if (fgets(str, sizeof str, fpi) != NULL) {
                    fputs(str, fpo);
                    if (!strncmp(str, "END", 3)) {
                        found_str = 2;
                        goto FINISHED;
                    }
                } else if (feof(fpi))
                    goto FINISHED;
                else
                    fatal("error reading PDB data file\n");
            }
        }
    }

  FINISHED:
    if (!found_str) {
        fprintf(stderr, "3DNA structure #: %ld not found in file %s\n", inum, inpfile);
        if (inum == 1) {  /* only for the 1st unit */
            fprintf(stderr, "    Copy the whole file: %s ===> %s\n", inpfile, outfile);
            copy_file(inpfile, outfile);
        }
    } else if (found_str != 2) {
        fprintf(stderr, "3DNA structure #: %ld has no END record in file %s\n", inum,
                inpfile);
        fputs("END\n", fpo);
    }

    close_file(fpi);
    close_file(fpo);
}

void extract_BIOU(char *inpfile, char *outfile, long inum)
/* Extract a biological unit from an x-ray crystal structure, based on
 * REMARK 350 record */
{
    char str[BUF512], *ChainID;
    char **AtomName, **ResName, **Miscs, **chlist;
    double **xyz;
    long j, k, nb, num, *ResSeq;
    FILE *fp;

    /* read in the PDB file */
    num = number_of_atoms(inpfile, 1, "*");
    AtomName = cmatrix(1, num, 0, 4);
    ResName = cmatrix(1, num, 0, 3);
    ChainID = cvector(1, num);
    ResSeq = lvector(1, num);
    xyz = dmatrix(1, num, 1, 3);
    Miscs = cmatrix(1, num, 0, NMISC);
    read_pdb(inpfile, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, 1, "*");

    chlist = cmatrix(1, BUF512, 0, MAXCH);
    nb = biounit_cid(inpfile, chlist);

    if (!nb) {
        fprintf(stderr, "File %s does NOT contain `REMARK 350' records\n", inpfile);
        fprintf(stderr, "    Copy the whole file: %s ===> %s\n", inpfile, outfile);
        copy_file(inpfile, outfile);
    } else {
        fprintf(stderr, "Number of biological units: %ld\n", nb);
        if (inum > nb) {
            fprintf(stderr, "%ld > total # of biological unit: %ld\n", inum, nb);
            fprintf(stderr, "    Please re-select your biological unit #\n");
        } else {
            fprintf(stderr, "Biological unit selected #%ld\n", inum);
            fp = open_file(outfile, "w");
            print_pdb_title(inpfile, chlist[inum], fp);
            fprintf(fp, "REMARK 350 BIOMOLECULE: %ld\n", inum);
            k = 0;
            for (j = 1; j <= num; j++)
                if (strchr(chlist[inum], ChainID[j])) {
                    strcpy(str, Miscs[j]);
                    fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n",
                            (str[0] == 'A') ? "ATOM  " : "HETATM", ++k, AtomName[j],
                            str[1], ResName[j], ChainID[j], ResSeq[j], str[2], xyz[j][1],
                            xyz[j][2], xyz[j][3], str + 3);
                }
            fprintf(fp, "END\n");
            close_file(fp);
        }
    }
    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_cmatrix(chlist, 1, BUF512, 0, MAXCH);
}

long biounit_cid(char *pdbfile, char **chlist)
/* chain id for each biounit */
{
    char str[BUF512], *item[MAXCH];
    long i, n, num_bu = 0;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        upperstr(str);
        if (!strncmp(str, "ATOM", 4) || !strncmp(str, "HETATM", 6) ||
            !strncmp(str, "END", 3))
            break;
        if (!strncmp(str, "REMARK 350 APPLY THE FOLLOWING TO CHAINS:", 41)) {
            if (++num_bu > BUF512)
                fatal("Number of biological units exceed upper limit\n");
            n = itemize(str + 41, item, MAXCH) + 1;
            for (i = 0; i < n; i++)
                chlist[num_bu][i] = item[i][0];
            chlist[num_bu][i] = '\0';
        }
    }
    close_file(fp);

    return num_bu;
}
