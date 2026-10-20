#include "x3dna.h"

int main(int argc, char *argv[])
/* data file for generating a regular DNA structure by "rebuild" */
{
    char user_input[BUF512], outfile[BUF512];
    static char *Wbase = CB_LIST, *Cbase = "TGCCAA";
    char **bp_seq;
    char *format = "%8.2f";
    char *options[2] = { "-STEP", "-HELICAL" };

    double bp_par[6], sp_par[6];

    long i, j, num_bp;
    long is_helical = 0;

    FILE *fp;

    strcpy(user_input, options[0]);  /* step parameters */
    regular_dna_cmdline(argc, argv, options, user_input, outfile);

    fprintf(stderr, "\nSix base-pair parameters (Dft: 0s) in the order of: \n"
            "Shear  Stretch  Stagger  Buckle  Propeller  Opening\n");
    get_six_pars(1, bp_par);
    if (!strcmp(user_input, options[0]))
        fprintf(stderr, "\nSix step parameters in the order of: \n"
                "Shift  Slide  Rise  Tilt  Roll  Twist\n");
    else {
        fprintf(stderr, "\nSix helical parameters in the order of: \n"
                "X-disp  Y-disp  h-Rise  Incl.  Tip  h-Twist\n");
        is_helical = 1;
    }
    get_six_pars(0, sp_par);

    bp_seq = cmatrix(1, MAXBASE, 1, 2);

    bp_sequence(Wbase, Cbase, &num_bp, bp_seq);

    fp = open_file(outfile, "w");

    fprintf(fp, "%4ld base-pairs\n", num_bp);
    if (is_helical) {
        fprintf(fp, "   1  ***local base-pair & helical parameters***\n");
        fprintf(fp, "       Shear  Stretch  Stagger Buckle Prop-Tw Opening  "
                "X-disp  Y-disp  h-Rise  Incl.    Tip   h-Twist\n");
    } else {
        fprintf(fp, "   0  ***local base-pair & step parameters***\n");
        fprintf(fp, "       Shear  Stretch  Stagger Buckle Prop-Tw Opening  "
                " Shift  Slide    Rise    Tilt    Roll   Twist\n");
    }

    /* 1st base-pair */
    fprintf(fp, "%c-%c ", bp_seq[1][1], bp_seq[1][2]);
    for (i = 0; i < 6; i++)
        fprintf(fp, format, bp_par[i]);
    for (i = 0; i < 6; i++)
        fprintf(fp, format, 0.0);
    fprintf(fp, "\n");

    for (i = 2; i <= num_bp; i++) {
        fprintf(fp, "%c-%c ", bp_seq[i][1], bp_seq[i][2]);
        for (j = 0; j < 6; j++)
            fprintf(fp, format, bp_par[j]);
        for (j = 0; j < 6; j++)
            fprintf(fp, format, sp_par[j]);
        fprintf(fp, "\n");
    }

    close_file(fp);

    free_cmatrix(bp_seq, 1, MAXBASE, 1, 2);

    return 0;
}

void regular_dna_usage(void)
{
    help3dna("regular_dna");
    contact_msg(0);
}

void regular_dna_cmdline(int argc, char *argv[], char *options[2], char *user_input,
                         char *outfile)
{
    long i, j;

    if (argc < 2)
        regular_dna_usage();
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        upperstr(argv[i]);
        for (j = 0; j < 2; j++)
            if (!strcmp(argv[i], options[j])) {
                strcpy(user_input, argv[i]);
                break;
            }
        if (j >= 2)
            regular_dna_usage();
    }
    if (argc == i + 1)
        strcpy(outfile, argv[i]);
    else
        regular_dna_usage();
}

void get_six_pars(long is_bp, double *pars)
{
    char str[BUF512];
    long i, ik;

    while (1)
        if (fgets(str, sizeof str, stdin) != NULL) {
            for (i = 0; i < (long) strlen(str); i++)
                if (str[i] == ',')
                    str[i] = ' ';  /* change comma to space */
            ik = sscanf(str, "%lf %lf %lf %lf %lf %lf", &pars[0], &pars[1],
                        &pars[2], &pars[3], &pars[4], &pars[5]);
            if (is_bp && (!ik || ik == EOF)) {
                for (i = 0; i < 6; i++)
                    pars[i] = 0.0;
                return;
            } else if (ik == 6)
                return;
            else
                fprintf(stderr, "please input six parameters\n");
        } else
            fatal("error reading six parameters\n");
}
