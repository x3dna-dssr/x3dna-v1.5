#include "x3dna.h"

int main(int argc, char *argv[])
/* ALCHEMY schematic presentation to PS[0], XFIG[1], Raster3D[2] image */
{
    char inpfile[BUF512], outfile[BUF512];
    double scale_factor = 0.0;
    long opts[BUF512] = { 0 };  /* all initialized to zeros */
    /* 0: 0 for PS, 1 (F) for XFIG, and 2 (R) for Raster3D
     * 1: (B) frame_box
     * 2: (C) is_color
     * 3: (L) oline
     * 4: peptide (NOT used here, for compatibility with PDB2IMG)
     * 5: (A) helix_axis
     * 6: (G) ref_axis
     * 7: (I) same_faces
     * 8: (N) no_header
     * 9: (U) updown
     * 10: (M) isomM */

    alc2img_cmdline(argc, argv, inpfile, outfile, &scale_factor, opts);
    process_alc(inpfile, outfile, scale_factor, opts);

    return 0;
}

void alc2img_usage(void)
{
    help3dna("alc2img");
    contact_msg(0);
}

void alc2img_cmdline(int argc, char *argv[], char *inpfile, char *outfile,
                     double *scale_factor, long *opts)
{
    long i, j;

    if (argc < 3)
        alc2img_usage();
    for (i = 1; i < argc; i++) {
        if (*argv[i] != '-')
            break;
        upperstr(argv[i]);
        for (j = 1; j < (long) strlen(argv[i]); j++)  /* skip - */
            if (argv[i][j] == 'F')
                opts[0] = 1;
            else if (argv[i][j] == 'R')
                opts[0] = 2;
            else if (argv[i][j] == 'B')
                opts[1] = 1;
            else if (argv[i][j] == 'C')
                opts[2] = 1;
            else if (argv[i][j] == 'L')
                opts[3] = 1;
            else if (argv[i][j] == 'A')
                opts[5] = 1;
            else if (argv[i][j] == 'G')
                opts[6] = 1;
            else if (argv[i][j] == 'I')
                opts[7] = 1;
            else if (argv[i][j] == 'N')
                opts[8] = 1;
            else if (argv[i][j] == 'U')
                opts[9] = 1;
            else if (argv[i][j] == 'M')
                opts[10] = 1;
            else if (argv[i][j] == 'S') {
                if (sscanf(argv[i], "-S=%lf", scale_factor) != 1) {
                    fprintf(stderr, "wrong format for setting scale factor\n");
                    alc2img_usage();
                } else
                    break;  /* -s=factor not combined with others */
            } else
                alc2img_usage();
    }
    if (argc == i + 2) {
        strcpy(inpfile, argv[i]);
        strcpy(outfile, argv[i + 1]);
    } else
        alc2img_usage();
    if (opts[10])  /* is 4 + minor/major faces, set -i option */
        opts[7] = 1;
}
