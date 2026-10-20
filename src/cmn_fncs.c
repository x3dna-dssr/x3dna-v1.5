#include "x3dna.h"

void print_used_time(time_t time0)
{
    char str[BUF512];
    double dtime;
    long scnds;

    dtime = difftime(time(NULL), time0);
    sprintf(str, "%.0f", dtime);
    if (sscanf(str, "%ld", &scnds) != 1)
        fatal("Wrong time format: %s\n", str);
    fprintf(stderr, "\nTime used: %ld second%s\n", scnds, (scnds == 1) ? "" : "s");
}

void parcat(char *str, double par, char *format, char *bstr)
{
    char temp[BUF512];

    if (par > EMPTY_CRITERION) {
        sprintf(temp, format, par);
        strcat(str, temp);
    } else
        strcat(str, bstr);
}

void print_bp_crit(double *BPRS, char *HB_ATOM, FILE * fp)
/* print out base-pair criteria used */
{
    long i;

    fprintf(fp, "Base-pair criteria used: ");
    for (i = 1; i <= 6; i++)
        fprintf(fp, "%6.2f", BPRS[i]);
    fprintf(fp, " %s\n", (HB_ATOM == NULL) ? "" : HB_ATOM);
}

static int endofline(FILE * fp, int c)
/* endofline: check for and consume \r, \n, \r\n, or EOF */
{
    int eol;

    eol = (c == '\r' || c == '\n');
    if (c == '\r') {
        c = getc(fp);
        if (c != '\n' && c != EOF)
            ungetc(c, fp);  /* read too far; put c back */
    }
    return eol;
}

static char *enlarge_cline(long *maxline, char *line)
/* for a 0-index char-array */
{
    char *newline;

    *maxline *= SFACTOR;

    if ((newline = (char *) realloc(line, (*maxline) * sizeof(char))) == NULL)
        fatal("realloc failure in enlarge_cline()\n");

    return newline;
}

char *my_getline(FILE * fp)
/* read line of arbitrary length from fp */
{
    int c;
    long i, maxline = BUF512;
    char *line = NULL;

    line = (char *) malloc(maxline * sizeof(char));
    if (line == NULL)
        fatal("out of memory in my_getline()\n");

    for (i = 0; (c = getc(fp)) != EOF && !endofline(fp, c); i++) {
        if (i >= maxline - 1)
            line = enlarge_cline(&maxline, line);
        line[i] = c;
    }
    line[i] = '\0';
    return (c == EOF && i == 0) ? NULL : line;
}

long csplit(char *str, char *item[], long itemsize, char sepc)
/* strtok is GREEDY, meaning empty field will be ignored:
 * e.g. \t\t. This function is created to fix the problem. */
{
    char *p0, *p;
    long nitem = 0;

    if (str[0] == '\0')
        return nitem;

    p = str;
    while (*p) {  /* change '"' to ' ' */
        if (*p == '"')
            *p = ' ';
        p++;
    }

    p0 = str;
    while ((p = strchr(p0, sepc)) != NULL) {
        *p = '\0';
        item[++nitem] = trim(p0);
        if (nitem >= itemsize)
            return itemsize;
        p0 = p + 1;
    }
    item[++nitem] = trim(p0);
    return nitem;
}

char *trim(char *a)
/* trim leading and trailing white spaces */
{
    int c;

    while (isspace(c = *a))
        a++;
    for (c = strlen(a) - 1; c >= 0; c--)
        if (!isspace(a[c]))
            break;
    a[c + 1] = '\0';

    return a;
}

char *ltrim(char *a)
/* trim leading white spaces */
{
    int c;

    while (isspace(c = *a))
        a++;

    return a;
}

char *rtrim(char *a)
/* trim trailing white spaces */
{
    int c;

    for (c = strlen(a) - 1; c >= 0; c--)
        if (!isspace(a[c]))
            break;
    a[c + 1] = '\0';

    return a;
}

long itemize(char *str, char *item[], long itemsize)
/* itemize a string into its components */
{
    static char *sep_chars = " ,\t\n";  /* space, comma, tab & newline */
    char *p;
    long nitem = 0;

    for (p = strtok(str, sep_chars); p != NULL; p = strtok(NULL, sep_chars))
        if (nitem > itemsize) {
            fprintf(stderr, "ignoring items following %ldth\n", nitem);
            return itemsize;
        } else
            item[nitem++] = p;
    if (!nitem)
        fatal("empty input line not allowed\n");
    return nitem - 1;  /* base-pair does not count */
}

void refs_right_left(long bnum, double **orien, double **org, double **r1, double *o1,
                     double **r2, double *o2)
/* right-left base in a pair */
{
    long ioffset3, ioffset9;

    ioffset3 = (bnum - 1) * 3;
    ioffset9 = (bnum - 1) * 9;

    cpxyz(org[2] + ioffset3, o1);
    cpxyz(org[1] + ioffset3, o2);
    orien2mst(orien[2], ioffset9, r1);
    orien2mst(orien[1], ioffset9, r2);
}

void refs_i_j(long b1, long b2, double *bp_orien, double *bp_org, double **r1,
              double *o1, double **r2, double *o2)
/* base (-pair) step */
{
    ref_frame_i(b1, bp_orien, bp_org, r1, o1);
    ref_frame_i(b2, bp_orien, bp_org, r2, o2);
}

void ref_frame_i(long bnum, double *bp_orien, double *bp_org, double **r, double *o)
/* get the reference frame for "bnum" */
{
    long ioffset3, ioffset9;

    ioffset3 = (bnum - 1) * 3;
    ioffset9 = (bnum - 1) * 9;

    cpxyz(bp_org + ioffset3, o);
    orien2mst(bp_orien, ioffset9, r);
}

void mst2orien(double *orien_vec, long ioffset, double **mst)
/* set orientation vector with ioffset from a 3-by-3 column-wise matrix */
{
    long i, ik, j;

    for (i = 1; i <= 3; i++) {
        ik = ioffset + (i - 1) * 3;
        for (j = 1; j <= 3; j++)
            orien_vec[ik + j] = mst[j][i];
    }
}

void orien2mst(double *orien_vec, long ioffset, double **mst)
/* get 3-by-3 column-wise matrix from orientation vector with ioffset */
{
    long i, ik, j;

    for (i = 1; i <= 3; i++) {
        ik = ioffset + (i - 1) * 3;
        for (j = 1; j <= 3; j++)
            mst[j][i] = orien_vec[ik + j];
    }
}

void x_y_z_2_mtx(double *x, double *y, double *z, double **mtx)
/* given 3 vectors (x, y & z), get matrix (column-wise) */
{
    long i;

    for (i = 1; i <= 3; i++) {
        mtx[i][1] = x[i];
        mtx[i][2] = y[i];
        mtx[i][3] = z[i];
    }
}

void cehs_average(long inum_base, long *ivec, double **orien, double **org, double **mst,
                  double *morg)
/* Averaging method 1 (CEHS way): always used for a single base-pair */
{
    double pars[7], orgn[4], **bi, **mstn;
    long i, ik, j;

    mstn = dmatrix(1, 3, 1, 3);
    bi = dmatrix(1, 3, 1, 3);

    ik = ivec[1];  /* middle frame initialized to 1st base */
    cpxyz(org[ik], morg);
    orien2mst(orien[ik], 0, mst);
    for (i = 1; i <= inum_base; i++) {  /* for each residue */
        ik = ivec[i];
        orien2mst(orien[ik], 0, bi);
        if (dot(&orien[ik][6], &orien[ivec[1]][6]) < 0.0)
            for (j = 1; j <= 3; j++) {  /* reverse y- and z-axes: column vector */
                bi[j][2] = -bi[j][2];
                bi[j][3] = -bi[j][3];
            }
        bpstep_par(bi, org[ik], mst, morg, pars, mstn, orgn);
        copy_matrix(mstn, 3, 3, mst);
        cpxyz(orgn, morg);
    }

    free_dmatrix(mstn, 1, 3, 1, 3);
    free_dmatrix(bi, 1, 3, 1, 3);
}

void geom_average(long inum_base, long *ivec, double **orien, double **org, double **mst,
                  double *morg)
/* Geometrical averaging (z-, x-, then y-axes): always used for > 2 bases */
{
    long ap, i, ik, j;
    double sx[4], sy[4], sz[4];

    for (i = 1; i <= 3; i++)
        morg[i] = sx[i] = sz[i] = 0.0;  /* y-axis follows x- and z-axes */

    for (i = 1; i <= inum_base; i++) {  /* for each residue */
        ik = ivec[i];
        ap = dot(&orien[ik][6], &orien[ivec[1]][6]) < 0.0;
        for (j = 1; j <= 3; j++) {
            morg[j] += org[ik][j];
            sx[j] += orien[ik][j];  /* x-axis */
            sz[j] += (ap) ? -orien[ik][6 + j] : orien[ik][6 + j];
        }
    }
    vec_norm(sz);
    vec_orth(sx, sz);  /* x-axis normalized */
    cross(sz, sx, sy);

    x_y_z_2_mtx(sx, sy, sz, mst);
    for (i = 1; i <= 3; i++)
        morg[i] /= inum_base;
}

void pair2mst(long inum_base, long *ivec, char **AtomName, char **ResName, char *ChainID,
              long *ResSeq, char **Miscs, double **xyz, double **orien, double **org,
              long **seidx, double *mst_orien, double *mst_org, long **htm_water,
              FILE * fp)
/* Write a base-pair or a multiplet with reference to its middle-frame:
 * There are 2 methods to calculate the middle-frame with > 2 bases:
 * [1] Get the middle-frame between 1 & 2 (m1_2) in the usual way
 *     then the middle-frame between m1_2 & 3 (m1_2_3), etc.
 *     The final result DEPENDS ON the order, i.e., m1_3_2 is different
 *     from m1_2_3, both in position and in orientation!
 * [2] Get the mean z-axis of ALL BASES as the middle-frame z-axis, and
 *     similarly for the x-axis but with a correction for orthogonality.
 *     y-axis follows a right-handed rule. Applying y-first will give
 *     slightly different result; for pan-anti tetrad which is symmetrical,
 *     averaging will give near zero values: ===> big error!
 * So [1] is used all the time ...... */
{
    double morg[4], **mst, **xyz_residue;
    long i, ik, j, m, tnum_res, inum = 0;
    long ivec2[BUF512];

    tnum_res = attached_residues(inum_base, ivec, ivec2, seidx, xyz, htm_water);

    mst = dmatrix(1, 3, 1, 3);
    init_dvector(morg, 1, 3, 0.0);
    cehs_average(inum_base, ivec, orien, org, mst, morg);

    if (mst_orien != NULL) {  /* return middle frame */
        cpxyz(morg, mst_org);
        mst2orien(mst_orien, 0, mst);
    }

    xyz_residue = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
    for (i = 1; i <= tnum_res; i++) {  /* for each residue + attachments */
        ik = ivec2[i];
        for (j = seidx[ik][1]; j <= seidx[ik][2]; j++) {
            m = j - seidx[ik][1] + 1;
            cpxyz(xyz[j], xyz_residue[m]);
        }
        change_xyz(0, morg, mst, seidx[ik][2] - seidx[ik][1] + 1, xyz_residue);
        pdb_record(seidx[ik][1], seidx[ik][2], &inum, 1, AtomName, ResName,
                   ChainID, ResSeq, xyz_residue, Miscs, fp);
    }

    free_dmatrix(mst, 1, 3, 1, 3);
    free_dmatrix(xyz_residue, 1, NUM_RESIDUE_ATOMS, 1, 3);
}

void get_chi_angle(long num_residue, long *RY, char *bseq, long **seidx, double **xyz,
                   char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                   char **Miscs, double *chi, long **idxCN)
/* get the chi torsion angle for each base residue */
{
    char idmsg[BUF512];
    long i, ib, ie, idx[5], j;
    double **xyz4;

    xyz4 = dmatrix(1, 4, 1, 3);
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0)
            continue;  /* no base residue */
        ib = seidx[i][1];
        ie = seidx[i][2];
        sprintf(idmsg, ": residue name %s, chain %c, number %4ld%c",
                ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2]);
        idx[1] = find_1st_atom(" O4'", AtomName, ib, ie, idmsg);
        idx[2] = find_1st_atom(" C1'", AtomName, ib, ie, idmsg);
        if (RY[i] == 1) {
            idx[3] = find_1st_atom(" N9 ", AtomName, ib, ie, idmsg);
            idx[4] = find_1st_atom(" C4 ", AtomName, ib, ie, idmsg);
        } else {
            if (bseq[i] == 'P' || bseq[i] == 'p') {
                idx[3] = find_1st_atom(" C5 ", AtomName, ib, ie, idmsg);
                idx[4] = find_1st_atom(" C4 ", AtomName, ib, ie, idmsg);
            } else {
                idx[3] = find_1st_atom(" N1 ", AtomName, ib, ie, idmsg);
                idx[4] = find_1st_atom(" C2 ", AtomName, ib, ie, idmsg);
            }
        }
        if (idxCN != NULL) {  /* get C1* and RN9/YN1 indices */
            idxCN[i][1] = idx[2];
            idxCN[i][2] = idx[3];
        }
        for (j = 1; j <= 4; j++) {
            if (!idx[j])
                break;
            cpxyz(xyz[idx[j]], xyz4[j]);
        }
        chi[i] = (j > 4) ? torsion(xyz4) : EMPTY_NUMBER;
    }
    free_dmatrix(xyz4, 1, 4, 1, 3);
}

FILE *open_file(char *filename, char *filemode)
{
    FILE *fp;

    errno = 0;
    if (filename == NULL)
        filename = "\0";
    if (!strcmp(filename, "stdin"))
        fp = stdin;
    else if (!strcmp(filename, "stdout"))
        fp = stdout;
    else if (!strcmp(filename, "stderr"))
        fp = stderr;
    else {
        fp = fopen(filename, filemode);
        if (fp == NULL)
            fatal("open_file <%s> failed: %s\n", filename, strerror(errno));
    }
    return fp;
}

long close_file(FILE * fp)
{
    long i;

    if (fp == NULL || fp == stdin || fp == stdout || fp == stderr)
        return 0;
    errno = 0;
    i = fclose(fp);
    if (i == EOF)
        fatal("close_file failed: %s\n", strerror(errno));
    return i;
}

long exist_file(char *filename)
{
    long status;
    FILE *fp;

    fp = fopen(filename, "r");
    status = (fp != NULL) ? 1 : 0;
    close_file(fp);
    return status;
}

void remove_file(char *filename)
{
    if (!exist_file(filename))
        return;

    if (remove(filename))
        fatal("can not remove file: %s\n", filename);
}

void rename_file(char *src, char *dst)
{
    if (!exist_file(src))
        fatal("file to be renamed <%s> does NOT exist\n", src);

    if (rename(src, dst))
        fatal("can not rename file <%s> to <%s>\n", src, dst);
}

void copy_file(char *src, char *dst)
/* copy file contents from src to dst */
{
    char str[BUF512];
    size_t num_bytes;
    FILE *fpi, *fpo;

    if (!strcmp(src, dst)) {
        fprintf(stderr, "same source/destination file name: <%s>\n", src);
        return;
    }

    fpi = open_file(src, "rb");
    fpo = open_file(dst, "wb");
    while (!feof(fpi)) {
        num_bytes = fread(str, 1, BUF512, fpi);
        if (fwrite(str, 1, num_bytes, fpo) != num_bytes)
            fatal("file copy error from <%s> to <%s>\n", src, dst);
    }

    close_file(fpi);
    close_file(fpo);
}

void cat_file(char *src, char *dst)
/* cat file contents from src to dst */
{
    char str[BUF512];
    size_t num_bytes;
    FILE *fpi, *fpo;

    fpi = open_file(src, "rb");
    fpo = open_file(dst, "ab");
    while (!feof(fpi)) {
        num_bytes = fread(str, 1, BUF512, fpi);
        if (fwrite(str, 1, num_bytes, fpo) != num_bytes)
            fatal("file concatenation error from <%s> to <%s>\n", src, dst);
    }

    close_file(fpi);
    close_file(fpo);
}

long upperstr(char *a)
/* change to upper case, and return string length */
{
    long nlen = 0;

    while (*a) {
        nlen++;
        if (islower((int) *a))
            *a = toupper(*a);
        a++;
    }
    return nlen;
}

long lowerstr(char *a)
/* change to lower case, and return string length */
{
    long nlen = 0;

    while (*a) {
        nlen++;
        if (isupper((int) *a))
            *a = tolower(*a);
        a++;
    }
    return nlen;
}

char *my_strdup(char *src)
{
    char *dst;

    dst = cvector(0, strlen(src));
    strcpy(dst, src);

    return dst;
}

void print_sep(FILE * fp, char x, long n)
/* print char 'x' n-times to stream fp */
{
    long i;

    for (i = 1; i <= n; i++)
        if (fputc(x, fp) == EOF)
            fatal("error writing characters to the stream\n");
    if (fputc('\n', fp) == EOF)
        fatal("error writing '\n' to the stream\n");
}

void check_slash(char *BDIR)
/* check if '/' exist at the end of the original BDIR */
{
    char *pchar;
    long n;

    pchar = strrchr(BDIR, '/');
    n = strlen(BDIR);
    if (pchar - BDIR != n - 1) {
        BDIR[n] = '/';
        BDIR[n + 1] = '\0';
    }
}

void delete_end_slash(char *str)
/* delete '/' at the end of the str */
{
    char *pchar;
    long n;

    pchar = strrchr(str, '/');
    n = strlen(str);
    if (pchar - str == n - 1)
        str[n - 1] = '\0';
}

char *basename(char *str)
/* return a pointer 1-char following the last slash or to str w/o '/' */
{
    char *p_lslash;

    p_lslash = strrchr(str, '/');
    if (p_lslash == NULL)
        return str;
    else
        return p_lslash + 1;
}

long lround(double d)
{
    return (long) ((d > 0.0) ? d + 0.5 : d - 0.5);
}

void del_extension(char *fullname, char *okname)
/* get ride of the extension in a file name */
{
    char *pchar;
    size_t i;

    pchar = strrchr(fullname, '.');
    if (pchar == NULL)
        strcpy(okname, fullname);
    else {
        i = pchar - fullname;
        strncpy(okname, fullname, i);
        okname[i] = '\0';
    }
}

void fatal(char *fmt, ...)
/* exit program with an error message */
{
    va_list args;

    if (strlen(fmt) > 0) {
        va_start(args, fmt);
        vfprintf(stderr, fmt, args);
        va_end(args);
    }

    exit(1);
}

void print_pdb_title(char *pdbfile, char *chain_list, FILE * fp)
/* print HEADER, TITLE etc records to parameter output file for easy reference */
{
    char str[BUF512];
    static char *titles[] =
        { "HEADER", "TITLE ", "AUTHOR", "COMPND", "SOURCE", "EXPDTA", "KEYWDS",
        "REVDAT", "JRNL  ", "HELIX ", "SHEET ", "TURN  "
    };
    long i, num, nlen;
    FILE *fpp;

    num = sizeof titles / sizeof titles[0];

    fpp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fpp) != NULL) {
        nlen = upperstr(str);
        if (!strncmp(str, "ATOM", 4) || !strncmp(str, "HETATM", 6) ||
            !strncmp(str, "END", 3))
            break;
        if (nlen >= 6) {  /* at least 6 characters */
            for (i = 0; i < num - 3; i++)
                if (!strncmp(str, titles[i], 6))
                    fprintf(fp, "%s", str);
            if ((!strncmp(str, "HELIX ", 6) || !strncmp(str, "TURN  ", 6)) &&
                (strchr(chain_list, '*') || strchr(chain_list, str[19])))
                fprintf(fp, "%s", str);
            if (!strncmp(str, "SHEET ", 6) &&
                (strchr(chain_list, '*') || strchr(chain_list, str[21])))
                fprintf(fp, "%s", str);
        }
    }
    close_file(fpp);
}

long number_of_atoms(char *pdbfile, long hetatm, char *ALT_LIST)
/* number of atom records in a PDB file */
{
    char str[BUF512], *pchar;
    long n = 0, nlen;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        if ((pchar = strrchr(str, '\n')) != NULL)
            str[pchar - str] = '\0';
        nlen = upperstr(str);
        if (!strcmp(str, "END"))
            break;
        if (nlen >= 54 && (!strncmp(str, "ATOM", 4)
                           || (hetatm && !strncmp(str, "HETATM", 6)))
            && (strchr(ALT_LIST, '*') || strchr(ALT_LIST, str[16])))
            n++;
    }
    close_file(fp);

    return n;
}

long read_pdb(char *pdbfile, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs, long hetatm,
              char *ALT_LIST)
/* read in a PDB file and do some processing
 * Miscs[][NMISC]: H/A, altLoc, iCode, occ./tempFac./segID/element/charge
 *           col#   0     1       2       3-28 [combined together]
 * 29-32, 4 characters for Model Numbers (in PDB MODEL record, col: 11-14)
 * Element symbol: column# 25-26, right justified
 * 80 - 54 = 26 characters from ATOM/HETATM records 26 + 2 = 28
 *      set the 29th column to '\0' for end of the ATOM/HETATM record
 * 30 -- 33 for the 4-character model number; 34 is set to '\0'
 *   so we have NMISC = 34
 */
{
    char rname[4], str[BUF512], temp[BUF512], *pchar;
    long i, n = 0, nlen, k, modelNum = 0;
    FILE *fp;

    fp = open_file(pdbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL) {
        if ((pchar = strrchr(str, '\n')) != NULL)
            str[pchar - str] = '\0';
        nlen = upperstr(str);
        if (!strncmp(str, "MODEL ", 6)  /* model record */
            &&(sscanf(str, "%*s %ld", &k) == 1))
            modelNum = k;
        if (!strcmp(str, "END"))
            break;
        if (nlen >= 54 &&  /* same as in number_of_atoms */
            (!strncmp(str, "ATOM", 4)
             || (hetatm && !strncmp(str, "HETATM", 6)))
            && (strchr(ALT_LIST, '*') || strchr(ALT_LIST, str[16]))) {
            n++;

            if (AtomSNum != NULL) {  /* read in original atom serial number */
                strncpy(temp, str + 6, 5);
                temp[5] = '\0';  /* 5 digits */
                if (sscanf(temp, "%5ld", &AtomSNum[n]) != 1) {
                    fprintf(stderr, "Atom serial number ? %s ?\n", str);
                    AtomSNum[n] = n;  /* sequential number */
                }
            }

            strncpy(AtomName[n], str + 12, 4);
            AtomName[n][4] = '\0';

            /* for very special case 'P   ', 'N9  ' etc. Fei XU */
            if ((AtomName[n][0] != ' ' && !isdigit(AtomName[n][0]))
                && (AtomName[n][1] == ' ' || isdigit(AtomName[n][1]))
                && AtomName[n][3] == ' ') {  /* isupper for 1st char? */
                strncpy(temp, AtomName[n], 3);
                temp[3] = '\0';
                sprintf(AtomName[n], " %s", temp);
            }

            if (AtomName[n][3] == '*')  /* * to ' */
                AtomName[n][3] = '\'';
            if (!strcmp(AtomName[n], " O1'"))  /* O1' to O4' */
                strcpy(AtomName[n], " O4'");
            if (!strcmp(AtomName[n], " OL "))  /* OL to O1P */
                strcpy(AtomName[n], " O1P");
            if (!strcmp(AtomName[n], " OR "))  /* OR to O2P */
                strcpy(AtomName[n], " O2P");
            if (!strcmp(AtomName[n], " C5A"))  /* C5A to C5M */
                strcpy(AtomName[n], " C5M");
            if (!strcmp(AtomName[n], " O5T"))  /* terminal O5' */
                strcpy(AtomName[n], " O5'");
            if (!strcmp(AtomName[n], " O3T"))  /* terminal O3' */
                strcpy(AtomName[n], " O3'");

            strncpy(rname, str + 17, 3);  /* residue name */
            rname[3] = '\0';
            /* delete ending spaces as in "C  " */
            for (i = 1; i <= 2; i++)
                if (rname[2] == ' ') {
                    rname[2] = rname[1];
                    rname[1] = rname[0];
                    rname[0] = ' ';
                }
            if (rname[2] == ' ')
                fatal("==> residue name [%s] field empty <==\n", str);
            strcpy(ResName[n], rname);
            ChainID[n] = str[21];
            strncpy(temp, str + 22, 4);
            temp[4] = '\0';
            if (sscanf(temp, "%4ld", &ResSeq[n]) != 1) {
                fprintf(stderr, "residue #? ==> %.54s\n", str);
                ResSeq[n] = 9999;
            }
            strncpy(temp, str + 30, 8);
            temp[8] = '\0';
            if (sscanf(temp, "%8lf", &xyz[n][1]) != 1)
                fatal("error reading x-coordinate\n");
            strncpy(temp, str + 38, 8);
            temp[8] = '\0';
            if (sscanf(temp, "%8lf", &xyz[n][2]) != 1)
                fatal("error reading y-coordinate\n");
            strncpy(temp, str + 46, 8);
            temp[8] = '\0';
            if (sscanf(temp, "%8lf", &xyz[n][3]) != 1)
                fatal("error reading z-coordinate\n");
            if (Miscs != NULL) {
                Miscs[n][0] = str[0];  /* H for HETATM, A for ATOM */
                Miscs[n][1] = str[16];  /* alternative location indicator */
                Miscs[n][2] = str[26];  /* code of insertion residues */
                strncpy(Miscs[n] + 3, str + 54, 26);  /* upto the 80th column */
                Miscs[n][29] = '\0';
                sprintf(Miscs[n] + 30, "%4.4ld", modelNum);
                Miscs[n][NMISC] = '\0';  /* just to make sure */
            }
        }
    }
    close_file(fp);

    if (!n)
        fprintf(stderr, "PDB file <%s> has NO ATOM/HETATM records\n", pdbfile);

    return n;
}

void free_pdb(long num, long *AtomSNum, char **AtomName, char **ResName, char *ChainID,
              long *ResSeq, double **xyz, char **Miscs)
/* free all PDB relevant arrays to make code clean and short */
{
    if (AtomSNum != NULL)
        free_lvector(AtomSNum, 1, num);
    if (AtomName != NULL)
        free_cmatrix(AtomName, 1, num, 0, 4);
    if (ResName != NULL)
        free_cmatrix(ResName, 1, num, 0, 3);
    if (ChainID != NULL)
        free_cvector(ChainID, 1, num);
    if (ResSeq != NULL)
        free_lvector(ResSeq, 1, num);
    if (xyz != NULL)
        free_dmatrix(xyz, 1, num, 1, 3);
    if (Miscs != NULL)
        free_cmatrix(Miscs, 1, num, 0, NMISC);
}

void reset_xyz(long num, double **xyz, char *fmt)
/* reset xyz coordinates for PDB and ALCHEMY formats */
{
    double ave_xyz[4], max_xyz[4], min_xyz[4];

    /* check the range of the coordinates */
    max_dmatrix(xyz, num, 3, max_xyz);
    min_dmatrix(xyz, num, 3, min_xyz);
    ave_dmatrix(xyz, num, 3, ave_xyz);

    if (max_dvector(max_xyz, 3) > 9999.99) {
        fprintf(stderr, "xyz coordinate over %s limit. reset origin"
                " to geometrical center\n", fmt);
        move_position(xyz, num, 3, ave_xyz);
    } else if (min_dvector(min_xyz, 3) < -999.99) {
        fprintf(stderr, "xyz coordinate under %s limit. reset origin"
                " to minimum xyz coordinates\n", fmt);
        move_position(xyz, num, 3, min_xyz);
    }
}

void pdb_record(long ib, long ie, long *inum, long idx, char **AtomName, char **ResName,
                char *ChainID, long *ResSeq, double **xyz, char **Miscs, FILE * fp)
/* write out ATOM and HETATM record: xyz could be 1 to [ie - ib + 1] */
{
    char str[BUF512];
    long i, j;

    for (i = ib; i <= ie; i++) {
        (Miscs == NULL) ? strcpy(str, "A  ") : strcpy(str, Miscs[i]);
        j = (idx) ? i - ib + 1 : i;
        fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n", (str[0] == 'A')
                ? "ATOM  " : "HETATM", ++*inum, AtomName[i], str[1],
                ResName[i], ChainID[i], ResSeq[i], str[2], xyz[j][1],
                xyz[j][2], xyz[j][3], str + 3);
    }
}

void write_pdb(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               double **xyz, char **Miscs, char *pdbfile)
{
    long inum = 0;
    FILE *fp;

    reset_xyz(num, xyz, "f8.3");

    fp = open_file(pdbfile, "w");
    fprintf(fp, "REMARK    %s\n", X3DNA_VER);
    pdb_record(1, num, &inum, 0, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, fp);
    fprintf(fp, "END\n");
    close_file(fp);
}

void write_pdbcnt(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                  double **xyz, long **connect, char *pdbfile)
{
    long i, inum = 0, j;
    FILE *fp;

    reset_xyz(num, xyz, "f8.3");

    fp = open_file(pdbfile, "w");
    fprintf(fp, "REMARK    %s\n", X3DNA_VER);
    pdb_record(1, num, &inum, 0, AtomName, ResName, ChainID, ResSeq, xyz, NULL, fp);

    for (i = 1; i <= num; i++)
        if (connect[i][7]) {
            fprintf(fp, "CONECT%5ld", i);
            for (j = 1; j <= connect[i][7]; j++)
                fprintf(fp, "%5ld", connect[i][j]);
            /* March-2-2005: Guanming.Wang.Adv04@Alum.Dartmouth.ORG */
            for (j = 6 + 5 * (connect[i][7] + 1); j <= 70; j++)
                fprintf(fp, " ");
            fprintf(fp, "\n");
        }
    fprintf(fp, "END\n");

    close_file(fp);
}

void move_position(double **d, long nr, long nc, double *mpos)
{
    long i, j;

    for (i = 1; i <= nr; i++)
        for (j = 1; j <= nc; j++)
            d[i][j] -= mpos[j];
}

long **residue_idx(long num, long *ResSeq, char **Miscs, char *ChainID, char **ResName,
                   long *num_residue)
/* number of residues, and starting-ending indexes for each */
{
    char iCode;
    char **bidx;
    long i, n, **seidx, *temp;

    bidx = cmatrix(1, num, 0, 12);  /* normally 9 */
    temp = lvector(1, num);

    for (i = 1; i <= num; i++) {
        iCode = (Miscs == NULL) ? ' ' : Miscs[i][2];
        sprintf(bidx[i], "%3s%c%4ld%c", ResName[i], ChainID[i], ResSeq[i], iCode);
    }
    for (i = 1; i < num; i++)
        temp[i] = strcmp(bidx[i + 1], bidx[i]) ? 1 : 0;
    temp[num] = 1;

    n = 0;  /* get number of residues */
    for (i = 1; i <= num; i++)
        if (temp[i])
            ++n;

    seidx = lmatrix(1, n, 1, 2);  /* allocate spaces */
    n = 0;
    for (i = 1; i <= num; i++)
        if (temp[i])
            seidx[++n][2] = i;
    for (i = 2; i <= n; i++)
        seidx[i][1] = seidx[i - 1][2] + 1;
    seidx[1][1] = 1;

    *num_residue = n;

    free_cmatrix(bidx, 1, num, 0, 12);
    free_lvector(temp, 1, num);

    return seidx;
}

long residue_ident(char **AtomName, double **xyz, char **Miscs, long ib, long ie)
/* identifying a residue as follows:
 *  R-base  Y-base  amino-acid, others [default]
 *   +1        0        -1        -2 [default]
 */
{
    double d1, d2, d3, dcrt = 2.0, dcrt2 = 3.0;
    long id = -2;
    long CA, C, N1, C2, C6, N9;

    N9 = find_1st_atom(" N9 ", AtomName, ib, ie, "");
    N1 = find_1st_atom(" N1 ", AtomName, ib, ie, "");
    C2 = find_1st_atom(" C2 ", AtomName, ib, ie, "");
    C6 = find_1st_atom(" C6 ", AtomName, ib, ie, "");
    if (N1 && C2 && C6) {
        d1 = p1p2_dist(xyz[C2], xyz[N1]);
        d2 = p1p2_dist(xyz[C6], xyz[N1]);
        d3 = p1p2_dist(xyz[C6], xyz[C2]);
        if (d1 <= dcrt && d2 <= dcrt && d3 <= dcrt2) {
            id = 0;
            if (N9) {
                d3 = p1p2_dist(xyz[N9], xyz[N1]);
                if (dval_in_range(d3, 3.15, 4.85))  /* ~4.0 */
                    id = 1;
            }
        }
        return id;
    }

    CA = find_1st_atom(" CA ", AtomName, ib, ie, "");
    C = find_1st_atom(" C  ", AtomName, ib, ie, "");
    if (!C)  /* if C does not exist, use N */
        C = find_1st_atom(" N  ", AtomName, ib, ie, "");
    if (CA && C && Miscs[CA][0] == 'A' && Miscs[C][0] == 'A' &&  /* ATOM record */
        within_limits(xyz[CA], xyz[C], 0, dcrt))
        id = -1;

    return id;  /* other cases */
}

char base_ident(char *rname, char *idmsg, char **baselist, long num_sb)
/* get one letter base name according to 3-letter to 1-letter conversion */
{
    long i;

    for (i = 1; i <= num_sb; i++)
        if (!strncmp(rname, baselist[i], 3))
            break;
    if (i > num_sb)
        fatal("unknown base [%s]: add it to file <baselist.dat>\n", idmsg);
    if (!isupper((int) baselist[i][3]) || baselist[i][3] == 'P')
        fprintf(stderr, "uncommon %s assigned to: %c\n", idmsg, baselist[i][3]);
    return baselist[i][3];
}

void get_baselist(char **baselist, long *num_sb)
/* get the correspondence between 3-letter and 1-letter base names */
{
    char BDIR[BUF512], str[BUF512], base3[BUF512], base1[BUF512];
    long n = 0;
    FILE *fp;

    get_BDIR(BDIR, BASE_FILE, 1);
    strcat(BDIR, BASE_FILE);
    fprintf(stderr, " ...... reading file: %s ...... \n", BASE_FILE);

    fp = open_file(BDIR, "r");
    while ((fgets(str, sizeof str, fp) != NULL)) {
        if (str[0] == '#')
            continue;
        if (sscanf(str, "%s %s", base3, base1) != 2)
            continue;
        if (base3[0] == '#' || base1[0] == '#')
            continue;
        if (strlen(base3) > 3 || strlen(base1) != 1) {
            fprintf(stderr, "ignoring unacceptable format: %s\n", str);
            continue;
        }
        if (++n > BUF512)
            fatal("too many base types\n");
        sprintf(baselist[n], "%3.3s%c", base3, base1[0]);
    }
    *num_sb = n;

    close_file(fp);
}

void get_seq(long num_residue, long **seidx, char **AtomName, char **ResName,
             char *ChainID, long *ResSeq, char **Miscs, double **xyz, char *bseq,
             long *RY)
/* get base sequence of DNA. similar to GET_BPSEQ */
{
    char idmsg[BUF512], **baselist;
    long i, ib, ie, num_sb;

    baselist = cmatrix(1, BUF512, 0, 4);
    get_baselist(baselist, &num_sb);

    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];
        RY[i] = residue_ident(AtomName, xyz, Miscs, ib, ie);
        if (RY[i] >= 0) {
            sprintf(idmsg, "residue %3s %4ld%c on chain %c [#%ld]",
                    ResName[ib], ResSeq[ib], Miscs[ib][2], ChainID[ib], i);
            bseq[i] = base_ident(ResName[ib], idmsg, baselist, num_sb);
        }
    }
    free_cmatrix(baselist, 1, BUF512, 0, 4);
}

void get_bpseq(long ds, long num_bp, long **pair_num, long **seidx, char **AtomName,
               char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
               char **bp_seq, long *RY)
/* get base (pair) sequence, change (+) modified residue to lower case */
{
    char idmsg[BUF512], **baselist;
    long i, ib, ie, j, num_sb, rnum;

    baselist = cmatrix(1, BUF512, 0, 4);
    get_baselist(baselist, &num_sb);

    for (i = 1; i <= ds; i++) {
        for (j = 1; j <= num_bp; j++) {
            rnum = pair_num[i][j];
            ib = seidx[rnum][1];
            ie = seidx[rnum][2];
            RY[rnum] = residue_ident(AtomName, xyz, Miscs, ib, ie);
            sprintf(idmsg, "residue %3s %4ld%c on chain %c [#%ld]",
                    ResName[ib], ResSeq[ib], Miscs[ib][2], ChainID[ib], rnum);
            if (RY[rnum] >= 0)
                bp_seq[i][j] = base_ident(ResName[ib], idmsg, baselist, num_sb);
            else
                fatal("Non-base: %s\n", idmsg);
        }
    }
    free_cmatrix(baselist, 1, BUF512, 0, 4);
}

long num_strmatch(char *str, char **strmat, long nb, long ne)
/*  return number of matchs of str in strmat */
{
    long i, num = 0;

    for (i = nb; i <= ne; i++)
        if (!strcmp(str, strmat[i]))
            num++;

    return num;
}

long find_1st_atom(char *str, char **strmat, long nb, long ne, char *idmsg)
/* return index of the first match, or 0 for no-match */
{
    long i, num;

    num = num_strmatch(str, strmat, nb, ne);

    if (!num) {
        if (strcmp(idmsg, ""))
            fprintf(stderr, "missing \"%s\" atom %s\n", str, idmsg);
        return 0;
    }
    if (num > 1 && strcmp(idmsg, "")) {
        fprintf(stderr, "more than one %s atoms %s\n", str, idmsg);
        fprintf(stderr, "   *****the first atom is used*****\n");
    }
    for (i = nb; i <= ne; i++)
        if (!strcmp(str, strmat[i]))
            break;
    return i;
}

double torsion(double **d)
/* get torsion angle a-b-c-d in degrees */
{
    double ang_deg, **vec3;
    long i;

    vec3 = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= 3; i++) {
        ddxyz(d[i], d[i + 1], vec3[i]);
        if (veclen(vec3[i]) > BOND_UPPER_LIMIT) {
            ang_deg = EMPTY_NUMBER;  /* not directly linked */
            goto RTN_HERE;
        }
    }
    negate_xyz(vec3[1]);  /* b-->a */
    ang_deg = vec_ang(vec3[1], vec3[3], vec3[2]);

  RTN_HERE:
    free_dmatrix(vec3, 1, 3, 1, 3);
    return ang_deg;
}

double torsion2(double **d)
/* get torsion angle a-b-c-d in degrees: no consideration of breaks */
{
    double ang_deg, **vec3;
    long i;

    vec3 = dmatrix(1, 3, 1, 3);
    for (i = 1; i <= 3; i++)
        ddxyz(d[i], d[i + 1], vec3[i]);
    negate_xyz(vec3[1]);  /* b-->a */
    ang_deg = vec_ang(vec3[1], vec3[3], vec3[2]);
    free_dmatrix(vec3, 1, 3, 1, 3);

    return ang_deg;
}

void get_BDIR(char *BDIR, char *filename, long prt_msg)
/* search the directory containing standard bases & parameter files:
   (1) current directory
   (2) directory defined by the environmental variable "X3DNA"
   (3) standard directory $HOME/X3DNA/BASEPARS/ */
{
    char *temp;
    long iscd = 0;
    FILE *fp;

    fp = fopen(filename, "r");  /* check current directory */
    if (fp != NULL)
        iscd = 1;
    close_file(fp);

    if (iscd)
        strcpy(BDIR, "./");  /* current directory */
    else if ((temp = getenv("X3DNA")) != NULL) {
        strcpy(BDIR, temp);
        check_slash(BDIR);
        strcat(BDIR, "BASEPARS/");
    } else if ((temp = getenv("HOME")) != NULL ||  /* Unix */
               (temp = getenv("HOMEDRIVE")) != NULL) {  /* PC */
        strcpy(BDIR, temp);
        check_slash(BDIR);
        strcat(BDIR, "X3DNA/BASEPARS/");
    } else
        fatal("cannot locate base geometry and parameter files\n");
    if (prt_msg)
        fprintf(stderr, "\n ...... %s ...... \n", BDIR);
}

void align2zaxis(long num, double *haxis, double **rotmat, double **xyz, double **xyzH)
/* transform xyz coordinates so that "haxis" is aligned to the z-axis */
{
    double z[4] = { EMPTY_NUMBER, 0.0, 0.0, 1.0 };
    double ang_deg, hinge[4], **rotmatT;

    rotmatT = dmatrix(1, 3, 1, 3);

    cross(haxis, z, hinge);
    ang_deg = magang(haxis, z);
    arb_rotation(hinge, ang_deg, rotmat);
    transpose_matrix(rotmat, 3, 3, rotmatT);
    multi_matrix(xyz, num, 3, rotmatT, 3, 3, xyzH);

    free_dmatrix(rotmatT, 1, 3, 1, 3);
}

void cov_matrix(double **a, double **b, long nr, long nc, double **cmtx)
/* calculate the covariance matrix between two matrices */
{
    double ave_a[4], ave_b[4];
    double **ta, **ta_x_b;
    long i, j;

    ave_dmatrix(a, nr, nc, ave_a);
    ave_dmatrix(b, nr, nc, ave_b);

    ta = dmatrix(1, nc, 1, nr);  /* transpose of a */
    ta_x_b = dmatrix(1, nc, 1, nc);  /* transpose-a multiply b */

    transpose_matrix(a, nr, nc, ta);
    multi_matrix(ta, nc, nr, b, nr, nc, ta_x_b);

    for (i = 1; i <= nc; i++)
        for (j = 1; j <= nc; j++)
            cmtx[i][j] = (ta_x_b[i][j] - ave_a[i] * ave_b[j] * nr) / (nr - 1);

    free_dmatrix(ta, 1, nc, 1, nr);
    free_dmatrix(ta_x_b, 1, nc, 1, nc);
}

double ls_fitting(double **sxyz, double **exyz, long n, double **fitted_xyz, double **R,
                  double *orgi)
/* least-squares fitting between two structures */
{
    double temp, rms_value;
    double ave_exyz[4], ave_sxyz[4], D[5];
    double **N, **U, **V;
    long i, j;

    if (n < 3)
        fatal("too few atoms for least-squares fitting\n");

    /* get the covariance matrix U */
    U = dmatrix(1, 3, 1, 3);
    cov_matrix(sxyz, exyz, n, 3, U);

    /* get 4-by-4 symmetric matrix N */
    N = dmatrix(1, 4, 1, 4);
    N[1][1] = U[1][1] + U[2][2] + U[3][3];
    N[2][2] = U[1][1] - U[2][2] - U[3][3];
    N[3][3] = -U[1][1] + U[2][2] - U[3][3];
    N[4][4] = -U[1][1] - U[2][2] + U[3][3];
    N[1][2] = U[2][3] - U[3][2];
    N[2][1] = N[1][2];
    N[1][3] = U[3][1] - U[1][3];
    N[3][1] = N[1][3];
    N[1][4] = U[1][2] - U[2][1];
    N[4][1] = N[1][4];
    N[2][3] = U[1][2] + U[2][1];
    N[3][2] = N[2][3];
    N[2][4] = U[3][1] + U[1][3];
    N[4][2] = N[2][4];
    N[3][4] = U[2][3] + U[3][2];
    N[4][3] = N[3][4];

    /* get N's eigenvalues and eigenvectors */
    V = dmatrix(1, 4, 1, 4);
    jacobi(N, 4, D, V);

    /* get the rotation matrix */
    for (i = 1; i <= 4; i++)
        for (j = 1; j <= 4; j++)
            N[i][j] = V[i][4] * V[j][4];
    R[1][1] = N[1][1] + N[2][2] - N[3][3] - N[4][4];
    R[1][2] = 2 * (N[2][3] - N[1][4]);
    R[1][3] = 2 * (N[2][4] + N[1][3]);
    R[2][1] = 2 * (N[3][2] + N[1][4]);
    R[2][2] = N[1][1] - N[2][2] + N[3][3] - N[4][4];
    R[2][3] = 2 * (N[3][4] - N[1][2]);
    R[3][1] = 2 * (N[4][2] - N[1][3]);
    R[3][2] = 2 * (N[4][3] + N[1][2]);
    R[3][3] = N[1][1] - N[2][2] - N[3][3] + N[4][4];

    ave_dmatrix(sxyz, n, 3, ave_sxyz);
    ave_dmatrix(exyz, n, 3, ave_exyz);

    /* fitted sxyz origin */
    for (i = 1; i <= 3; i++)
        orgi[i] = ave_exyz[i] - dot(ave_sxyz, R[i]);

    /* fitted sxyz coordinates */
    for (i = 1; i <= n; i++)
        for (j = 1; j <= 3; j++)
            fitted_xyz[i][j] = dot(sxyz[i], R[j]) + orgi[j];

    /* rms deviation */
    temp = 0.0;
    for (i = 1; i <= n; i++) {
        ddxyz(fitted_xyz[i], exyz[i], D);
        temp += dot(D, D);
    }
    rms_value = sqrt(temp / n);

    free_dmatrix(U, 1, 3, 1, 3);
    free_dmatrix(N, 1, 4, 1, 4);
    free_dmatrix(V, 1, 4, 1, 4);

    return rms_value;
}

void ls_plane(double **bxyz, long n, double *pnormal, double *ppos, double *odist,
              double *adist)
/* fit a plane to a set of points by least squares */
{
    double D[4];
    double **cov_mtx, **identityV, **V;
    long i, j, nml = 0;

    if (n < 3)
        fatal("too few atoms for least-squares fitting\n");

    cov_mtx = dmatrix(1, 3, 1, 3);
    V = dmatrix(1, 3, 1, 3);
    identityV = dmatrix(1, 3, 1, 3);

    cov_matrix(bxyz, bxyz, n, 3, cov_mtx);
    jacobi(cov_mtx, 3, D, V);

    identity_matrix(identityV, 3);
    for (i = 1; i <= 3 && !nml; i++)
        for (j = 1; j <= 3 && !nml; j++)
            if (fabs(V[i][j] - identityV[i][j]) > XEPS)
                nml = 1;
    if (nml)
        for (i = 1; i <= 3; i++)
            pnormal[i] = V[i][1];
    else {  /* V is an identity matrix */
        pnormal[1] = 0.0;
        pnormal[2] = 0.0;
        pnormal[3] = 1.0;
    }

    ave_dmatrix(bxyz, n, 3, ppos);

    /* make the z-component of pnormal to be positive */
    if (pnormal[3] < 0)
        negate_xyz(pnormal);

    /* distance from the origin to the plane */
    *odist = dot(ppos, pnormal);

    /* distance from each point to the plane */
    for (i = 1; i <= n; i++)
        adist[i] = dot(bxyz[i], pnormal) - *odist;

    free_dmatrix(cov_mtx, 1, 3, 1, 3);
    free_dmatrix(V, 1, 3, 1, 3);
    free_dmatrix(identityV, 1, 3, 1, 3);
}

void arb_rotation(double *va, double ang_deg, double **rot_mtx)
/* get the arbitrary rotation matrix */
{
    double c, dc, s, vlen;
    long i;

    vlen = veclen(va);
    if (vlen < XEPS)  /* [0 0 0] */
        identity_matrix(rot_mtx, 3);
    else {
        for (i = 1; i <= 3; i++)
            va[i] /= vlen;  /* unit vector */
        ang_deg = deg2rad(ang_deg);
        c = cos(ang_deg);
        s = sin(ang_deg);
        dc = 1 - c;
        rot_mtx[1][1] = c + dc * va[1] * va[1];
        rot_mtx[1][2] = va[1] * va[2] * dc - va[3] * s;
        rot_mtx[1][3] = va[1] * va[3] * dc + va[2] * s;
        rot_mtx[2][1] = va[1] * va[2] * dc + va[3] * s;
        rot_mtx[2][2] = c + dc * va[2] * va[2];
        rot_mtx[2][3] = va[2] * va[3] * dc - va[1] * s;
        rot_mtx[3][1] = va[1] * va[3] * dc - va[2] * s;
        rot_mtx[3][2] = va[2] * va[3] * dc + va[1] * s;
        rot_mtx[3][3] = c + dc * va[3] * va[3];
    }
}

double vec_ang(double *va, double *vb, double *vref)
/* angle in degrees between va and vb with vref for sign control
   va & vb are unchanged by making an addition copy of each
   all three vectors are 1-by-3 */
{
    double ang_deg, vc[4], va_cp[4], vb_cp[4];

    /* make a copy of va and vb */
    cpxyz(va, va_cp);
    cpxyz(vb, vb_cp);

    /* get orthogonal components */
    vec_orth(va_cp, vref);
    vec_orth(vb_cp, vref);

    /* angle in absolute sense */
    ang_deg = magang(va_cp, vb_cp);

    /* sign control */
    cross(va_cp, vb_cp, vc);
    if (dot(vc, vref) < 0)
        ang_deg = -ang_deg;

    return ang_deg;
}

void get_vector(double *va, double *vref, double deg_ang, double *vo)
/* get the vector which has certain angle with another vector */
{
    double va_cp[4], **temp;
    long i;

    cpxyz(va, va_cp);  /* make a copy of va: <ana_fncs.c> */
    if (dot(va_cp, vref) > XEPS) {
        fprintf(stderr, "Angle between va/vref: %.3f degrees\n", magang(va_cp, vref));
        vec_orth(va_cp, vref);
    }

    temp = dmatrix(1, 3, 1, 3);

    arb_rotation(vref, deg_ang, temp);
    for (i = 1; i <= 3; i++)
        vo[i] = dot(temp[i], va_cp);
    vec_norm(vo);

    free_dmatrix(temp, 1, 3, 1, 3);
}

void rotate(double **a, long i, long j, long k, long l, double *g, double *h, double s,
            double tau)
{
    *g = a[i][j];
    *h = a[k][l];
    a[i][j] = *g - s * (*h + *g * tau);
    a[k][l] = *h + s * (*g - *h * tau);
}

void eigsrt(double *d, double **v, long n)
/* sort eigenvalues into ascending order and rearrange eigenvectors */
{
    double p;
    long i, j, k;

    for (i = 1; i < n; i++) {
        p = d[k = i];
        for (j = i + 1; j <= n; j++)
            if (d[j] < p)
                p = d[k = j];
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 1; j <= n; j++)
                dval_swap(&v[j][i], &v[j][k]);
        }
    }
}

void jacobi(double **a, long n, double *d, double **v)
{
    long i, j, iq, ip;
    double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

    b = dvector(1, n);
    z = dvector(1, n);
    identity_matrix(v, n);
    for (ip = 1; ip <= n; ip++) {
        b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
    }
    for (i = 1; i <= 100; i++) {
        sm = 0.0;
        for (ip = 1; ip <= n - 1; ip++) {
            for (iq = ip + 1; iq <= n; iq++)
                sm += fabs(a[ip][iq]);
        }
        if (sm < XEPS) {
            free_dvector(z, 1, n);
            free_dvector(b, 1, n);
            eigsrt(d, v, n);
            return;
        }
        if (i < 4)
            tresh = 0.2 * sm / (n * n);
        else
            tresh = 0.0;
        for (ip = 1; ip <= n - 1; ip++) {
            for (iq = ip + 1; iq <= n; iq++) {
                g = 100.0 * fabs(a[ip][iq]);
                if (i > 4 && (fabs(d[ip]) + g) == fabs(d[ip])
                    && (fabs(d[iq]) + g) == fabs(d[iq]))
                    a[ip][iq] = 0.0;
                else if (fabs(a[ip][iq]) > tresh) {
                    h = d[iq] - d[ip];
                    if ((fabs(h) + g) == fabs(h))
                        t = a[ip][iq] / h;
                    else {
                        theta = 0.5 * h / a[ip][iq];
                        t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0)
                            t = -t;
                    }
                    c = 1.0 / sqrt(1 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 1; j <= ip - 1; j++)
                        rotate(a, j, ip, j, iq, &g, &h, s, tau);
                    for (j = ip + 1; j <= iq - 1; j++)
                        rotate(a, ip, j, j, iq, &g, &h, s, tau);
                    for (j = iq + 1; j <= n; j++)
                        rotate(a, ip, j, iq, j, &g, &h, s, tau);
                    for (j = 1; j <= n; j++)
                        rotate(v, j, ip, j, iq, &g, &h, s, tau);
                }
            }
        }
        for (ip = 1; ip <= n; ip++) {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    fatal("too many iterations\n");
}

void dludcmp(double **a, long n, long *indx, double *d)
{
    double big, dum, sum, temp;
    double *vv;
    long i, j, k;
    long imax = 0;  /* initialization */

    vv = dvector(1, n);
    *d = 1.0;
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++)
            if ((temp = fabs(a[i][j])) > big)
                big = temp;
        if (big == 0.0)
            fatal("singular matrix in routine dludcmp\n");
        vv[i] = 1.0 / big;
    }
    for (j = 1; j <= n; j++) {
        for (i = 1; i < j; i++) {
            sum = a[i][j];
            for (k = 1; k < i; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 1; k <= n; k++)
                dval_swap(&a[imax][k], &a[j][k]);
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j][j] == 0.0)
            a[j][j] = XEPS;
        if (j != n) {
            dum = 1.0 / a[j][j];
            for (i = j + 1; i <= n; i++)
                a[i][j] *= dum;
        }
    }
    free_dvector(vv, 1, n);
}

void dlubksb(double **a, long n, long *indx, double *b)
{
    double sum;
    long i, ii = 0, ip, j;

    for (i = 1; i <= n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii; j <= i - 1; j++)
                sum -= a[i][j] * b[j];
        else if (sum)
            ii = i;
        b[i] = sum;
    }
    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i + 1; j <= n; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

void dinverse(double **a, long n, double **y)
{
    double d, *col;
    long i, j, *indx;

    col = dvector(1, n);
    indx = lvector(1, n);

    dludcmp(a, n, indx, &d);

    for (j = 1; j <= n; j++) {
        for (i = 1; i <= n; i++)
            col[i] = 0.0;
        col[j] = 1.0;
        dlubksb(a, n, indx, col);
        for (i = 1; i <= n; i++)
            y[i][j] = col[i];
    }

    free_lvector(indx, 1, n);
    free_dvector(col, 1, n);
}

void rotx(double ang_deg, double **rotmat)
{
    double ang, c, s;

    ang = deg2rad(ang_deg);
    c = cos(ang);
    s = sin(ang);

    rotmat[1][1] = 1.0;
    rotmat[1][2] = 0.0;
    rotmat[1][3] = 0.0;
    rotmat[2][1] = 0.0;
    rotmat[2][2] = c;
    rotmat[2][3] = -s;
    rotmat[3][1] = 0.0;
    rotmat[3][2] = s;
    rotmat[3][3] = c;
}

void roty(double ang_deg, double **rotmat)
{
    double ang, c, s;

    ang = deg2rad(ang_deg);
    c = cos(ang);
    s = sin(ang);

    rotmat[1][1] = c;
    rotmat[1][2] = 0.0;
    rotmat[1][3] = s;
    rotmat[2][1] = 0.0;
    rotmat[2][2] = 1.0;
    rotmat[2][3] = 0.0;
    rotmat[3][1] = -s;
    rotmat[3][2] = 0.0;
    rotmat[3][3] = c;
}

void rotz(double ang_deg, double **rotmat)
{
    double ang, c, s;

    ang = deg2rad(ang_deg);
    c = cos(ang);
    s = sin(ang);

    rotmat[1][1] = c;
    rotmat[1][2] = -s;
    rotmat[1][3] = 0.0;
    rotmat[2][1] = s;
    rotmat[2][2] = c;
    rotmat[2][3] = 0.0;
    rotmat[3][1] = 0.0;
    rotmat[3][2] = 0.0;
    rotmat[3][3] = 1.0;
}

void get_alc_nums(char *alcname, long *num, long *nbond)
{
    char str[BUF512];
    FILE *fp;

    fp = open_file(alcname, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld %*s %ld", num, nbond) != 2)
        fatal("can not read atom & bond numbers: %s\n", str);

    if (!*num)
        fprintf(stderr, "ALCHEMY file <%s> has NO atoms\n", alcname);

    close_file(fp);
}

void read_alc(char *alcname, long *num, long *nbond, char **AtomName, double **xyz,
              long *ibase, long **linkage)
{
    char str[BUF512], temp[BUF512];
    double c;
    long i;
    FILE *fp;

    fp = open_file(alcname, "r");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld %*s %ld", num, nbond) != 2)
        fatal("can not read atom & bond numbers: %s\n", str);

    for (i = 1; i <= *num; i++)
        if (fgets(str, sizeof str, fp) != NULL) {
            strncpy(AtomName[i], str + 6, 2);
            AtomName[i][2] = '\0';
            upperstr(AtomName[i]);  /* to upper case */
            strncpy(temp, str + 11, 9);
            temp[9] = '\0';
            if (sscanf(temp, "%lf", &xyz[i][1]) != 1)
                fatal("error reading x-coordinate\n");
            strncpy(temp, str + 20, 9);
            temp[9] = '\0';
            if (sscanf(temp, "%lf", &xyz[i][2]) != 1)
                fatal("error reading y-coordinate\n");
            strncpy(temp, str + 29, 9);
            temp[9] = '\0';
            if (sscanf(temp, "%lf", &xyz[i][3]) != 1)
                fatal("error reading z-coordinate\n");
            strncpy(temp, str + 40, 9);
            temp[9] = '\0';
            if (sscanf(temp, "%lf", &c) == 1)
                ibase[i] = lround(10.0 * c);
            else
                ibase[i] = NON_WC_IDX;  /* default */
        } else
            fatal("error in reading atom records\n");

    for (i = 1; i <= *nbond; i++)
        if (fgets(str, sizeof str, fp) != NULL) {
            strncpy(temp, str + 6, 5);
            temp[5] = '\0';
            if (sscanf(temp, "%ld", &linkage[i][1]) != 1)
                fatal("error reading linkage atom 1\n");
            strncpy(temp, str + 12, 5);
            temp[5] = '\0';
            if (sscanf(temp, "%ld", &linkage[i][2]) != 1)
                fatal("error reading linkage atom 2\n");
        } else
            fatal("error in reading linkage information\n");

    close_file(fp);
}

void write_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
               long **linkage, char *alcfile)
{
    char aname[5];
    long i;
    FILE *fp;

    reset_xyz(num, xyz, "f9.4");

    fp = open_file(alcfile, "w");
    fprintf(fp, "%5ld ATOMS, %5ld BONDS\n", num, nbond);
    for (i = 1; i <= num; i++) {
        if (sscanf(AtomName[i], "%s", aname) == EOF)
            strcpy(aname, "C ");
        if (isdigit((int) aname[0])) {  /* like "5H1 " etc */
            aname[0] = aname[1];
            aname[1] = aname[2];
        }
        if (isdigit((int) aname[1]) || aname[1] == '\0')
            aname[1] = ' ';
        aname[2] = '\0';
        fprintf(fp, "%5ld %-2s   %9.4f%9.4f%9.4f  %9.4f\n", i, aname, xyz[i][1],
                xyz[i][2], xyz[i][3], (ibase == NULL) ? 0.0 : ibase[i] / 10.0);
    }
    for (i = 1; i <= nbond; i++)
        fprintf(fp, "%5ld %5ld %5ld\n", i, linkage[i][1], linkage[i][2]);
    close_file(fp);
}

void free_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
              long zero_1, long **linkage)
/* free all ALCHEMY relevant arrays to make code clean and short */
{
    free_cmatrix(AtomName, 1, num, 0, 2);
    free_dmatrix(xyz, 1, num, 1, 3);
    free_lvector(ibase, 1, num);
    free_lmatrix(linkage, 1, nbond, zero_1, 2);
}

void cnct_org(long num_bp, long ia, long ib, char **tAtomName, double **txyz,
              long *tibase, long **tlinkage, double **org_xyz)
/* add lines connecting consecutive origins to ALCHEMY file */
{
    long i, ik;

    for (i = 1; i <= num_bp; i++) {
        ik = ia + i;
        strcpy(tAtomName[ik], "O");
        cpxyz(org_xyz[i], txyz[ik]);
        tibase[ik] = NON_WC_IDX;  /* non-common base atoms */
    }

    for (i = 1; i < num_bp; i++) {
        ik = ib + i;
        tlinkage[ik][1] = ia + i;
        tlinkage[ik][2] = ia + i + 1;
    }
}

long get_num_bp(char *inpfile)
/* get number of base-pairs */
{
    char str[BUF512];
    long num_bp = 0;
    FILE *fp;

    fp = open_file(inpfile, "r");
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", &num_bp) != 1)
        fatal("error reading parameter file (number of base-pairs): %s\n", str);
    if (num_bp <= 0)
        fatal("number of base-pairs [%ld] <= 0\n", num_bp);
    close_file(fp);

    return num_bp;
}

void dsort(long n, double *a, long *idx)
/* sort a double vector into ascending order & keep the index
   use Shell's method as in NR in C book Ed. 2, pp.331-332
   cf. lsort below */
{
    double v;
    long i, inc, iv, j;

    inc = 1;
    do {
        inc *= 3;
        inc++;
    } while (inc <= n);

    for (i = 1; i <= n; i++)
        idx[i] = i;

    do {
        inc /= 3;
        for (i = inc + 1; i <= n; i++) {
            v = a[i];
            iv = idx[i];
            j = i;
            while (a[j - inc] > v) {
                a[j] = a[j - inc];
                idx[j] = idx[j - inc];
                j -= inc;
                if (j <= inc)
                    break;
            }
            a[j] = v;
            idx[j] = iv;
        }
    } while (inc > 1);
}

void lsort(long n, long *a, long *idx)
/* sort a long vector into ascending order & keep the index
   use Shell's method as in NR in C book Ed. 2, pp.331-332 */
{
    long v;
    long i, inc, iv, j;

    inc = 1;
    do {
        inc *= 3;
        inc++;
    } while (inc <= n);

    for (i = 1; i <= n; i++)
        idx[i] = i;

    do {
        inc /= 3;
        for (i = inc + 1; i <= n; i++) {
            v = a[i];
            iv = idx[i];
            j = i;
            while (a[j - inc] > v) {
                a[j] = a[j - inc];
                idx[j] = idx[j - inc];
                j -= inc;
                if (j <= inc)
                    break;
            }
            a[j] = v;
            idx[j] = iv;
        }
    } while (inc > 1);
}

void lreverse(long ia, long n, long *lvec)
/* reverse a long vector: output replaces the original */
{
    long i, *ltmp;

    ltmp = lvector(1, n);
    for (i = 1; i <= n; i++)
        ltmp[i] = lvec[n + ia - i];
    for (i = 1; i <= n; i++)
        lvec[ia + i - 1] = ltmp[i];

    free_lvector(ltmp, 1, n);
}

void fig_title(FILE * fp)
/* write a title to XFIG 3.2 format */
{
    fprintf(fp,
            "#FIG 3.2  Creator: %s\n"
            "Portrait\n"
            "Flush Left\n"
            "Inches\n" "Letter\n" "100.00\n" "Single\n" "-2\n" "1200 2\n\n", X3DNA_VER);
}

void ps_title_cmds(FILE * fp, char *imgfile, long *bbox)
{
    char BDIR[BUF512], str[BUF512];
    char *ps_image_par = "ps_image.par";
    long i;
    time_t run_time;
    FILE *fpp;

    run_time = time(NULL);

    fprintf(fp, "%%!PS-Adobe-3.0\n");
    fprintf(fp, "%%%%Title: (%s)\n", imgfile);
    fprintf(fp, "%%%%Creator: (%s)\n", X3DNA_VER);
    strcpy(str, ctime(&run_time));
    str[strlen(str) - 1] = '\0';  /* change \n to \0 */
    fprintf(fp, "%%%%CreationDate: (%s)\n", str);
    fprintf(fp, "%%%%Orientation: Portrait\n");
    fprintf(fp, "%%%%BoundingBox: ");
    for (i = 1; i <= 4; i++)
        fprintf(fp, "%6ld", bbox[i]);
    fprintf(fp, "\n\n");

    /* read in color parameter file */
    get_BDIR(BDIR, ps_image_par, 1);
    strcat(BDIR, ps_image_par);
    fprintf(stderr, " ...... reading file: %s ...... \n", ps_image_par);

    fpp = open_file(BDIR, "r");
    while (fgets(str, sizeof str, fpp) != NULL)
        fprintf(fp, "%s", str);
    close_file(fpp);
}

void get_fig_xy(long num, double **xyz, long nO, double **oxyz, long *urxy,
                long frame_box, FILE * fp)
/* reset x/y coordinates to XFIG units */
{
    char *format = "%6ld%6ld";
    double paper_size[2] = { 8.5, 11.0 };  /* US letter size */
    double max_xy[3];
    long i, j, llxy[3];

    /* change y-axis according to XFIG coordinate system */
    for (i = 1; i <= num; i++)
        xyz[i][2] = urxy[2] - xyz[i][2];
    if (nO)
        for (i = 1; i <= nO; i++)
            oxyz[i][2] = urxy[2] - oxyz[i][2];

    max_dmatrix(xyz, num, 2, max_xy);
    for (i = 1; i <= 2; i++)
        urxy[i] = lround(max_xy[i]);  /* new urxy */

    /* 1/4 offset the figure on a US letter (8.5in-by-11in) */
    for (i = 1; i <= 2; i++)
        llxy[i] = lround(0.25 * (paper_size[i - 1] * 1200 - urxy[i]));
    for (i = 1; i <= num; i++)
        for (j = 1; j <= 2; j++)
            xyz[i][j] += llxy[j];
    if (nO)
        for (i = 1; i <= nO; i++)
            for (j = 1; j <= 2; j++)
                oxyz[i][j] += llxy[j];

    max_dmatrix(xyz, num, 2, max_xy);
    for (i = 1; i <= 2; i++)
        urxy[i] = lround(max_xy[i]);

    fig_title(fp);

    if (frame_box) {
        /* draw a box around the figure at the deepest depth (999) */
        fprintf(fp, "# draw a boundary box here\n");
        fprintf(fp, "2 3 0 1 0 0 999 0 -1 0.0 2 1 0 0 0 5\n");
        fprintf(fp, format, llxy[1] - FIG_BOUND, llxy[2] - FIG_BOUND);
        fprintf(fp, format, urxy[1] + FIG_BOUND, llxy[2] - FIG_BOUND);
        fprintf(fp, format, urxy[1] + FIG_BOUND, urxy[2] + FIG_BOUND);
        fprintf(fp, format, llxy[1] - FIG_BOUND, urxy[2] + FIG_BOUND);
        fprintf(fp, format, llxy[1] - FIG_BOUND, llxy[2] - FIG_BOUND);
        fprintf(fp, "\n\n");
    }
}

void get_pxy(double *xy1, double *xy2, double r, double *px, double *py)
/* get the point position of x-, y-, or z-axis */
{
    double dx, dy, dd;

    dx = xy2[1] - xy1[1];
    dy = xy2[2] - xy1[2];
    dd = sqrt(dx * dx + dy * dy);

    if (dd > XEPS) {
        *px = xy2[1] + r * dx / dd;
        *py = xy2[2] + r * dy / dd;
    } else {
        *px = xy2[1];
        *py = xy2[2];
    }
}

void alc2fig(long nobj, long *idx, long *depth, long **allobj, double **blkxyz,
             double **oxyz, long *ibase, long faces[][5], long *opts, FILE * fp)
{
    char *format = "%6.0f%6.0f";
    double dot_sep, msat, Msat, px, py;
    long is_color, same_faces, updown, mgroove;
    long dlcol, dwidth, w1, w2, line_width, join_style, cap_style, o_sides, mfcol;
    long bcol_code, i, ib, ie, ioffset8, ip, j, k, Mfill, mfill;
    long bc_idx[2][7];

    is_color = opts[2];  /* decomposed for clarity */
    same_faces = opts[7];
    updown = opts[9];
    mgroove = opts[10];

    get_fig_pars(&dot_sep, &dlcol, &dwidth, &w1, &w2, bc_idx, &msat, &Msat,
                 &o_sides, &line_width, &join_style, &cap_style, &mfcol);

    /* render each object from inside to outside */
    for (i = 1; i <= nobj; i++) {
        j = idx[i];  /* object number */
        if (allobj[1][j] > 0) {  /* blocks */
            bcol_code = bc_idx[is_color][ibase[allobj[1][j]]];
            k = lround(msat * 20);
            if (mgroove)  /* could still be in color */
                mfill = 20 - k;
            else
                mfill = (bcol_code) ? 20 + k : 20 - k;
            k = lround(Msat * 20);
            if (mgroove)  /* could still be in color */
                Mfill = k;
            else
                Mfill = (bcol_code) ? 40 - k : k;

            if (!same_faces) {
                if ((!updown && !allobj[2][j]) ||  /* minor groove side */
                    (updown && allobj[2][j] == 4))  /* upper face */
                    fprintf(fp, "2 3 0 %2ld %2ld %2ld %4ld 0 %2ld",
                            line_width, bcol_code, bcol_code, depth[i], mfill);
                else if ((!updown && allobj[2][j] == 1) ||  /* major groove side: back */
                         (updown && allobj[2][j] == 5))  /* lower face */
                    fprintf(fp, "2 3 0 %2ld %2ld %2ld %4ld 0 %2ld",
                            line_width, bcol_code, bcol_code, depth[i], Mfill);
                else {  /* other sides */
                    if (bcol_code)  /* color */
                        fprintf(fp, "2 3 0 %2ld %2ld %2ld %4ld 0 20",
                                line_width, bcol_code, o_sides, depth[i]);
                    else  /* black & white */
                        fprintf(fp, "2 3 0 %2ld 0 0 %4ld 0 0", line_width, depth[i]);
                }
            } else {
                if (mgroove && !allobj[2][j])  /* minor groove side */
                    fprintf(fp, "2 3 0 %2ld 0 %2ld %4ld 0 %2ld", line_width,
                            bc_idx[0][0], depth[i], mfcol);
                else
                    fprintf(fp, "2 3 0 %2ld 0 %2ld %4ld 0 %2ld", line_width, bcol_code,
                            depth[i], (bcol_code) ? mfill : Mfill);
            }

            fprintf(fp, " 0.0 %2ld %2ld 0 0 0 5\n", join_style, cap_style);
            ioffset8 = (allobj[1][j] - 1) * 8;
            for (k = 0; k < 5; k++) {
                ip = ioffset8 + faces[allobj[2][j]][k];
                fprintf(fp, format, blkxyz[ip][1], blkxyz[ip][2]);
            }
            fprintf(fp, "\n");
        } else {  /* O connection lines */
            if (!is_color)
                dlcol = 0;  /* black */
            if (!allobj[1][j]) {  /* origin lines */
                fprintf(fp, "2 1 2 %2ld %2ld 0 %4ld 0 -1 %5.1f",
                        dwidth, dlcol, depth[i], dot_sep);
                fprintf(fp, " %2ld %2ld 0 0 0 2\n", join_style, cap_style);
            } else if (allobj[1][j] == -1) {  /* helix axis */
                fprintf(fp, "2 1 0 %2ld %2ld 0 %4ld 0 -1 0.0", w2, dlcol, depth[i]);
                fprintf(fp, " %2ld %2ld 0 0 0 2\n", join_style, cap_style);
            } else {  /* reference axis */
                fprintf(fp, "2 1 0 %2ld %2ld 0 %4ld 0 -1 0.0", w1, dlcol, depth[i]);
                fprintf(fp, " %2ld %2ld 0 1 0 2\n", join_style, cap_style);
                fprintf(fp, "  2 0 2 36 66\n");
            }
            ib = allobj[2][j] / 10000;
            fprintf(fp, format, oxyz[ib][1], oxyz[ib][2]);
            ie = allobj[2][j] % 10000;
            fprintf(fp, format, oxyz[ie][1], oxyz[ie][2]);
            fprintf(fp, "\n");
            if (allobj[1][j] <= -2) {
                get_pxy(oxyz[ib], oxyz[ie], 10.0, &px, &py);
                fprintf(fp, "4 1 %2ld %4ld 0 18 18 0.0 4 165 165", dlcol, depth[i]);
                /* Helvetica-Bold: #18 */
                fprintf(fp, format, px, py);
                if (allobj[1][j] == -2)
                    fprintf(fp, " x\\001\n");
                else if (allobj[1][j] == -3)
                    fprintf(fp, " y\\001\n");
                else
                    fprintf(fp, " z\\001\n");
            }
        }
    }
}

void get_ps_xy(char *imgfile, long *urxy, long frame_box, FILE * fp)
/* reset x/y coordinates to PS units */
{
    char *format = "%6ld%6ld";
    double paper_size[2] = { 8.5, 11.0 };  /* US letter size */
    long i;
    long bbox[5], llxy[3];

    /* centralize the figure on a US letter (8.5in-by-11in) */
    for (i = 1; i <= 2; i++)
        llxy[i] = lround(0.5 * (paper_size[i - 1] * 72 - urxy[i]));

    /* boundary box */
    for (i = 1; i <= 2; i++) {
        bbox[i] = llxy[i] - PS_BOUND;
        bbox[i + 2] = urxy[i] + llxy[i] + PS_BOUND;
    }

    ps_title_cmds(fp, imgfile, bbox);

    fprintf(fp, "%6ld%6ld translate\n\n", llxy[1], llxy[2]);

    if (frame_box) {
        /* draw a box around the figure */
        fprintf(fp, "NP ");
        fprintf(fp, format, -PS_BOUND, -PS_BOUND);
        fprintf(fp, format, urxy[1] + PS_BOUND, -PS_BOUND);
        fprintf(fp, format, urxy[1] + PS_BOUND, urxy[2] + PS_BOUND);
        fprintf(fp, format, -PS_BOUND, urxy[2] + PS_BOUND);
        fprintf(fp, " DB stroke\n\n");
    }
}

void alc2ps(long nobj, long *idx, long **allobj, double **blkxyz, double **oxyz,
            long *ibase, long faces[][5], long *opts, FILE * fp)
{
    char bname;
    char *format = "%7.1f%7.1f";
    char bc_idx[2][7] = { "XXXXXXX",  /* black & white style */
        CX_LIST
    };  /* color-coded for ACGITUX */
    double px, py;
    long is_color, same_faces, updown, mgroove;
    long i, ib, ie, ioffset8, ip, j, k;

    is_color = opts[2];  /* decomposed for clarity */
    same_faces = opts[7];
    updown = opts[9];
    mgroove = opts[10];

    /* render each object from inside to outside */
    for (i = 1; i <= nobj; i++) {
        j = idx[i];  /* object number */
        fprintf(fp, "NP ");
        if (allobj[1][j] > 0) {  /* blocks */
            bname = bc_idx[is_color][ibase[allobj[1][j]]];
            fprintf(fp, "%cl ", (!same_faces) ? bname : 'X');
            ioffset8 = (allobj[1][j] - 1) * 8;
            for (k = 0; k < 4; k++) {
                ip = ioffset8 + faces[allobj[2][j]][k];
                fprintf(fp, format, blkxyz[ip][1], blkxyz[ip][2]);
            }
            fprintf(fp, " DB\n");
            if (!same_faces) {
                if ((!updown && !allobj[2][j]) ||  /* minor groove side */
                    (updown && allobj[2][j] == 4))  /* upper face */
                    fprintf(fp, "  gsave %cm grestore stroke\n", bname);
                else if ((!updown && allobj[2][j] == 1) ||  /* major groove side: back */
                         (updown && allobj[2][j] == 5))  /* lower face */
                    fprintf(fp, "  gsave %cM grestore stroke\n", bname);
                else {  /* other sides */
                    if (bname == 'X')  /* black & white */
                        fprintf(fp, "  gsave 1.0 FB grestore stroke\n");
                    else
                        fprintf(fp, "  gsave OTHER_SIDES grestore stroke\n");
                }
            } else {
                if (mgroove && !allobj[2][j])  /* minor groove side */
                    fprintf(fp, "  gsave Sm grestore stroke\n");
                else {
                    if (bname == 'X')  /* black & white */
                        fprintf(fp, "  gsave XM grestore stroke\n");
                    else
                        fprintf(fp, "  gsave %cm grestore stroke\n", bname);
                }
            }
        } else {  /* O connection lines */
            (is_color) ? fprintf(fp, "Dl ") : fprintf(fp, "Xl ");
            ib = allobj[2][j] / 10000;
            fprintf(fp, format, oxyz[ib][1], oxyz[ib][2]);
            ie = allobj[2][j] % 10000;
            fprintf(fp, format, oxyz[ie][1], oxyz[ie][2]);
            if (!allobj[1][j])  /* origin lines */
                fprintf(fp, "  gsave Dw Ds LN grestore\n");
            else if (allobj[1][j] == -1)  /* helix axis */
                fprintf(fp, "  gsave W2 LN grestore\n");
            else {  /* reference axis */
                fprintf(fp, "  gsave W1 LN grestore\n");

                get_pxy(oxyz[ib], oxyz[ie], 8.0, &px, &py);
                fprintf(fp, format, px, py);
                fprintf(fp, "  moveto");
                if (allobj[1][j] == -2)
                    fprintf(fp, " (x) SCENTER\n");
                else if (allobj[1][j] == -3)
                    fprintf(fp, " (y) SCENTER\n");
                else
                    fprintf(fp, " (z) SCENTER\n");
            }
        }
    }
    fprintf(fp, "\nshowpage\n");
}

void bring_atoms(long ib, long ie, long ra_num, char **AtomName, long *nmatch,
                 long *batom)
/* get base ring atom index in one residue */
{
    static char *RingAtom[] = { RA_LIST };
    long i, j;

    *nmatch = 0;

    for (i = 0; i < ra_num; i++) {
        j = find_1st_atom(RingAtom[i], AtomName, ib, ie, "in base ring atoms");
        if (j)
            batom[++*nmatch] = j;
    }
}

void all_bring_atoms(long num_residue, long *RY, long **seidx, char **AtomName,
                     long *num_ring, long **ring_atom)
/* get base ring atom index for all residues: num_ring */
{
    long i, j, nmatch;

    *num_ring = 0;
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0) {  /* non-base residue */
            ring_atom[i][10] = -1;
            continue;
        }
        j = (RY[i] == 1) ? 9 : 6;
        bring_atoms(seidx[i][1], seidx[i][2], j, AtomName, &nmatch, ring_atom[i]);
        if (nmatch == j) {
            ring_atom[i][10] = j;
            ++*num_ring;
        }
    }
}

void base_idx(long num, char *bseq, long *ibase, long single)
/* get base index for coloring purpose */
{
    static char *cmn_base = CB_LIST, *pchar;
    long i;

    if (single) {  /* for a single case */
        if ((pchar = strchr(cmn_base, toupper(*bseq))) != NULL)
            *ibase = pchar - cmn_base;
        else
            *ibase = NON_WC_IDX;
    } else {
        for (i = 1; i <= num; i++)
            if ((pchar = strchr(cmn_base, toupper(bseq[i]))) != NULL)
                ibase[i] = pchar - cmn_base;
            else
                ibase[i] = NON_WC_IDX;
    }
}

long basepair_idx(char *bpi)
/* get basepair index for coloring purpose */
{
    static char *WC[9] = { WC_LIST };
    long i, bidx;

    i = find_1st_atom(bpi, WC, 1, 8, "");
    if (!i)
        bidx = NON_WC_IDX;  /* non A-T/G-C pair */
    else if (i <= 2)  /* AT or AU */
        bidx = 0;  /* 0 for A */
    else if (i <= 4)  /* TA or UA */
        bidx = 4;  /* 4 for T */
    else if (i <= 6)  /* GC or IC */
        bidx = 2;  /* 2 for G */
    else  /* CG or CI */
        bidx = 1;  /* 1 for C */

    return bidx;
}

void plane_xyz(long num, double **xyz, double *ppos, double *nml, double **nxyz)
/* given plane normal and its center, project all coordinates onto it */
{
    long i, j;
    double temp, d[4];

    for (i = 1; i <= num; i++) {
        ddxyz(ppos, xyz[i], d);
        temp = dot(d, nml);
        for (j = 1; j <= 3; j++)
            nxyz[i][j] = ppos[j] + d[j] - temp * nml[j];
    }
}

void prj2plane(long num, long ra_num, char **AtomName, double **xyz, double z0,
               double **nxyz)
/* project base atoms onto its least-squares plane defined by ring atoms */
{
    double ang, temp;
    double zaxis[4] = { EMPTY_NUMBER, 0.0, 0.0, 1.0 };
    double adist[10], hinge[4], ppos[4], z[4];
    double **bxyz, **rmtx;

    long i, j, nmatch;
    long batom[10];

    /* find base ring atoms */
    bring_atoms(1, num, ra_num, AtomName, &nmatch, batom);

    /* base least-squares plane */
    bxyz = dmatrix(1, nmatch, 1, 3);

    for (i = 1; i <= nmatch; i++)
        cpxyz(xyz[batom[i]], bxyz[i]);
    ls_plane(bxyz, nmatch, z, ppos, &temp, adist);

    /* get the new set of coordinates */
    plane_xyz(num, xyz, ppos, z, nxyz);

    /* reorient the structure to make z-coordinate zero */
    rmtx = dmatrix(1, 3, 1, 3);
    if (z0) {
        cross(z, zaxis, hinge);
        ang = magang(z, zaxis);
        arb_rotation(hinge, ang, rmtx);
        for (i = 1; i <= num; i++)
            for (j = 1; j <= 3; j++)
                nxyz[i][j] = dot(xyz[i], rmtx[j]);
        for (i = 1; i <= num; i++)
            nxyz[i][3] -= nxyz[1][3];
    }
    free_dmatrix(bxyz, 1, nmatch, 1, 3);
    free_dmatrix(rmtx, 1, 3, 1, 3);
}

void adjust_xy(long num, double **xyz, long nO, double **oxyz, double scale_factor,
               long default_size, long *urxy)
/* reset x & y coordinates to fit the scale */
{
    long i, j;
    double temp;
    double max_xy[3], min_xy[3];

    max_dmatrix(xyz, num, 2, max_xy);
    min_dmatrix(xyz, num, 2, min_xy);

    /* get maximum dx or dy */
    temp = dval_max(max_xy[1] - min_xy[1], max_xy[2] - min_xy[2]);

    scale_factor = fabs(scale_factor);
    if (scale_factor < XEPS)
        scale_factor = default_size / temp;

    fprintf(stderr, "\n ...... scale factor: %.2f ...... \n", scale_factor);

    move_position(xyz, num, 2, min_xy);
    for (i = 1; i <= num; i++)
        for (j = 1; j <= 2; j++)
            xyz[i][j] *= scale_factor;
    if (nO) {
        move_position(oxyz, nO, 2, min_xy);
        for (i = 1; i <= nO; i++)
            for (j = 1; j <= 2; j++)
                oxyz[i][j] *= scale_factor;
    }
    max_dmatrix(xyz, num, 2, max_xy);
    for (i = 1; i <= 2; i++)
        urxy[i] = lround(max_xy[i]);
}

void get_depth(long nobj, long *zval, long *depth)
{
    double temp;
    long depth_low = 991, depth_up = 11;  /* depth level */
    long i, j;

    /* reset zval to [depth_low -- depth_up] for depth level */
    j = depth_low - depth_up;
    temp = zval[nobj] - zval[1];
    for (i = 1; i <= nobj; i++)
        depth[i] = lround(depth_low - j * (zval[i] - zval[1]) / temp);
}

void raster3d_header(long num, double **xyz, double scale_factor, long no_header,
                     long frame_box, FILE * fp)
{
    char BDIR[BUF512], str[BUF512], *header_file = "my_header.r3d";
    double temp, ave_xyz[4], min_xyz[4], max_xyz[4];
    double rad = 0.06, rgbv[4] = { 0.0, 0.25, 0.25, 0.25 };
    long i, itype = 3;
    FILE *fpp;

    min_dmatrix(xyz, num, 3, min_xyz);
    max_dmatrix(xyz, num, 3, max_xyz);
    avexyz(max_xyz, min_xyz, ave_xyz);
    temp = dval_max(max_xyz[1] - min_xyz[1], max_xyz[2] - min_xyz[2]);

    if (!no_header) {  /* write header section */
        if (scale_factor < XEPS)
            scale_factor = 6.0 + temp;
        fprintf(stderr, "\n ...... scale factor: %.2f ...... \n", scale_factor);

        get_BDIR(BDIR, header_file, 1);
        strcat(BDIR, header_file);
        fprintf(stderr, " ...... reading file: %s ...... \n", header_file);

        fpp = open_file(BDIR, "r");
        for (i = 1; i <= 20; i++) {
            if (fgets(str, sizeof str, fpp) == NULL)
                fatal("error reading header.r3d\n");
            (i != 16) ? fputs(str, fp) :
                fprintf(fp, "%9.3f%9.3f%9.3f%9.3f\n", -ave_xyz[1],
                        -ave_xyz[2], -ave_xyz[3], scale_factor);
        }
        close_file(fpp);
    }

    if (frame_box) {
        fprintf(fp, "### Section of the frame box: 4 lines\n");
        ave_xyz[3] = min_xyz[3] = max_xyz[3];
        /* box strongly influenced by z-coordinate */
        ave_xyz[1] = max_xyz[1];  /* lower-right corner */
        ave_xyz[2] = min_xyz[2];
        r3d_rod(itype, min_xyz, ave_xyz, rad, rgbv, fp);
        r3d_rod(itype, ave_xyz, max_xyz, rad, rgbv, fp);
        ave_xyz[1] = min_xyz[1];  /* upper-left corner */
        ave_xyz[2] = max_xyz[2];
        r3d_rod(itype, min_xyz, ave_xyz, rad, rgbv, fp);
        r3d_rod(itype, ave_xyz, max_xyz, rad, rgbv, fp);
    }
}

void get_r3dpars(double **base_col, double *hb_col, double *width3, double **atom_col,
                 char *label_style)
/* read in parameters for Raster3D input */
{
    char BDIR[BUF512], str[BUF512], *raster3d_par = "raster3d.par";
    char *format = "%lf %lf %lf", *format4 = "%lf %lf %lf %lf";
    long i;
    FILE *fp;

    get_BDIR(BDIR, raster3d_par, 1);
    strcat(BDIR, raster3d_par);
    fp = open_file(BDIR, "r");
    fprintf(stderr, " ...... reading file: %s ...... \n", raster3d_par);

    if (fgets(str, sizeof str, fp) == NULL)  /* skip one line */
        fatal("error in reading comment line\n");
    for (i = 0; i <= NBASECOL; i++)
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &base_col[i][1], &base_col[i][2], &base_col[i][3]) != 3)
            fatal("error reading base residue RGB color\n");
    if (fgets(str, sizeof str, fp) == NULL)  /* skip one line */
        fatal("error in reading comment line\n");
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, format4, &hb_col[1], &hb_col[2], &hb_col[3], &hb_col[4]) != 4)
        fatal("error reading H-bond RGB color\n");
    if (fgets(str, sizeof str, fp) == NULL)  /* skip one line */
        fatal("error in reading comment line\n");
    if (fgets(str, sizeof str, fp) == NULL ||
        sscanf(str, "%lf %lf %lf", &width3[1], &width3[2], &width3[3]) != 3)
        fatal("error cylinder radius for bp-center line & bp 1 & 2\n");
    if (fgets(str, sizeof str, fp) == NULL)  /* skip one line */
        fatal("error in reading comment line\n");
    for (i = 0; i <= NATOMCOL; i++)
        if (fgets(str, sizeof str, fp) == NULL ||
            sscanf(str, format, &atom_col[i][1], &atom_col[i][2], &atom_col[i][3]) != 3)
            fatal("error reading atom RGB color\n");
    if (fgets(str, sizeof str, fp) == NULL)  /* skip one line */
        fatal("error in reading comment line\n");
    if (fgets(str, sizeof str, fp) == NULL)
        fatal("error reading label style\n");
    strcpy(label_style, str);

    close_file(fp);
}

void r3d_rod(long itype, double *xyz1, double *xyz2, double rad, double *rgbv, FILE * fp)
/* write a record of round-ended cylinder (itype = 3) or flat-ended (5) or comments */
{
    static char *format = "%9.3f";
    long i;

    if (itype == 3 || itype == 5)
        fprintf(fp, "%ld\n", itype);
    else
        fprintf(fp, "#5\n#");  /* comments and default to type 5 */
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz1[i]);
    fprintf(fp, format, rad);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz2[i]);
    fprintf(fp, format, rad);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, rgbv[i]);
    fprintf(fp, "\n");
}

void r3d_dash(double *xyz1, double *xyz2, double hb_width, double *hb_col, FILE * fp)
/* draw dashed line for H-bonds, and lines connecting base-pair centers */
{
    double distance, dnum, dxyz[4], m1[4], m2[4];
    long i, j, num;

    ddxyz(xyz1, xyz2, dxyz);
    distance = veclen(dxyz);
    vec_norm(dxyz);  /* normal vector */

    num = lround(hb_col[4] * distance);
    if (num % 2 == 0)
        num++;  /* odd "num" for ending on last point */
    if (num <= 1)  /* solid line */
        r3d_rod(5, xyz1, xyz2, hb_width, hb_col, fp);
    else {
        dnum = (double) num;
        for (i = 0; i < num; i += 2) {
            for (j = 1; j <= 3; j++) {
                m1[j] = xyz1[j] + dxyz[j] * distance * i / dnum;
                m2[j] = xyz1[j] + dxyz[j] * distance * (i + 1) / dnum;
            }
            r3d_rod(5, m1, m2, hb_width, hb_col, fp);
        }
    }
}

void r3d_sphere(double *xyz1, double rad, double *rgbv, FILE * fp)
/* write a record of sphere */
{
    static char *format = "%9.3f";
    long i;

    fprintf(fp, "2\n");
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz1[i]);
    fprintf(fp, format, rad);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, rgbv[i]);
    fprintf(fp, "\n");
}

void cpk_model(long num, long *idx, double **xyz, double ballrad, double **colrgb,
               FILE * fp)
{
    double vdw_radii[NELE];
    long i;

    if (ballrad <= 0.0)
        return;

    atom_info(3, NULL, NULL, vdw_radii);
    fprintf(fp, "###\n### The following section is for ball/CPK model\n");
    for (i = 1; i <= num; i++)
        r3d_sphere(xyz[i], vdw_radii[idx[i]] * ballrad, colrgb[i], fp);
}

void r3d_tripln(long itype, double *xyz1, double *xyz2, double *xyz3, double *rgbv,
                FILE * fp)
/* write a record of triangle (itype = 1, default) or plane (itype = 6) */
{
    static char *format = "%9.3f";
    long i;

    fprintf(fp, "%d\n", (itype != 6) ? 1 : 6);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz1[i]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz2[i]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, xyz3[i]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, rgbv[i]);
    fprintf(fp, "\n");
}

void r3d_block_edge(double *rgbv, long ioffset8, double **blkxyz, double w1, FILE * fp)
{
    static long blk_lkg[12][2] = {
        {1, 2},
        {2, 3},
        {3, 4},
        {1, 4},
        {1, 5},
        {2, 6},
        {3, 7},
        {4, 8},
        {5, 6},
        {6, 7},
        {7, 8},
        {5, 8}
    };
    long j;

    for (j = 0; j < 12; j++)
        r3d_rod(3, blkxyz[ioffset8 + blk_lkg[j][0]],
                blkxyz[ioffset8 + blk_lkg[j][1]], w1, rgbv, fp);
}

void base_label(double **rxyz, char *label_style, double *rgbv, char *bname_num,
                FILE * fp)
/* label the base in the center of six-membered ring */
{
    static char *format = "%9.3f";
    double cxyz[4];
    long i;

    avexyz(rxyz[1], rxyz[4], cxyz);  /* N1 + C4 */
    fprintf(fp, "10\n%s", label_style);  /* label_style has \n */
    fprintf(fp, "11\n");
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, cxyz[i]);
    for (i = 1; i <= 3; i++)
        fprintf(fp, format, rgbv[i]);
    fprintf(fp, "\n%s\n", bname_num);
}

void fill_base_ring(long num_residue, long num_ring, long **ring_atom, double **xyz,
                    long *ibase, char *bseq, double **base_col, char *label_style,
                    long label_ring, long *ResSeq, FILE * fp)
{
    char bname_num[BUF512];
    double **rxyz;
    long i, j;

    rxyz = dmatrix(1, 9, 1, 3);  /* base ring atoms */

    fprintf(fp, "###\n### The following section is for %ld filled base rings\n",
            num_ring);
    for (i = 1; i <= num_residue; i++) {
        if (ring_atom[i][10] <= 0)
            continue;
        for (j = 1; j <= ring_atom[i][10]; j++)
            cpxyz(xyz[ring_atom[i][j]], rxyz[j]);
        r3d_tripln(1, rxyz[3], rxyz[4], rxyz[6], base_col[ibase[i]], fp);
        r3d_tripln(1, rxyz[1], rxyz[2], rxyz[3], base_col[ibase[i]], fp);
        r3d_tripln(1, rxyz[4], rxyz[5], rxyz[6], base_col[ibase[i]], fp);
        r3d_tripln(1, rxyz[1], rxyz[3], rxyz[6], base_col[ibase[i]], fp);
        if (ring_atom[i][10] == 9) {
            r3d_tripln(1, rxyz[6], rxyz[7], rxyz[8], base_col[ibase[i]], fp);
            r3d_tripln(1, rxyz[1], rxyz[8], rxyz[9], base_col[ibase[i]], fp);
            r3d_tripln(1, rxyz[1], rxyz[6], rxyz[8], base_col[ibase[i]], fp);
        }
        if (label_ring) {
            sprintf(bname_num, "%c%ld", bseq[i], ResSeq[ring_atom[i][1]]);
            base_label(rxyz, label_style, base_col[ibase[i]], bname_num, fp);
        }
    }

    free_dmatrix(rxyz, 1, 9, 1, 3);
}

void process_alc(char *alcfile, char *imgfile, double scale_factor, long *opts)
{
    char **AtomName;
    double **blkxyz, **oxyz, **xyz;
    long i, j, k, nbond, nobj, npoint_blk, num, num_blk = 0, num_blk8;
    long ioffset8, ioffset_blk, nN = 0, nO = 0, nO_all = 0, nO_lkg = 0, nO_tmp;
    long Nidx[6], *blkibase, *ibase, **linkage;

    /* read in ALCHEMY file */
    get_alc_nums(alcfile, &num, &nbond);
    AtomName = cmatrix(1, num, 0, 2);
    xyz = dmatrix(1, num, 1, 3);
    ibase = lvector(1, num);
    linkage = lmatrix(1, nbond, 0, 2);
    read_alc(alcfile, &num, &nbond, AtomName, xyz, ibase, linkage);

    /* assume the first 8 atoms are 4 Ns followed by 4 Cs in
       a specific numbering scheme */
    for (i = 1; i <= num; i++)
        if (AtomName[i][0] == 'N')
            ++nN;
    if (!nN || nN % 4)
        fatal("wrong type of standard block file\n");
    else
        num_blk = nN / 4;

    nN = 0;
    for (i = 1; i <= num && nN < 5; i++)
        if (AtomName[i][0] == 'N')
            Nidx[++nN] = i;
    npoint_blk = (num_blk > 1) ? Nidx[5] - Nidx[1] : num;

    /* get the xyz coordinates of 8-point-per-block for each */
    num_blk8 = num_blk * 8;
    blkxyz = dmatrix(1, num_blk8, 1, 3);
    blkibase = lvector(1, num_blk);
    for (i = 1; i <= num_blk; i++) {
        ioffset8 = (i - 1) * 8;
        ioffset_blk = (i - 1) * npoint_blk;
        blkibase[i] = ibase[ioffset_blk + 1];
        if (!lval_in_range(blkibase[i], 0, 5)) {
            fprintf(stderr, "uncommon base. set to default color coding\n");
            blkibase[i] = NON_WC_IDX;
        }
        for (j = 1; j <= 8; j++)
            cpxyz(xyz[ioffset_blk + j], blkxyz[ioffset8 + j]);
    }
    nobj = num_blk * 6;  /* blocks */

    /* get the total number of O atoms:
       "O ", "OH", "OO", "OX", "OY", "OZ" */
    for (i = 1; i <= num; i++)
        if (AtomName[i][0] == 'O')
            nO_all++;
    oxyz = dmatrix(1, nO_all, 1, 3);

    if (opts[3]) {  /* draw a line to link neighbor origins */
        nO_tmp = 0;
        for (i = 1; i <= num; i++)
            if (!strcmp(AtomName[i], "O ")) {
                ++nO_tmp;
                cpxyz(xyz[i], oxyz[nO_tmp]);
            }
        if (nO_tmp) {
            k = nO_tmp - 1;
            for (i = 1; i <= k; i++) {
                linkage[i][0] = 0;
                linkage[i][1] = i;
                linkage[i][2] = i + 1;
            }
            nO_lkg = k;
            nO = nO_tmp;
        }
    }
    if (opts[5]) {  /* draw helix-axis */
        nO_tmp = 0;
        for (i = 1; i <= num; i++)
            if (!strcmp(AtomName[i], "OH")) {
                ++nO_tmp;
                cpxyz(xyz[i], oxyz[nO + nO_tmp]);
            }
        if (nO_tmp) {
            if (nO_tmp % 2)
                fatal("wrong number of helix axis origins\n");
            k = nO_tmp / 2;
            for (i = 1; i <= k; i++) {
                j = nO_lkg + i;
                linkage[j][0] = -1;
                linkage[j][1] = nO + 2 * i - 1;
                linkage[j][2] = nO + 2 * i;
            }
            nO_lkg += k;
            nO += nO_tmp;
        }
    }
    if (opts[6]) {  /* draw global reference axis */
        nO_tmp = 0;
        for (i = 1; i <= num; i++)
            if (!strcmp(AtomName[i], "OO") || !strcmp(AtomName[i], "OX") ||
                !strcmp(AtomName[i], "OY") || !strcmp(AtomName[i], "OZ")) {
                ++nO_tmp;
                cpxyz(xyz[i], oxyz[nO + nO_tmp]);
            }
        if (nO_tmp) {
            k = nO_tmp - 1;
            for (i = 1; i <= k; i++) {
                j = nO_lkg + i;
                linkage[j][0] = -(i + 1);
                linkage[j][1] = nO + 1;
                linkage[j][2] = nO + i + 1;
            }
            nO_lkg += k;
            nO += nO_tmp;
        }
    }

    alc_3images(opts, nobj, num_blk, num_blk8, nO_lkg, nO, oxyz, blkxyz, blkibase,
                linkage, scale_factor, imgfile);

    free_dmatrix(oxyz, 1, nO_all, 1, 3);
    free_alc(num, nbond, AtomName, xyz, ibase, 0, linkage);
    free_dmatrix(blkxyz, 1, num_blk8, 1, 3);
    free_lvector(blkibase, 1, num_blk);
}

void alc_3images(long *opts, long nobj, long num_blk, long num_blk8, long nO_lkg,
                 long nO, double **oxyz, double **blkxyz, long *blkibase, long **linkage,
                 double scale_factor, char *imgfile)
/* From ALCHEMY format to three image formats: PS, XFIG, Raster3D */
{
    char label_style[BUF512];
    double width3[4], hb_col[5], *pd, **atom_col, **base_col;
    static long faces[6][5] = {
        {1, 2, 3, 4, 1},  /* 0 front: minor groove */
        {5, 6, 7, 8, 5},  /* 1 back */
        {1, 4, 8, 5, 1},  /* 2 left */
        {2, 3, 7, 6, 2},  /* 3 right */
        {1, 2, 6, 5, 1},  /* 4 upper */
        {3, 4, 8, 7, 3}  /*  5 lower */
    };
    long default_size[2] = { PS_DFTSIZE, FIG_DFTSIZE };  /* PS & XFIG */
    long is_color, same_faces, updown, mgroove;
    long i, j, k, ioffset8, urxy[3];
    long *depth, *idx, **allobj;
    FILE *fp;

    is_color = opts[2];  /* decomposed for clarity: cf. alc2ps/alc2fig */
    same_faces = opts[7];
    updown = opts[9];
    mgroove = opts[10];

    fp = open_file(imgfile, "w");
    if (opts[0] == 2) {  /* Raster3D */
        atom_col = dmatrix(0, NATOMCOL, 1, 3);
        base_col = dmatrix(0, NBASECOL, 1, 3);
        raster3d_header(num_blk8, blkxyz, scale_factor, opts[8], opts[1], fp);
        /* only base_col is used here */
        get_r3dpars(base_col, hb_col, width3, atom_col, label_style);
        for (i = 1; i <= num_blk; i++) {
            ioffset8 = (i - 1) * 8;
            k = is_color ? blkibase[i] : NON_WC_IDX;  /* if color image */
            r3d_block_edge(same_faces ? base_col[NON_WC_IDX] : base_col[k],  /* if same-face */
                           ioffset8, blkxyz, width3[2], fp);  /* bp1 width */
            for (j = 0; j < 6; j++) {
                if (same_faces)
                    pd = (mgroove && !j) ? base_col[NON_WC_IDX] : base_col[k];
                else {
                    if ((!updown && !j) || (updown && j == 4))
                        pd = base_col[k];
                    else
                        pd = base_col[NBASECOL];
                }
                r3d_tripln(1, blkxyz[ioffset8 + faces[j][0]], blkxyz[ioffset8 + faces[j][1]], blkxyz[ioffset8 + faces[j][2]], pd, fp);  /* 0-1-2 */
                r3d_tripln(1, blkxyz[ioffset8 + faces[j][0]], blkxyz[ioffset8 + faces[j][2]], blkxyz[ioffset8 + faces[j][3]], pd, fp);  /* 0-2-3 */
            }
        }
        if (!is_color)
            cpxyz(base_col[NON_WC_IDX], hb_col);
        for (i = 1; i <= nO_lkg; i++) {
            if (!linkage[i][0])  /* bp-center connecting line */
                r3d_dash(oxyz[linkage[i][1]], oxyz[linkage[i][2]], width3[1], hb_col, fp);
            else  /* helix or reference frame */
                r3d_rod(3, oxyz[linkage[i][1]], oxyz[linkage[i][2]],
                        (linkage[i][0] == -1) ? width3[3] : width3[2], hb_col, fp);
        }

        free_dmatrix(atom_col, 0, NATOMCOL, 1, 3);
        free_dmatrix(base_col, 0, NBASECOL, 1, 3);
    } else {
        nobj += nO_lkg;  /* total number of objects */

        /* combine all objects (faces + origin lines) together */
        allobj = lmatrix(1, 3, 1, nobj);
        get_alc_objs(num_blk, blkxyz, nO, oxyz, nO_lkg, linkage, faces, allobj);

        idx = lvector(1, nobj);
        lsort(nobj, allobj[3], idx);

        adjust_xy(num_blk8, blkxyz, nO, oxyz, scale_factor, default_size[opts[0]], urxy);

        if (opts[0] == 1) {  /* for XFIG */
            depth = lvector(1, nobj);
            get_depth(nobj, allobj[3], depth);
            get_fig_xy(num_blk8, blkxyz, nO, oxyz, urxy, opts[1], fp);
            alc2fig(nobj, idx, depth, allobj, blkxyz, oxyz, blkibase, faces, opts, fp);
            free_lvector(depth, 1, nobj);
        } else {  /* PS */
            get_ps_xy(imgfile, urxy, opts[1], fp);
            alc2ps(nobj, idx, allobj, blkxyz, oxyz, blkibase, faces, opts, fp);
        }
        free_lmatrix(allobj, 1, 3, 1, nobj);
        free_lvector(idx, 1, nobj);
    }
    close_file(fp);
}

void get_alc_objs(long num_blk, double **blkxyz, long nO, double **oxyz, long nO_lkg,
                  long **linkage, long faces[][5], long **allobj)
/* "allobj" has three row as follows
   (a) for blocks        [1] block index (1 -- num_blk)
   +                     [2] face index (0 -- 5)
   +                     [3] average z-coord x 1000 for each face
   (b) for O linkages    [1] 0: origin, -1: helix, -2, -3, -4: frame x, y, z
   +                     [2] 10000 x id1 + id2
   +                     [3] average z-coord x 1000 for the line
 */
{
    long i, ioffset8, ip, j, k;

    for (i = 1; i <= num_blk; i++) {
        ioffset8 = (i - 1) * 8;  /* ith block atom index offset */
        k = (i - 1) * 6 + 1;  /* ith block face offset */
        for (j = 0; j < 6; j++) {
            ip = k + j;
            allobj[1][ip] = i;
            allobj[2][ip] = j;
            /* only two diagonal atoms are needed for z-coordinates */
            allobj[3][ip] = lround(500.0 * (blkxyz[ioffset8 + faces[j][0]][3] +
                                            blkxyz[ioffset8 + faces[j][2]][3]));
        }
    }

    if (nO) {
        k = num_blk * 6;  /* current object offset */
        for (i = 1; i <= nO_lkg; i++) {
            j = k + i;
            allobj[1][j] = linkage[i][0];
            allobj[2][j] = linkage[i][1] * 10000 + linkage[i][2];
            allobj[3][j] = lround(500.0 * (oxyz[linkage[i][1]][3] +
                                           oxyz[linkage[i][2]][3]));
        }
    }
}

void get_fig_pars(double *dot_sep, long *dlcol, long *dwidth, long *bp1width,
                  long *bp2width, long bc_idx[][7], double *msat, double *Msat,
                  long *o_sides, long *line_width, long *join_style, long *cap_style,
                  long *mfcol)
/* read in rendering parameters for XFIG file */
{
    char BDIR[BUF512], str[BUF512];
    char *fig_image_par = "fig_image.par";
    long i;
    FILE *fp;

    get_BDIR(BDIR, fig_image_par, 1);
    strcat(BDIR, fig_image_par);
    fprintf(stderr, " ...... reading file: %s ...... \n", fig_image_par);

    fp = open_file(BDIR, "r");
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%lf", dot_sep) != 1)
        fatal("error in reading dot separation\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", dlcol) != 1)
        fatal("error in reading dot line color\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", dwidth) != 1)
        fatal("error in reading dot line width\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", bp1width) != 1)
        fatal("error in reading bp1 width\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", bp2width) != 1)
        fatal("error in reading bp2 width\n");

    for (i = 0; i < 7; i++)
        if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", &bc_idx[1][i]) != 1)
            fatal("error in reading color-code\n");
    for (i = 0; i < 7; i++)  /* black & white option */
        bc_idx[0][i] = 0;

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%lf", msat) != 1)
        fatal("error in reading minor groove color saturation\n");
    *msat = fabs(*msat);
    if (*msat > 1.0)
        *msat = 0.9;
    *msat = 1.0 - *msat;
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%lf", Msat) != 1)
        fatal("error in reading major groove color saturation\n");
    *Msat = fabs(*Msat);
    if (*Msat > 1.0)
        *Msat = 0.1;
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", o_sides) != 1)
        fatal("error in reading other-side color code\n");

    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", line_width) != 1)
        fatal("error in reading line width\n");
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", join_style) != 1)
        fatal("error in reading join style\n");
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", cap_style) != 1)
        fatal("error in reading cap style\n");
    if (fgets(str, sizeof str, fp) == NULL || sscanf(str, "%ld", mfcol) != 1)
        fatal("error in reading minor groove filling color\n");

    close_file(fp);
}

void change_xyz(long side_view, double *morg, double **mst, long num, double **xyz)
{
    double **tmpxyz;

    tmpxyz = dmatrix(1, num, 1, 3);

    move_position(xyz, num, 3, morg);
    multi_matrix(xyz, num, 3, mst, 3, 3, tmpxyz);

    copy_matrix(tmpxyz, num, 3, xyz);

    if (side_view)
        get_side_view(1, num, xyz);

    free_dmatrix(tmpxyz, 1, num, 1, 3);
}

void get_side_view(long ib, long ie, double **xyz)
/* adjust orientation: xyz*rotz(-90)*rotx(90) ==> [-y z -x] */
{
    long i;
    double temp;

    for (i = ib; i <= ie; i++) {
        temp = xyz[i][1];
        xyz[i][1] = -xyz[i][2];
        xyz[i][2] = xyz[i][3];
        xyz[i][3] = -temp;
    }
}

void get_CNidx(long ds, long num_bp, long **chi, long **idx, long *nvec, long *C1b,
               long *C1e)
/* get C1* and RN9/YN1 atom pair index for finding helix */
{
    long i, ia, ib, ioffset, j, joffset, k;

    *nvec = 0;
    for (i = 1; i <= ds; i++)
        for (j = 1; j <= num_bp - 1; j++) {
            ioffset = (j - 1) * 4;
            joffset = j * 4;
            for (k = 2; k <= 3; k++) {
                ia = chi[i][ioffset + k];
                ib = chi[i][joffset + k];
                if (ia && ib) {
                    idx[++*nvec][1] = ia;
                    idx[*nvec][2] = ib;
                }
            }
        }

    ioffset = (num_bp - 1) * 4;
    for (i = 1; i <= ds; i++) {
        C1b[i] = chi[i][2];
        C1e[i] = chi[i][ioffset + 2];
    }
}

void add_3axes(long *num, char **AtomName, long *ibase, double **xyz, long *nbond,
               long **linkage, long side_view)
/* attach reference frame axes */
{
    static char *ref_symbol[4] = { "OO", "OX", "OY", "OZ" };
    static double ref_xyz[4][3] = {
        {0.0, 0.0, 0.0},  /* origin */
        {1.0, 0.0, 0.0},  /* x-axis */
        {0.0, 1.0, 0.0},  /* y-axis */
        {0.0, 0.0, 1.0}  /* z-axis */
    };
    long i, j, k;

    for (i = 1; i <= 4; i++) {
        k = *num + i;
        strcpy(AtomName[k], ref_symbol[i - 1]);
        for (j = 1; j <= 3; j++)
            xyz[k][j] = AXIS_LENGTH * ref_xyz[i - 1][j - 1];
        ibase[k] = NON_WC_IDX;  /* non-common base atoms */
    }
    if (side_view)
        get_side_view(*num + 1, *num + 4, xyz);
    for (i = 1; i <= 3; i++) {
        k = *nbond + i;
        linkage[k][1] = *num + 1;
        linkage[k][2] = *num + i + 1;
    }
    *num += 4;
    *nbond += 3;
}

void bp_sequence(char *Wbase, char *Cbase, long *num_bp, char **bp_seq)
{
    char c, *pchar;
    char str[BUF512], crepeat[BUF512];
    long i, j, k, ich = 0, n = 0, nbp = 0;
    FILE *fp;

    fprintf(stderr, "\nInput your base sequence with only A,C,G & T:\n");
    fprintf(stderr, "1. From a data file (complete sequence)\n");
    fprintf(stderr, "2. From keyboard (enter only the repeating sequence)\n");
    fprintf(stderr, "Your choice (1 or 2, Dft: 2): ");
    if (fgets(str, sizeof str, stdin) != NULL) {
        j = sscanf(str, "%ld", &ich);
        if (!j || j == EOF)
            ich = 2;
    } else
        fatal("error in reading your choice\n");
    fprintf(stderr, "\n");

    if (ich == 1) {
        fprintf(stderr, "Name of your base sequence file: ");
        if (fgets(str, sizeof str, stdin) == NULL)
            fatal("error in reading your sequence file name\n");
        str[strlen(str) - 1] = '\0';  /* change \n to \0 */
        fp = open_file(str, "r");
        while (fgets(str, sizeof str, fp) != NULL)
            for (i = 0; i < (long) strlen(str); i++) {
                c = (islower((int) str[i])) ? toupper(str[i]) : str[i];
                if ((pchar = strchr(Wbase, c)) != NULL) {
                    if (++nbp > MAXBASE)
                        fatal("number of bases out of limit\n");
                    bp_seq[nbp][1] = c;
                } else if (!isspace((int) c))
                    fprintf(stderr, "ignore uncommon base %c\n", c);
            }
        close_file(fp);
        if (!nbp)
            fatal("your input file contains no common bases\n");
    } else {
        fprintf(stderr, "Repeating unit (Dft: polyA): ");
        if (fgets(str, sizeof str, stdin) == NULL)
            fatal("error in reading your repeating unit\n");
        for (i = 0; i < (long) strlen(str); i++) {
            c = (islower((int) str[i])) ? toupper(str[i]) : str[i];
            if ((pchar = strchr(Wbase, c)) != NULL)
                crepeat[nbp++] = c;
            else if (!isspace((int) c))
                fprintf(stderr, "ignore uncommon base %c\n", c);
        }
        if (!nbp) {
            strcpy(crepeat, "A");  /* default for ployA */
            nbp = 1;
        } else
            crepeat[nbp] = '\0';
        fprintf(stderr, "Repeating unit: %s\n", crepeat);

        fprintf(stderr, "Number of repeats (Dft: 10): ");
        if (fgets(str, sizeof str, stdin) != NULL) {
            j = sscanf(str, "%ld", &n);
            if (!j || j == EOF)
                n = 10;
        } else
            fatal("error in reading number of repeats\n");
        if (n < 0) {
            fprintf(stderr, "number of repeats < 0. reset to positive\n");
            n = -n;
        }
        if (n * nbp > MAXBASE)
            fatal("number of bases out of limit\n");
        for (i = 1; i <= n; i++) {
            k = (i - 1) * nbp + 1;
            for (j = 0; j < nbp; j++)
                bp_seq[k + j][1] = crepeat[j];
        }
        nbp *= n;
    }

    for (i = 1; i <= nbp; i++) {
        j = strchr(Wbase, bp_seq[i][1]) - Wbase;
        bp_seq[i][2] = *(Cbase + j);
    }

    *num_bp = nbp;
}

void combine_pstnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                    char **tAtomName, char **tResName, char *tChainID, long *tResSeq,
                    double **txyz, char **tAtomName2, double **txyz2)
/* combine strand II with strand I for a parallel duplex */
{
    long i, ik, j;

    for (i = 1; i <= num_bp; i++) {
        for (j = s2idx[i][1] + 1; j <= s2idx[i][1] + s2idx[i][2]; j++) {
            ik = *tnum + j - s2idx[i][1];
            strcpy(tAtomName[ik], tAtomName2[j]);
            sprintf(tResName[ik], "  %c", bp_seq[i][2]);
            tChainID[ik] = 'B';
            tResSeq[ik] = num_bp + i;
            cpxyz(txyz2[j], txyz[ik]);
        }
        *tnum += s2idx[i][2];
    }
}

void reverse_stnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                   char **tAtomName, char **tResName, char *tChainID, long *tResSeq,
                   double **txyz, char **tAtomName2, double **txyz2, long basep)
/* reverse strand II to be in 5'-->3' direction & combined with I */
{
    long i, ik, j;

    for (i = num_bp; i >= 1; i--) {
        for (j = s2idx[i][1] + 1; j <= s2idx[i][1] + s2idx[i][2]; j++) {
            ik = *tnum + j - s2idx[i][1];
            strcpy(tAtomName[ik], tAtomName2[j]);
            sprintf(tResName[ik], "  %c", bp_seq[i][2]);
            tChainID[ik] = 'B';
            tResSeq[ik] = 2 * num_bp - i + 1;
            if (basep && !strcmp(tAtomName[ik], " P  "))
                ++tResSeq[ik];  /* re-position P on strand II */
            cpxyz(txyz2[j], txyz[ik]);
        }
        *tnum += s2idx[i][2];
    }
}

void pair_checking(long ip, long ds, long num_residue, char *pdbfile, long *num_bp,
                   long **pair_num)
{
    long i, j;

    if (!ip) {  /* ideal base-pairing */
        if (ds == 1) {
            if (*num_bp > num_residue) {
                *num_bp = num_residue;
                fprintf(stderr, "processing all %ld residues\n\n", *num_bp);
            }
        } else if (num_residue % 2) {  /* ds = 2 */
            fprintf(stderr, "%s has odd number (%ld) of residues\n", pdbfile,
                    num_residue);
            fatal("Please specify base-pairing residue numbers\n");
        } else {
            i = num_residue / 2;
            if (*num_bp < i)
                fatal("Please specify base-pairing residue numbers\n");
            if (*num_bp > i) {
                *num_bp = i;
                fprintf(stderr, "processing all %ld base-pairs\n\n", i);
            }
        }

        for (i = 1; i <= *num_bp; i++) {
            pair_num[1][i] = i;
            if (ds == 2)
                pair_num[2][i] = 2 * *num_bp - i + 1;
        }
    } else {
        if (ds * *num_bp > num_residue)
            fprintf(stderr, "some residue has more than one pair\n");
        for (i = 1; i <= ds; i++)
            for (j = 1; j <= *num_bp; j++)
                if (pair_num[i][j] > num_residue) {
                    fprintf(stderr, "residue index %ld too big (> %ld)\n",
                            pair_num[i][j], num_residue);
                    fatal("please check your input file\n");
                }
    }
}

void drct_checking(long ds, long num_bp, long **pair_num, long **seidx, char **AtomName,
                   double **xyz, long *parallel, long *bbexist, long **o3p_brk, FILE * fp)
/* check for proper O3'-P linkage, strand I 5'--> 3' direction, parallel duplex */
{
    long i, ib, ie, ioffset, j, jp1, rnum, n = 0;
    long direction[7], **O3, **P;

    O3 = lmatrix(1, 2, 1, num_bp);
    P = lmatrix(1, 2, 1, num_bp);

    for (i = 1; i <= ds; i++)
        for (j = 1; j <= num_bp; j++) {
            rnum = pair_num[i][j];
            ib = seidx[rnum][1];
            ie = seidx[rnum][2];
            P[i][j] = find_1st_atom(" P  ", AtomName, ib, ie, "");
            O3[i][j] = find_1st_atom(" O3'", AtomName, ib, ie, "");
            if (P[i][j] && O3[i][j]) {
                ++*bbexist;
                if (within_limits(xyz[O3[i][j]], xyz[P[i][j]], 0, O3P_UPPER))
                    n++;
            }
        }

    if (n) {
        print_sep(fp, '*', 76);
        fprintf(stderr, "WARNING: %ld out of %ld bases have O3'[i] wrongly connected"
                " to P[i]\n", n, ds * num_bp);
        fprintf(fp, "WARNING: %ld out of %ld bases have O3'[i] wrongly connected"
                " to P[i]\n", n, ds * num_bp);
    }

    init_lvector(direction, 1, 6, 0);
    for (i = 1; i <= ds; i++) {
        ioffset = (i - 1) * 3;
        for (j = 1; j < num_bp; j++) {
            jp1 = j + 1;
            if (O3[i][j] && P[i][jp1] &&  /* j to j + 1 */
                within_limits(xyz[P[i][jp1]], xyz[O3[i][j]], 0, O3P_UPPER)) {
                ++direction[ioffset + 1];
                continue;
            }
            if (O3[i][jp1] && P[i][j] &&  /* j + 1 to j */
                within_limits(xyz[P[i][j]], xyz[O3[i][jp1]], 0, O3P_UPPER)) {
                ++direction[ioffset + 2];
                continue;
            }
            ++direction[ioffset + 3];
            if (*bbexist)
                o3p_brk[i][j] = 1;
        }
    }

    if ((direction[1] && direction[2]) || (direction[4] && direction[5])) {
        fprintf(stderr, "This structure contains intra-chain direction reverse\n");
        fprintf(fp, "This structure contains intra-chain direction reverse\n");
    } else {
        if (*bbexist && (direction[3] || direction[6])) {
            fprintf(stderr, "This structure has broken O3'[i] to P[i+1] linkages\n");
            fprintf(fp, "This structure has broken O3'[i] to P[i+1] linkages\n");
        }
        if (direction[2] && !direction[1]) {
            fprintf(stderr, "WARNING: Strand I in 3'-->5' direction!!!\n");
            fprintf(fp, "WARNING: Strand I in 3'-->5' direction!!!\n");
        }
        if (direction[4] && !direction[5]) {
            fprintf(stderr, "This is a parallel duplex\n");
            *parallel = 1;
        }
    }
    print_sep(fp, '*', 76);

    free_lmatrix(O3, 1, ds, 1, num_bp);
    free_lmatrix(P, 1, ds, 1, num_bp);
}

void base_frame(long num_residue, char *bseq, long **seidx, long *RY, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                char *BDIR, double **orien, double **org)
/* get the local reference frame for each base. only the ring atoms are
   included in least-squares fitting. similar to REF_FRAMES */
{
    static char *RingAtom[] = { RA_LIST };
    char idmsg[BUF512], sidmsg[BUF512], spdb[BUF512];
    char *sChainID, **sAtomName, **sResName, **sMiscs;

    double orgi[4];
    double **eRing_xyz, **fitted_xyz, **sRing_xyz, **sxyz, **R;

    long i, ib, ie, j, RingAtom_num;
    long exp_katom, nmatch, snum, std_katom;
    long *sResSeq;

    eRing_xyz = dmatrix(1, 9, 1, 3);
    sRing_xyz = dmatrix(1, 9, 1, 3);
    sAtomName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 4);
    sResName = cmatrix(1, NUM_RESIDUE_ATOMS, 0, 3);
    sChainID = cvector(1, NUM_RESIDUE_ATOMS);
    sResSeq = lvector(1, NUM_RESIDUE_ATOMS);
    sxyz = dmatrix(1, NUM_RESIDUE_ATOMS, 1, 3);
    sMiscs = cmatrix(1, NUM_RESIDUE_ATOMS, 0, NMISC);
    fitted_xyz = dmatrix(1, 9, 1, 3);
    R = dmatrix(1, 3, 1, 3);

    for (i = 1; i <= num_residue; i++) {
        if (RY[i] < 0)
            continue;  /* non-bases */

        ib = seidx[i][1];
        ie = seidx[i][2];
        sprintf(idmsg, ": residue name %s, chain %c, number %4ld%c",
                ResName[ib], ChainID[ib], ResSeq[ib], Miscs[ib][2]);
        RingAtom_num = (RY[i] == 1) ? 9 : 6;
        sprintf(spdb, "%sAtomic_%c.pdb", BDIR, bseq[i]);
        snum = read_pdb(spdb, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
                        sMiscs, 1, "*");
        sprintf(sidmsg, "in standard base: %s", spdb);

        nmatch = 0;
        for (j = 0; j < RingAtom_num; j++) {
            exp_katom = find_1st_atom(RingAtom[j], AtomName, ib, ie, idmsg);
            std_katom = find_1st_atom(RingAtom[j], sAtomName, 1, snum, sidmsg);
            if (exp_katom && std_katom) {
                ++nmatch;
                cpxyz(xyz[exp_katom], eRing_xyz[nmatch]);
                cpxyz(sxyz[std_katom], sRing_xyz[nmatch]);
            }
        }

        ls_fitting(sRing_xyz, eRing_xyz, nmatch, fitted_xyz, R, orgi);

        cpxyz(orgi, org[i]);
        mst2orien(orien[i], 0, R);
    }

    free_dmatrix(eRing_xyz, 1, 9, 1, 3);
    free_dmatrix(sRing_xyz, 1, 9, 1, 3);
    free_pdb(NUM_RESIDUE_ATOMS, NULL, sAtomName, sResName, sChainID, sResSeq, sxyz,
             sMiscs);
    free_dmatrix(fitted_xyz, 1, 9, 1, 3);
    free_dmatrix(R, 1, 3, 1, 3);
}

void base_str(char chain_id, long res_seq, char *misc, char *rname, char bcode,
              long stnd, char *idmsg)
/* get base name: all information to uniquely identify a base residue */
{
    char snum[10], rname_cp[10], modelNum[10], iCode = misc[2];
    long i, nlen = 4;

    strncpy(modelNum, misc + 30, nlen);
    modelNum[nlen] = '\0';
    for (i = 0; i < nlen; i++)
        if (modelNum[i] == '0')
            modelNum[i] = '.';
    sprintf(snum, "%4ld", res_seq);
    for (i = 0; i < 4; i++)
        if (snum[i] == ' ')
            snum[i] = '.';
    if (chain_id == ' ')
        chain_id = '-';
    if (iCode == ' ')
        iCode = '_';
    strcpy(rname_cp, rname);
    for (i = 0; i < 3; i++)
        if (rname_cp[i] == ' ')
            rname_cp[i] = '.';
    if (stnd == 1)  /* strand I */
        sprintf(idmsg, "%4s>%c:%4s%c:[%s]%c", modelNum, chain_id, snum, iCode,
                rname_cp, bcode);
    else  /* strand II */
        sprintf(idmsg, "%c[%s]:%4s%c:%c<%4s", bcode, rname_cp, snum, iCode,
                chain_id, modelNum);
}

void write_lkglist(long nbond, long **linkage, char **AtomName, char **ResName,
                   char *ChainID, long *ResSeq, char **Miscs)
/* write out a bond list for checking and easy modification */
{
    char b1[BUF512], b2[BUF512];
    long i, ia, ib;
    FILE *fp;

    fp = open_file(LKG_FILE, "w");
    for (i = 1; i <= nbond; i++) {
        ia = linkage[i][1];
        ib = linkage[i][2];
        base_str(ChainID[ia], ResSeq[ia], Miscs[ia], ResName[ia], '\0', 1, b1);
        base_str(ChainID[ib], ResSeq[ib], Miscs[ib], ResName[ib], '\0', 1, b2);
        fprintf(fp, "%5ld %5ld  # %5ld %s {%s} with %s {%s}\n", ia, ib, i,
                AtomName[ia], b1, AtomName[ib], b2);
    }
    close_file(fp);
}

void hbond_info(long **pair_num, char *bseq, long **seidx, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                long *RY, long *num_hbond, long **hb_linkage)
/* get all possible H-bond information with upper/lower layer information */
{
    char atom1[10], HB_ATOM[BUF512];
    double HB_UPPER[2];
    long i, j, k, ki, kj, num_bp = 2;
    long **num_list;
    FILE *fph;

    hb_crt_alt(HB_UPPER, HB_ATOM, atom1);

    num_list = lmatrix(1, BUF512, 1, 2);
    fph = open_file(HB_FILE, "w");
    for (k = 1; k <= num_bp; k++)
        for (i = 1; i < pair_num[k][0]; i++)
            for (j = i + 1; j <= pair_num[k][0]; j++) {  /* pairwise */
                ki = pair_num[k][i];
                kj = pair_num[k][j];
                if (RY[ki] < 0 || RY[kj] < 0)  /* ki could be moved out */
                    continue;  /* got to be a base residue: bd0053/1J8L */
                hbond_list(ki, kj, AtomName, ResName, ChainID, ResSeq, xyz, Miscs,
                           bseq, seidx, hb_linkage, HB_UPPER, HB_ATOM, num_list,
                           num_hbond, k, fph);
            }
    free_lmatrix(num_list, 1, BUF512, 1, 2);
    close_file(fph);
}

void hbond_pdb(long num, long num_residue, char *bseq, long **seidx, long *RY,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, long *num_hbond, long **hb_linkage, long pwise)
/* get H-bonds betwen base-pairs for a general PDB data file */
{
    char atom1[10], HB_ATOM[BUF512];
    double BPRS[7], **orien, **org, **NC1xyz, **o3_p;
    double HB_UPPER[2], rtn_val[RTNNUM];
    long bpid, i, j, *idx, **num_list, **ring_atom;
    FILE *fph;

    orien = dmatrix(1, num_residue, 1, 9);
    org = dmatrix(1, num_residue, 1, 3);
    NC1xyz = dmatrix(1, num_residue, 1, 7);  /* RN9/YN1 & C1' atomic coordinates */
    o3_p = dmatrix(1, num_residue, 1, 8);  /* O3'/P atomic coordinates */

    base_info(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID, ResSeq, Miscs,
              xyz, orien, org, NC1xyz, o3_p, BPRS);

    /* 1-9 ring atom index, 10 # of ring atoms, 11-19 first level */
    ring_atom = lmatrix(1, num_residue, 1, 19);
    ring_oidx(num, num_residue, RY, seidx, AtomName, xyz, ring_atom);

    hb_crt_alt(HB_UPPER, HB_ATOM, atom1);
    num_list = lmatrix(1, BUF512, 1, 2);

    idx = lvector(1, num);
    atom_idx(num, AtomName, NULL, idx);

    fph = open_file(HB_FILE, "w");
    for (i = 1; i < num_residue; i++) {
        if (!pwise && RY[i] < 0)
            continue;
        for (j = i + 1; j <= num_residue; j++) {
            if (!pwise) {  /* check for base pairing */
                if (RY[j] < 0)  /* only considering hydration for the 1st layer */
                    continue;
                check_pair(i, j, bseq, seidx, xyz, NC1xyz, orien, org, idx,
                           AtomName, BPRS, rtn_val, &bpid, ring_atom, 0);
                if (!bpid)
                    continue;  /* not a base-pair */
            }
            hbond_list(i, j, AtomName, ResName, ChainID, ResSeq, xyz, Miscs, bseq, seidx,
                       hb_linkage, HB_UPPER, HB_ATOM, num_list, num_hbond, 1, fph);
        }
    }

    free_dmatrix(orien, 1, num_residue, 1, 9);
    free_dmatrix(org, 1, num_residue, 1, 3);
    free_dmatrix(NC1xyz, 1, num_residue, 1, 7);
    free_dmatrix(o3_p, 1, num_residue, 1, 8);
    free_lmatrix(ring_atom, 1, num_residue, 1, 19);
    free_lmatrix(num_list, 1, BUF512, 1, 2);
    free_lvector(idx, 1, num);
    close_file(fph);
}

void hbond_list(long i, long j, char **AtomName, char **ResName, char *ChainID,
                long *ResSeq, double **xyz, char **Miscs, char *bseq, long **seidx,
                long **hb_linkage, double *HB_UPPER, char *HB_ATOM, long **num_list,
                long *num_hbond, long ilayer, FILE * fph)
/* write out h-bond listing information: <num_hbond> already a pointer */
{
    char b1[BUF512], b2[BUF512];
    long ir, jr, k, num_hb;

    ir = seidx[i][1];
    jr = seidx[j][1];
    base_str(ChainID[ir], ResSeq[ir], Miscs[ir], ResName[ir], bseq[i], 1, b1);
    base_str(ChainID[jr], ResSeq[jr], Miscs[jr], ResName[jr], bseq[j], 1, b2);
    hb_numlist(i, j, bseq[i], bseq[j], seidx, AtomName, xyz, HB_UPPER, HB_ATOM,
               &num_hb, num_list);
    for (k = 1; k <= num_hb; k++) {
        hb_linkage[++*num_hbond][1] = ilayer;
        hb_linkage[*num_hbond][2] = num_list[k][1];
        hb_linkage[*num_hbond][3] = num_list[k][2];
        fprintf(fph, "%5ld %5ld  # %5ld", num_list[k][1], num_list[k][2], *num_hbond);
        fprintf(fph, "%s {%s} with %s {%s}\n", AtomName[num_list[k][1]], b1,
                AtomName[num_list[k][2]], b2);
    }
}

void hb_numlist(long i, long j, char basei, char basej, long **seidx, char **AtomName,
                double **xyz, double *HB_UPPER, char *HB_ATOM, long *num_hb,
                long **num_list)
/* determine atom # list for each base-pair from HB_INFO */
{
    char atom1[10], atom2[10], HB_INFO[BUF512], *pchar;
    long k, num;

    get_hbond_ij(i, j, basei, basej, HB_UPPER, seidx, AtomName, HB_ATOM, xyz, HB_INFO);

    if (sscanf(HB_INFO, "[%ld]", &num) != 1)
        num = 0;
    if (num) {
        pchar = strchr(HB_INFO, ' ');
        for (k = 1; k <= num; k++) {
            strncpy(atom1, pchar + 1, 4);
            atom1[4] = '\0';
            strncpy(atom2, pchar + 6, 4);
            atom2[4] = '\0';
            pchar += 15;  /* point to next H-bond */
            num_list[k][1] = find_1st_atom(atom1, AtomName, seidx[i][1], seidx[i][2], "");
            num_list[k][2] = find_1st_atom(atom2, AtomName, seidx[j][1], seidx[j][2], "");
        }
    }
    *num_hb = num;
}

void hb_information(long num_bp, long **pair_num, char **bp_seq, long **seidx,
                    char **AtomName, double **xyz, long *WC_info, FILE * fp)
/* get detailed H-bonding information for each pair */
{
    char HB_INFO[BUF512], HB_ATOM[BUF512], temp[10];
    double HB_UPPER[2];
    long i, j, k;

    hb_crt_alt(HB_UPPER, HB_ATOM, temp);

    print_sep(fp, '*', 76);
    fprintf(fp, "Detailed H-bond information: atom-name pair and length" " [%s]\n",
            HB_ATOM);

    for (k = 1; k <= num_bp; k++) {
        i = pair_num[1][k];
        j = pair_num[2][k];
        get_hbond_ij(i, j, bp_seq[1][k], bp_seq[2][k], HB_UPPER, seidx, AtomName,
                     HB_ATOM, xyz, HB_INFO);
        fprintf(fp, "%4ld %c-%c%c%c-%c  %s\n", k, bp_seq[1][k],
                (WC_info[k] == 2) ? '-' : '*', WC_info[k] ? '-' : '*',
                bp_seq[0][k], bp_seq[2][k], HB_INFO);
    }
}

long good_hbatoms(char *HB_ATOM, char *atom1, char *atom2)
/* check if two atom pair could form a H-bond */
{
    static char *PO[] = { " O1P", " O2P", " O3'", " O4'", " O5'", " N7 " };
    static long numPO = sizeof PO / sizeof PO[0] - 1;

    if (num_strmatch(atom1, PO, 0, numPO) && num_strmatch(atom2, PO, 0, numPO))
        return 0;  /* no H-bond within PO4 group */
    else {
        if ((strchr("ON", atom1[1]) || strchr("ON", atom2[1]))
            /* at least one O/N atom */
            && (strchr(HB_ATOM, '*') ||
                (strchr(HB_ATOM, atom1[1]) && strchr(HB_ATOM, atom2[1]))))
            return 1;
        else
            return 0;
    }
}

void read_lkginfo(char *lkgfile, long num, long *nbond, long **linkage)
/* read in bond linkage information from "lkgfile": r3d_atom & stack2img generate
 * a file called "bonds_lkg.dat", which can be modified to fit any purpose
 */
{
    char str[BUF512];
    long ia, ib;
    FILE *fp;

    fp = open_file(lkgfile, "r");
    while (fgets(str, sizeof str, fp) != NULL)
        if (sscanf(str, "%ld %ld", &ia, &ib) == 2) {
            if (!lval_in_range(ia, 1, num) || !lval_in_range(ib, 1, num)) {
                fprintf(stderr, "%s", str);
                fprintf(stderr, "         has atom serial number out of range\n");
                continue;
            }
            linkage[++*nbond][1] = ia;
            linkage[*nbond][2] = ib;
        }
    close_file(fp);
}

void read_hbinfo(char *hbfile, long num, long *num_hbond, long **hb_linkage)
/* read in H-bonding information from "hbfile": r3d_atom & stack2img generate
 * a file called "hbonds_info.dat", which can be modified to fit any purpose
 */
{
    char str[BUF512];
    long ia, ib;
    FILE *fp;

    fp = open_file(hbfile, "r");
    while (fgets(str, sizeof str, fp) != NULL)
        if (sscanf(str, "%ld %ld", &ia, &ib) == 2) {
            if (!lval_in_range(ia, 1, num) || !lval_in_range(ib, 1, num)) {
                fprintf(stderr, "[%s] has atom serial number out of range\n", str);
                continue;
            }
            hb_linkage[++*num_hbond][1] = 1;
            hb_linkage[*num_hbond][2] = ia;
            hb_linkage[*num_hbond][3] = ib;
        }
    close_file(fp);
}

void hb_atompair(long num_hbonds, char **hb_atom1, char **hb_atom2, double *hb_dist,
                 double *HB_UPPER)
/* get the H-bonding atom-pairs */
{
    double dtmp[3];
    long ddidx[3], *matched_idx, **idx2;
    long k, m = 0, n, num_iter = 1;

    if (!num_hbonds)
        return;

    matched_idx = lvector(1, num_hbonds);

    while (1) {
        if (matched_idx[num_iter]) {
            num_iter++;
            continue;
        }
        for (k = 1; k <= 2; k++) {
            dtmp[k] = hb_dist[num_iter];
            ddidx[k] = num_iter;
        }
        for (n = 1; n <= num_hbonds; n++) {
            if (n == num_iter || matched_idx[n])
                continue;
            if (!strcmp(hb_atom1[n], hb_atom1[num_iter]) && hb_dist[n] < dtmp[1]) {
                dtmp[1] = hb_dist[n];
                ddidx[1] = n;
            }
            if (!strcmp(hb_atom2[n], hb_atom2[num_iter]) && hb_dist[n] < dtmp[2]) {
                dtmp[2] = hb_dist[n];
                ddidx[2] = n;
            }
        }
        if (ddidx[1] == ddidx[2]) {
            hb_dist[ddidx[1]] = -hb_dist[ddidx[1]];
            num_iter = 1;  /* reset the iterator */
            for (n = 1; n <= num_hbonds; n++) {
                if (matched_idx[n])
                    continue;
                if (!strcmp(hb_atom1[n], hb_atom1[ddidx[1]]) ||
                    !strcmp(hb_atom2[n], hb_atom2[ddidx[2]])) {
                    matched_idx[n] = 1;
                    m++;
                }
            }
            if (m >= num_hbonds)
                break;
        } else
            num_iter++;
    }

    /* === further processing by adding more possible H-bonds === */
    idx2 = lmatrix(1, num_hbonds, 1, 2);

    for (k = 1; k <= num_hbonds; k++) {
        if (hb_dist[k] > 0.0)
            continue;
        idx2[k][1] = 9;
        idx2[k][2] = 9;
        for (m = 1; m <= num_hbonds; m++) {
            if (m == k || hb_dist[m] < 0.0)
                continue;
            if (!strcmp(hb_atom1[m], hb_atom1[k]))
                idx2[m][1] = 1;
            if (!strcmp(hb_atom2[m], hb_atom2[k]))
                idx2[m][2] = 1;
        }
    }

    /* Note that at this point, there are only 3 possibilities:
       9 + 9 = 18 means already located H-bonds;
       1 + 1 = 2 means the two atoms in two separately H-bonds;
       1 + 0 = 0 + 1 = 1 means only one of atoms in H-bonds */

    for (k = 1; k <= num_hbonds; k++)
        if (idx2[k][1] + idx2[k][2] == 1 && hb_dist[k] <= HB_UPPER[1])
            hb_dist[k] = -hb_dist[k];

    free_lmatrix(idx2, 1, num_hbonds, 1, 2);
    free_lvector(matched_idx, 1, num_hbonds);

    /* === End of further processing by adding more possible H-bonds === */
}

void get_hbond_ij(long i, long j, char basei, char basej, double *HB_UPPER, long **seidx,
                  char **AtomName, char *HB_ATOM, double **xyz, char *HB_INFO)
/* get H-bond length information between residue i and j */
{
    char donacc, **hb_atom1, **hb_atom2, stmp[20];
    double *hb_dist;
    double hbhigh = HB_UPPER[0], hblow = 2.0;
    long k, m, n, num_hbonds = 0;

    hb_atom1 = cmatrix(1, BUF512, 0, 4);
    hb_atom2 = cmatrix(1, BUF512, 0, 4);
    hb_dist = dvector(1, BUF512);

    for (m = seidx[i][1]; m <= seidx[i][2]; m++) {
        for (n = seidx[j][1]; n <= seidx[j][2]; n++) {
            if (good_hbatoms(HB_ATOM, AtomName[m], AtomName[n]) &&
                within_limits(xyz[n], xyz[m], hblow, hbhigh)) {
                if (++num_hbonds > BUF512)
                    fatal("Too many possible H-bonds between two bases\n");
                strcpy(hb_atom1[num_hbonds], AtomName[m]);
                strcpy(hb_atom2[num_hbonds], AtomName[n]);
                hb_dist[num_hbonds] = p1p2_dist(xyz[n], xyz[m]);
            }
        }
    }

    if (num_hbonds)
        hb_atompair(num_hbonds, hb_atom1, hb_atom2, hb_dist, HB_UPPER);

    m = 0;
    for (k = 1; k <= num_hbonds; k++)
        if (hb_dist[k] < 0.0)
            m++;
    sprintf(HB_INFO, "[%ld]", m);
    for (k = 1; k <= num_hbonds; k++)
        if (hb_dist[k] < 0.0) {
            donacc = donor_acceptor(basei, basej, hb_atom1[k], hb_atom2[k]);
            sprintf(stmp, " %s%c%s %4.2f", hb_atom1[k], donacc, hb_atom2[k], -hb_dist[k]);
            strcat(HB_INFO, stmp);
        }

    free_cmatrix(hb_atom1, 1, BUF512, 0, 4);
    free_cmatrix(hb_atom2, 1, BUF512, 0, 4);
    free_dvector(hb_dist, 1, BUF512);
}

char donor_acceptor(char basei, char basej, char *hb_atom1, char *hb_atom2)
/* Based on canonical form of ACGITU, decide if a HB-pair is donor-acceptor */
{
    char da[3], hbatom_type = '*', *pchar;
    char ia = '\0', ja = '\0';
    /* _A: acceptor; _D: donor; _X: both acceptor & donor; _?: not sure */
    static char *cmn_base = CB_LIST;
    static char *da_type[7] = { "AD", "AX", "XD", "XX", "DA", "DX", "XA" };
    static char *bb_da[6] = { " O1P_A", " O2P_A", " O5'_A", " O4'_A", " O3'_A",
        " O2'_X"
    };
    static char *base_da[6][6] = {
        {" N9 _?", " N7 _A", " N6 _D", " N1 _A", " N3 _A"},  /* A */
        {" N1 _?", " O2 _A", " N3 _A", " N4 _D"},  /* C */
        {" N9 _?", " N7 _A", " O6 _A", " N1 _D", " N2 _D", " N3 _A"},  /* G */
        {" N9 _?", " N7 _A", " O6 _A", " N1 _D", " N3 _A"},  /* I */
        {" N1 _?", " O2 _A", " N3 _D", " O4 _A"},  /* T */
        {" N1 _?", " O2 _A", " N3 _D", " O4 _A"}  /* U */
    };
    long i, inum = -1, jnum = -1, num = 6;

    if ((pchar = strchr(cmn_base, toupper(basei))) != NULL)
        inum = pchar - cmn_base;
    if ((pchar = strchr(cmn_base, toupper(basej))) != NULL)
        jnum = pchar - cmn_base;
    if (inum >= 0 && jnum >= 0) {  /* both are one of the six bases */
        for (i = 0; i < num; i++) {  /* check if backbone atom */
            if (!strncmp(bb_da[i], hb_atom1, 4))
                ia = bb_da[i][5];
            if (!strncmp(bb_da[i], hb_atom2, 4))
                ja = bb_da[i][5];
        }
        if (!ia)
            for (i = 0; i < num; i++)
                if (base_da[inum][i] && !strncmp(base_da[inum][i], hb_atom1, 4)) {
                    ia = *(base_da[inum][i] + 5);
                    break;
                }
        if (!ja)
            for (i = 0; i < num; i++)
                if (base_da[jnum][i] && !strncmp(base_da[jnum][i], hb_atom2, 4)) {
                    ja = *(base_da[jnum][i] + 5);
                    break;
                }
        if (ia && ja) {
            sprintf(da, "%c%c", ia, ja);
            if (num_strmatch(da, da_type, 0, 6))
                hbatom_type = '-';
        }
    }

    return hbatom_type;
}

void hb_crt_alt(double *HB_UPPER, char *HB_ATOM, char *ALT_LIST)
/* read in H-bond length upper limit etc from <misc_3dna.par> */
{
    char BDIR[BUF512], str[BUF512];
    FILE *fp;

    /* read in H-bond length upper limit */
    get_BDIR(BDIR, PAR_FILE, 1);
    strcat(BDIR, PAR_FILE);
    fprintf(stderr, " ...... reading file: %s ...... \n", PAR_FILE);

    fp = open_file(BDIR, "r");
    if ((fgets(str, sizeof str, fp) == NULL) ||
        (sscanf(str, "%lf %lf %s %s", &HB_UPPER[0], &HB_UPPER[1], HB_ATOM,
                ALT_LIST) != 4))
        fatal("error reading upper HB criteria/alternative indicator\n");
    upperstr(HB_ATOM);
    upperstr(ALT_LIST);
    strcat(ALT_LIST, " ");  /* attach space */

    close_file(fp);
}

void atom_info(long idx, char atoms_list[NELE][3], double *covalence_radii,
               double *vdw_radii)
{
    static char *ALIST[NELE] =
        { "XX", " C", " O", " H", " N", " S", " P", " F", "CL", "BR", " I",
        "SI"
    };
    static double CRADII[NELE] =
        { 1.666, 0.762, 0.646, 0.352, 0.689, 1.105, 1.000, 0.619, 1.022,
        1.183, 1.378, 1.105  /* XX original 1.200 */
    };
    static double VRADII[NELE] =
        { 2.00, 1.70, 1.52, 1.20, 1.55, 1.80, 1.80, 1.47, 1.75, 1.85, 1.98,
        2.10
    };
    long i;

    if (idx == 1)
        for (i = 0; i < NELE; i++)
            strcpy(atoms_list[i], ALIST[i]);
    else if (idx == 2)
        for (i = 0; i < NELE; i++)
            covalence_radii[i] = CRADII[i];
    else if (idx == 3)
        for (i = 0; i < NELE; i++)
            vdw_radii[i] = VRADII[i];
    else
        fatal("wrong options for <atom_info>: should be 1, 2 or 3\n");
}

void aname2asym(char *aname0, char *my_asym)
/* get 2-letter atomic symbol from atom name in PDB file */
{
    char aname[5];
    long j;

    strcpy(aname, aname0);  /* make a copy for later reference */
    strcpy(my_asym, "XX");
    for (j = 0; j <= 3; j++)
        if (!isalpha((int) aname[j]))
            aname[j] = ' ';  /* no alphabets changed to ' ' */
    if (aname[0] != ' ' && aname[1] != ' ')  /* done */
        for (j = 0; j <= 1; j++)
            my_asym[j] = aname[j];
    else if (aname[0] != ' ' && aname[1] == ' ') {
        my_asym[0] = ' ';
        my_asym[1] = aname[0];  /* H2'1 or O1P' in RIA64 of trna12 */
    } else if (aname[0] == ' ' && aname[1] != ' ') {  /* most common */
        my_asym[0] = ' ';
        my_asym[1] = aname[1];  /*  */
    } else  /* unusual situation */
        fprintf(stderr, "Atomic name <%s>: anything wrong?\n", aname0);
    my_asym[2] = '\0';
}

void atom_idx(long num, char **AtomName, char **Miscs, long *idx)
/* get atom index for calculating bond linkage */
{
    static char *atom_names[] = {
        " H", "HE",
        "LI", "BE", " B", " C", " N", " O", " F", "NE",
        "NA", "MG", "AL", "SI", " P", " S", "CL", "AR",
        " K", "CA", "SC", "TI", " V", "CR", "MN", "FE", "CO",
        "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
        "RB", "SR", " Y", "ZR", "NB", "MO", "TC", "RU", "RH",
        "PD", "AG", "CD", "IN", "SN", "SB", "TE", " I", "XE",
        "CS", "BA", "LA", "HF", "TA", " W", "RE", "OS", "IR",
        "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",
        "FR", "RA", "AC", "RF", "DB", "SG", "BH", "HS", "MT",
        "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB",
        "DY", "HO", "ER", "TM", "YB", "LU",
        "TH", "PA", " U", "NP", "PU", "AM", "CM", "BK",
        "CF", "ES", "FM", "MD", "NO", "LR"
    };
    char pdb_asym[3], my_asym[3], atoms_list[NELE][3];
    long i, j, k, num_atoms;

    num_atoms = sizeof atom_names / sizeof atom_names[0] - 1;  /* minus - 1 */
    atom_info(1, atoms_list, NULL, NULL);

    for (i = 1; i <= num; i++) {
        /* atomic symbol from PDB data file */
        strcpy(pdb_asym, "XX");
        k = 0;  /* atomic symbol */
        if (Miscs && strlen(Miscs[i]) >= 27) {  /* with this info. */
            strncpy(pdb_asym, Miscs[i] + 25, 2);
            pdb_asym[2] = '\0';
            if (num_strmatch(pdb_asym, atom_names, 0, num_atoms))
                k = 1;
            else
                fprintf(stderr, "Columns #77-78 <%s>: illegal atomic symbol\n", pdb_asym);
        }

        /* atomic symbol deduced from atom name */
        aname2asym(AtomName[i], my_asym);
        if (num_strmatch(my_asym, atom_names, 0, num_atoms)) {
            if (k && strcmp(pdb_asym, my_asym))
                fprintf(stderr, "Atomic symbol <%s> (#77-78) and <%s> are not equal!\n",
                        pdb_asym, my_asym);
        } else {
            fprintf(stderr, "Atom name <%s> converted to <%s>: illegal atomic symbol\n",
                    AtomName[i], my_asym);
            fprintf(stderr, "Maybe your PDB file does not conform to the standard\n");
            if (k)  /* use that extracted from PDB file */
                strcpy(my_asym, pdb_asym);
        }

        for (j = 0; j < NELE; j++)
            if (!strcmp(my_asym, atoms_list[j])) {
                idx[i] = j;
                break;
            }
        if (j >= NELE)
            idx[i] = 0;  /* not found: use default */
    }

}

void get_bonds(long num, char **AtomName, double **xyz, long num_residue, long *RY,
               long **seidx, long **connect)
/* get linkage information based on inter-atomic distance */
{
    long i, nbond, nbond_estimated, *idx, **linkage;

    idx = lvector(1, num);
    atom_idx(num, AtomName, NULL, idx);

    nbond_estimated = lround(NBOND_FNUM * NUM_RESIDUE_ATOMS);
    linkage = lmatrix(1, nbond_estimated, 1, 2);
    for (i = 1; i <= num_residue; i++)  /* linkage within each residue */
        if (RY[i] >= 0) {  /* a base residue */
            nbond = 0;
            atom_linkage(seidx[i][1], seidx[i][2], idx, xyz, NULL, NULL,
                         nbond_estimated, &nbond, linkage);
            lkg2connect(seidx[i][1], seidx[i][2], nbond, linkage, connect);
        }

    free_lvector(idx, 1, num);
    free_lmatrix(linkage, 1, nbond_estimated, 1, 2);
}

void atom_linkage(long ib, long ie, long *idx, double **xyz, char **Miscs, char *ChainID,
                  long nbond_estimated, long *nbond, long **linkage)
/* get atom linkage using covalent radii criterion. alt. position & chain info. */
{
    char a1, a2, c1, c2;
    double dst, covalence_radii[NELE];
    long j, k;

    /* Bond criteria: 1.15 * (rA + rB)
       RasMol: 0.56 + (rA + rB)   MacroModel: 1.25 * (rA + rB) */

    atom_info(2, NULL, covalence_radii, NULL);
    for (j = ib; j <= ie - 1; j++) {
        a1 = (Miscs == NULL) ? ' ' : Miscs[j][1];  /* atom alternative position */
        c1 = (ChainID == NULL) ? '_' : ChainID[j];
        for (k = j + 1; k <= ie; k++) {
            a2 = (Miscs == NULL) ? ' ' : Miscs[k][1];  /* atom alternative position */
            c2 = (ChainID == NULL) ? '_' : ChainID[k];
            if (a1 != ' ' && a2 != ' ' && a1 != a2)  /* same model: alt. position */
                continue;
            if (c1 != c2)  /* same chain ID */
                continue;
            dst = BOND_FACTOR * (covalence_radii[idx[j]] + covalence_radii[idx[k]]);
            if (within_limits(xyz[k], xyz[j], 0, dst)) {
                if (++*nbond > nbond_estimated)
                    fatal("too many linkages\n");
                else {
                    linkage[*nbond][1] = j;
                    linkage[*nbond][2] = k;
                }
            }
        }
    }
}

void lkg2connect(long ib, long ie, long nbond, long **linkage, long **connect)
/* from bond linkage to PDB connection table: SLOW for large structure rr0033 */
{
    long i, ilink, j;
    long idx[7];

    for (i = ib; i <= ie; i++) {
        ilink = 0;
        for (j = 1; j <= nbond; j++) {
            if (i == linkage[j][1]) {
                if (++ilink > 6) {
                    fprintf(stderr, "atom %ld has over 6 bonds\n", i);
                    break;
                } else
                    idx[ilink] = linkage[j][2];
            }
            if (i == linkage[j][2]) {
                if (++ilink > 6) {
                    fprintf(stderr, "atom %ld has over 6 bonds\n", i);
                    break;
                } else
                    idx[ilink] = linkage[j][1];
            }
        }
        if (ilink > 6)
            ilink = 6;  /* maximum six bonds */
        if (ilink > 0) {
            for (j = 1; j <= ilink; j++)
                connect[i][j] = idx[j];
            connect[i][7] = ilink;
            lsort(ilink, connect[i], idx);  /* ignore idx */
        }
    }
}

void init_htm_water(long waters, long num, long num_residue, long *idx, long **htm_water)
/* initialize matrix htm_water as follows:
 * row#1: c0 for # of atoms; c1-num for sequential residue # in the PDB file
 * row#2: c0 for # of residues; c1-num for connected base sequential residue #
 * row#3: c0 for hydration indication; c1-num for atom index */
{
    long i;

    htm_water[1][0] = num;  /* total # of atoms */
    htm_water[2][0] = num_residue;  /* total # of residue */
    htm_water[3][0] = waters;  /* if to check for hydration */

    for (i = 1; i <= num; i++)  /* make a copy of atom index */
        htm_water[3][i] = idx[i];
}

void identify_htw(long num_residue, long **seidx, long *RY, char **AtomName,
                  char **ResName, char *ChainID, long *ResSeq, char **Miscs,
                  double **xyz, long **htm_water)
/* make each HETATM record, including water, to its connected base residue */
{
    static char *WATER[] = WATER_LIST;
    char a1, a2, c1, c2;
    double dst, covalence_radii[NELE];
    long num_wat = sizeof WATER / sizeof WATER[0] - 1;
    long i, ib, ie, id, j, k = 0, m, num_H2O = 0;

    /* identify ligands and water molecules */
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] >= 0)  /* skip base residues */
            continue;
        ib = seidx[i][1];
        ie = seidx[i][2];
        id = ie - ib + 1;
        for (j = ib; j <= ie; j++)
            if (Miscs[j][0] == 'H') {
                htm_water[1][j] = i;
                htm_water[2][j] = -8;
            }

        if (num_strmatch(ResName[ib], WATER, 0, num_wat) || (id == 1 &&  /* one atom in this residue */
                                                             (!strcmp
                                                              (AtomName[ib], " O  ")
                                                              || !strcmp(AtomName[ib],
                                                                         " OW ")))) {
            for (j = ib; j <= ie; j++) {  /* possibly with Hs, e.g. NMR */
                htm_water[1][j] = i;
                htm_water[2][j] = -1;
            }
            num_H2O++;
            htm_water[4][num_H2O] = i;
        }
    }
    htm_water[4][0] = num_H2O;

    /* This section checks for covalent connection of a non-base residue,
     * e.g., drug molecule, to a base residue. They have to be on the same
     * alternative position and chain ID. Pure H2Os do not count here. */
    atom_info(2, NULL, covalence_radii, NULL);  /* cf. <atom_linkage> */
    for (i = 1; i <= num_residue; i++) {
        if (RY[i] >= 0)  /* skip base residues */
            continue;
        id = 0;  /* check if residue "i" is connected to a base residue */
        for (j = seidx[i][1]; j <= seidx[i][2]; j++) {
            /* assuming NO more than 1 connection to base; only HETATM/H2O */
            if (htm_water[2][j] >= 0)  /* -8 or -1 */
                break;
            a1 = Miscs[j][1];
            c1 = ChainID[j];
            for (k = 1; k <= num_residue; k++) {
                if (RY[k] < 0)  /* check only for connections to base residue */
                    continue;
                for (m = seidx[k][1]; m <= seidx[k][2]; m++) {
                    a2 = Miscs[m][1];
                    c2 = ChainID[m];
                    if (a1 != ' ' && a2 != ' ' && a1 != a2)
                        continue;  /* same model: alt. position */
                    if (c1 != c2)  /* same chain ID */
                        break;  /* restricted to residue m */
                    dst = BOND_FACTOR * (covalence_radii[htm_water[3][j]] +
                                         covalence_radii[htm_water[3][m]]);
                    if (within_limits(xyz[m], xyz[j], 0, dst)) {
                        id = 1;
                        goto ID_CNCT;
                    }
                }
            }
        }

      ID_CNCT:
        if (id) {  /* residue i connected to base residue k */
            for (j = seidx[i][1]; j <= seidx[i][2]; j++)
                htm_water[2][j] = k;
            ib = seidx[i][1];  /* HETATM residue */
            ie = seidx[k][1];  /* contacted base residue */
            if (ResSeq[ib] != ResSeq[ie] || Miscs[ib][2] != Miscs[ie][2])
                fprintf(stderr, "%ld(%ld,%ld) and %ld(%ld,%ld) have different residue"
                        " # / insertion code!\n", ResSeq[ib], ib, i, ResSeq[ie], ie, k);
        }
    }
}

long attached_residues(long inum_base, long *ivec, long *ivec2, long **seidx,
                       double **xyz, long **htm_water)
/* attach covalently connected residues as in +A etc */
{
    char str[BUF512];
    double cmb_par, H2O_HBDIST;
    long atype;  /* atoms to be included in hydration calculation */
    long i, iw, ir, j, k, m, n, tnum_res = inum_base;
    long num_residue = htm_water[2][0], num_H2O = htm_water[4][0];

    /* check if there is a covalent bond linkage between hetero
     * group (e.g., a drug molecule) with the base residue. */
    for (i = 1; i <= inum_base; i++) {
        ivec2[i] = ivec[i];
        for (j = 1; j <= num_residue; j++) {
            k = seidx[j][1];
            if (htm_water[2][k] == ivec[i])
                ivec2[++tnum_res] = htm_water[1][k];
        }
    }

    if (!htm_water[3][0])
        return tnum_res;

    cmb_par = get_rowpar(8);
    sprintf(str, "%.0f", cmb_par);
    atype = str[strlen(str) - 1] - '0';
    H2O_HBDIST = (cmb_par - atype) / MFACTOR;

    inum_base = tnum_res;
    n = inum_base + 1;  /* starting index for hydration residues */
    for (ir = 1; ir <= num_H2O; ir++) {  /* time-consuming part: best way! */
        iw = htm_water[4][ir];  /* hydration index */
        for (m = seidx[iw][1]; m <= seidx[iw][2]; m++) {
            if (htm_water[3][m] != 2)  /* not O */
                continue;
            for (i = 1; i <= inum_base; i++) {
                k = ivec2[i];
                for (j = seidx[k][1]; j <= seidx[k][2]; j++) {
                    if (atype == 0 &&  /* only for N and O */
                        htm_water[3][j] != 2 && htm_water[3][j] != 4)
                        continue;
                    if (atype == 1 &&  /* N, O and C */
                        htm_water[3][j] != 2 && htm_water[3][j] != 4 && htm_water[3][j] != 1)  /* following Pascal Auffinger */
                        continue;
                    if (within_limits(xyz[m], xyz[j], 0, H2O_HBDIST) &&
                        !lval_in_set(iw, n, tnum_res, ivec2)) {
                        ivec2[++tnum_res] = iw;
                    }
                }
            }
        }
    }

    if (tnum_res > inum_base + 1) {  /* sort the hydration residues into order */
        long *idx;

        k = tnum_res - inum_base;
        idx = lvector(1, k);
        lsort(k, ivec2 + inum_base, idx);

        free_lvector(idx, 1, k);
    }

    return tnum_res;
}

void print_pairinfo(long i, long j, char basei, char basej, double *rtn_val, double *chi,
                    char *HB_ATOM, double *HB_UPPER, long **seidx, char **AtomName,
                    double **xyz, char *bseq, long detailed, FILE * fp)
/* print out base-pairing information */
{
    char antip, bptype[4], HB_INFO[BUF512];
    long k;

    antip = (rtn_val[35] < 0.0) ? '-' : '+';
    sprintf(bptype, "%c%c%c", bseq[i], antip, bseq[j]);
    fprintf(fp, "             %s  ", bptype);
    for (k = 27; k <= 32; k++)  /* bp parameters */
        fprintf(fp, "%8.2f", rtn_val[k]);
    fprintf(fp, "\n");
    get_hbond_ij(i, j, basei, basej, HB_UPPER, seidx, AtomName, HB_ATOM, xyz, HB_INFO);
    fprintf(fp, "             %s %s\n", bptype, HB_INFO);

    if (!detailed)
        return;

    fprintf(fp, "                    %s %s %s\n", rtn_val[35] < 0.0 ? "anti-parallel" : "parallel",  /* z-axis */
            rtn_val[33] < 0.0 ? "trans" : "cis",  /* x-axis */
            rtn_val[36] < 0.0 ? "trans" : "cis");  /* C1-N vector */
    fprintf(fp, "            ");
    for (k = 33; k <= 35; k++)
        fprintf(fp, "%c", rtn_val[k] > 0.0 ? '+' : '-');
    for (k = 33; k <= 36; k++)
        fprintf(fp, "%9.1f", dot2ang(rtn_val[k]));
    fprintf(fp, "%9.1f%9.1f\n", chi[i], chi[j]);

    fprintf(fp, "     ");
    for (k = 1; k <= 8; k++)
        fprintf(fp, "%8.2f", rtn_val[k]);
    fprintf(fp, "\n");
}

void check_pair(long i, long j, char *bseq, long **seidx, double **xyz, double **NC1xyz,
                double **orien, double **org, long *idx, char **AtomName, double *BPRS,
                double *rtn_val, long *bpid, long **ring_atom, long network)
/* Checking if two bases form a pair according to several criteria
 * rtn_val[RTNNUM]:
 *       d, dv, angle, dNN, dsum, bp-org,  x1,  y1,   z1,   x2,   y2,   z2
 * col#  1   2    3     4     5    6-8    9-11 12-14 15-17 18-20 21-23 24-26
 * bpid: 0: not-paired; +1: WC geometry; +2: WC pair; -1: other cases
 *       bp-pars,   bp_relative_orientation, C1_N_relative_orientation
 * col#   27-32        33-35 (dx, dy, dz)              36
 *
 * There are also possibilities when the pair is mostly maintained by
 * backbone + backbone; backbone + base H-bonds. e.g. A512 + C637
 * in rr0033 (chain 0). No direct base-to-base H-bonds **************
 */
{
    double OVLMT = 0.01, dd, dir_x, dir_y;
    double dorg[4], oave[4], zave[4], dNN_vec[4];
    long cdns, m, n;

    *bpid = 0;  /* default as not-paired */
    if (i == j)
        return;  /* same residue */

    get_bp_zoave(i, j, orien, org, oave, zave);

    ddxyz(org[i], org[j], dorg);
    ddxyz(NC1xyz[i], NC1xyz[j], dNN_vec);

    rtn_val[1] = veclen(dorg);  /* distance between origins */
    dir_x = dot(&orien[i][0], &orien[j][0]);  /* relative x direction */
    dir_y = dot(&orien[i][3], &orien[j][3]);  /* relative y direction */
    dd = dot(&orien[i][6], &orien[j][6]);  /* relative z direction */

    rtn_val[2] = fabs(dot(dorg, zave));  /* dv: projection onto mean normal */
    rtn_val[3] = 90.0 - fabs(dot2ang(dd) - 90.0);  /* angle between base normals */
    rtn_val[4] = veclen(dNN_vec);  /* RN9-YN1 distance */
    rtn_val[5] = rtn_val[1] + 2.0 * rtn_val[2];

    if (network) {  /* check if two bases in pairing network */
        if (  /* rtn_val[2] <= (BPRS[3] * DV_FACTOR) &&
                   [this condition can be marked out: better result] */
               rtn_val[3] <= BPRS[4] && rtn_val[4] >= BPRS[5] &&
               get_oarea(i, j, ring_atom, oave, zave, xyz, 0) < OVLMT)
            *bpid = 1;
        return;
    }

    cdns = rtn_val[1] <= BPRS[2] && rtn_val[2] <= BPRS[3]
        && rtn_val[3] <= BPRS[4] && rtn_val[4] >= BPRS[5];
    if (cdns && get_oarea(i, j, ring_atom, oave, zave, xyz, 0) < OVLMT)  /* no-overlap */
        for (m = seidx[i][1]; m <= seidx[i][2]; m++)
            for (n = seidx[j][1]; n <= seidx[j][2]; n++)
                if (is_baseatom(AtomName[m]) && is_baseatom(AtomName[n]) &&
                    (idx[m] == 2 || idx[m] == 4) && (idx[n] == 2 || idx[n] == 4) &&
                    /* O[2] or N[4] atoms on base: maybe exclude RN9/YN1 */
                    within_limits(xyz[m], xyz[n], 0, BPRS[1])) {  /* H-bond upper limit */
                    more_bppars(i, j, dir_x, dir_y, dd, orien, org, bseq, NC1xyz, rtn_val, bpid);  /* bpid alread a long *pointer* */
                    return;
                }  /* at least one pair of base O & N is within range */
}

void more_bppars(long i, long j, double dir_x, double dir_y, double dd, double **orien,
                 double **org, char *bseq, double **NC1xyz, double *rtn_val, long *bpid)
/* once a pair is found, calculate more parameters */
{
    char bpi[3];
    static char *WC[9] = { WC_LIST };
    long k, l, koffset;
    double zave[4], dNN_vec[4], pars[7], **r1, **r2, **mst;

    r1 = dmatrix(1, 3, 1, 3);
    r2 = dmatrix(1, 3, 1, 3);
    mst = dmatrix(1, 3, 1, 3);

    for (k = 1; k <= 9; k++) {
        rtn_val[k + 8] = orien[i][k];  /* base I xyz axes */
        rtn_val[k + 17] = orien[j][k];  /* base II xyz axes */
    }
    for (k = 1; k <= 3; k++) {
        koffset = (k - 1) * 3;
        for (l = 1; l <= 3; l++) {
            r1[l][k] = orien[i][koffset + l];
            r2[l][k] = (k == 1 || dd > 0) ?  /* keep x, reverse y & z if anti-parallel */
                orien[j][koffset + l] : -orien[j][koffset + l];
        }
    }
    bpstep_par(r2, org[j], r1, org[i], pars, mst, &rtn_val[5]);
    for (k = 1; k <= 6; k++)  /* bp parameters in columns 27-32 */
        rtn_val[26 + k] = pars[k];

    sprintf(bpi, "%c%c", toupper(bseq[i]), toupper(bseq[j]));
    *bpid = -1;  /* assumed to be non-WC */
    if (dir_x > 0.0 && dir_y < 0.0 && dd < 0.0) {
        *bpid = 1;  /* with WC geometry */
        if (fabs(pars[1]) <= WC_SHEAR && num_strmatch(bpi, WC, 1, 8)) {
            *bpid = 2;  /* WC */
            rtn_val[5] -= 1.5;  /* bonus for WC pair: Luger pd0171 */
        }
    }

    rtn_val[33] = dir_x;
    rtn_val[34] = dir_y;
    rtn_val[35] = dd;

    if (NC1xyz[i][7] > 0 && NC1xyz[j][7] > 0) {
        ddxyz(NC1xyz[i], NC1xyz[i] + 3, zave);
        ddxyz(NC1xyz[j], NC1xyz[j] + 3, dNN_vec);
        vec_norm(zave);
        vec_norm(dNN_vec);
        rtn_val[36] = dot(zave, dNN_vec);
    } else
        rtn_val[36] = EMPTY_NUMBER;

    free_dmatrix(r1, 1, 3, 1, 3);
    free_dmatrix(r2, 1, 3, 1, 3);
    free_dmatrix(mst, 1, 3, 1, 3);
}

void o3_p_xyz(long ib, long ie, char *aname, char **AtomName, double **xyz,
              double *o3_or_p, long idx)
/* get the xyz coordinates of O3' and P atoms */
{
    long i;

    i = find_1st_atom(aname, AtomName, ib, ie, "");
    if (i) {  /* with O3'/P atom */
        cpxyz(xyz[i], o3_or_p + idx - 4);
        o3_or_p[idx] = 1.0;
    } else  /* without O3'/P atom */
        o3_or_p[idx] = -1.0;
}

long is_baseatom(char *atomname)
/* Check if atomname is a NA base atom: ASSUMING NA residue */
{
    return (atomname[0] == ' ' && strchr("HP", atomname[1]) == NULL
            && isdigit((int) atomname[2]) && strchr(" M", atomname[3])) ? 1 : 0;
}

void base_info(long num_residue, char *bseq, long **seidx, long *RY, char **AtomName,
               char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
               double **orien, double **org, double **NC1xyz, double **o3_p, double *BPRS)
/* get base information for locating possible pairs later */
{
    char BDIR[BUF512], str[BUF512];
    long i, ib, ie, j;
    FILE *fpar;

    get_BDIR(BDIR, "Atomic_A.pdb", 1);

    /* get the reference frame for each base */
    base_frame(num_residue, bseq, seidx, RY, AtomName, ResName, ChainID,
               ResSeq, Miscs, xyz, BDIR, orien, org);

    /* read in base-pair criteria parameters */
    get_BDIR(BDIR, PAR_FILE, 1);
    strcat(BDIR, PAR_FILE);
    fpar = open_file(BDIR, "r");
    fprintf(stderr, " ...... reading file: %s ...... \n", PAR_FILE);
    for (i = 1; i <= 6; i++)
        if ((fgets(str, sizeof str, fpar) == NULL)
            || (sscanf(str, "%lf", &BPRS[i]) != 1))
            fatal("error reading base-pair criterion parameters\n");
    if (fabs(BPRS[4]) > 90.0)  /* make sure maximum 90 degrees */
        BPRS[4] = 90.0;
    BPRS[4] = 90.0 - fabs(fabs(BPRS[4]) - 90.0);
    if (BPRS[6] < 0.0)
        BPRS[6] = 7.5;  /* default to 7.5 A */
    close_file(fpar);

    for (i = 1; i <= num_residue; i++) {
        ib = seidx[i][1];
        ie = seidx[i][2];
        if (RY[i] >= 0) {  /* a base */
            j = (RY[i] == 1) ? find_1st_atom(" N9 ", AtomName, ib, ie, "") :
                (bseq[i] == 'P'
                 || bseq[i] == 'p') ? find_1st_atom(" C5 ", AtomName, ib, ie, "") :
                find_1st_atom(" N1 ", AtomName, ib, ie, "");
            cpxyz(xyz[j], NC1xyz[i]);  /* RN9/YN1 surely exist */
            j = find_1st_atom(" C1'", AtomName, ib, ie, "");
            if (j) {
                NC1xyz[i][7] = 1.0;  /* C1' atom exists */
                cpxyz(xyz[j], NC1xyz[i] + 3);  /* 4-6 */
            }
            o3_p_xyz(ib, ie, " O3'", AtomName, xyz, o3_p[i], 4);
            o3_p_xyz(ib, ie, " P  ", AtomName, xyz, o3_p[i], 8);
        }
    }
}

void bpstep_par(double **rot1, double *org1, double **rot2, double *org2, double *pars,
                double **mst_orien, double *mst_org)
/* calculate step or base-pair parameters (CEHS scheme) */
{
    double phi, rolltilt;
    double hinge[4], mstx[4], msty[4], mstz[4], t1[4], t2[4];
    double **para_bp1, **para_bp2, **temp;
    long i, j;

    for (i = 1; i <= 3; i++) {
        t1[i] = rot1[i][3];  /* z1 */
        t2[i] = rot2[i][3];  /* z2 */
    }

    cross(t1, t2, hinge);
    rolltilt = magang(t1, t2);

    /* Handle the special cases where t1 and t2 are perfectly
       anti-parallel, thus hinge can not be defined by t1 and t2. The
       vector would be (0 0 0), corresponding to tilt^2 + roll^2 = 180
       Date: July 20, 2005 */
    if (veclen(hinge) < XEPS && (fabs(rolltilt - 180) < XEPS || rolltilt < XEPS))
        for (i = 1; i <= 3; i++)
            hinge[i] = rot1[i][1] + rot2[i][1] + rot1[i][2] + rot2[i][2];

    para_bp1 = dmatrix(1, 3, 1, 3);
    para_bp2 = dmatrix(1, 3, 1, 3);
    temp = dmatrix(1, 3, 1, 3);

    arb_rotation(hinge, -0.5 * rolltilt, temp);
    multi_matrix(temp, 3, 3, rot2, 3, 3, para_bp2);
    arb_rotation(hinge, 0.5 * rolltilt, temp);
    multi_matrix(temp, 3, 3, rot1, 3, 3, para_bp1);

    for (i = 1; i <= 3; i++) {
        mstz[i] = para_bp2[i][3];  /* also para_bp1(:,3) */
        t1[i] = para_bp1[i][2];  /* y1 */
        t2[i] = para_bp2[i][2];  /* y2 */
    }

    /* twist is the angle between the two y- or x-axes */
    pars[6] = vec_ang(t1, t2, mstz);

    /* get y- and x-axes */
    get_vector(t1, mstz, 0.5 * pars[6], msty);
    cross(msty, mstz, mstx);

    avexyz(org1, org2, mst_org);
    ddxyz(org1, org2, t1);
    x_y_z_2_mtx(mstx, msty, mstz, mst_orien);

    /* get the xyz displacement parameters */
    for (i = 1; i <= 3; i++) {
        pars[i] = 0.0;
        for (j = 1; j <= 3; j++)
            pars[i] += t1[j] * mst_orien[j][i];
    }

    /* phi angle is defined by hinge and msty */
    phi = deg2rad(vec_ang(hinge, msty, mstz));

    /* get roll and tilt angles */
    pars[5] = rolltilt * cos(phi);
    pars[4] = rolltilt * sin(phi);

    free_dmatrix(para_bp1, 1, 3, 1, 3);
    free_dmatrix(para_bp2, 1, 3, 1, 3);
    free_dmatrix(temp, 1, 3, 1, 3);
}

void help3dna(char *program_name)
/* Easy 3DNA help message read from a text file */
{
    char BDIR[BUF512], str[BUF512], pname[BUF512];
    long nlen, found_help = 0;
    FILE *fp;

    get_BDIR(BDIR, HELP3DNA, 0);
    strcat(BDIR, HELP3DNA);
    fp = open_file(BDIR, "r");

    strcpy(pname, program_name);  /* program_name is a constant */
    nlen = upperstr(pname);  /* to upper case */
    while (fgets(str, sizeof str, fp) != NULL) {
        if (str[0] != '<')
            continue;
        upperstr(str);
        if (strncmp(str + 1, pname, nlen))
            continue;  /* not match */
        found_help = 1;
        while (fgets(str, sizeof str, fp) != NULL) {
            if (str[0] == '<') {
                upperstr(str);
                if (str[1] != '/' || (strncmp(str + 2, pname, nlen)))
                    fatal("error in help format: <tag> ... </tag>\n");
                found_help = 2;
                goto FINISHED;
            } else {
                if (str[0] != '#')  /* not-comment */
                    fprintf(stderr, "%s", str);  /* str includes \n */
            }
        }
    }

  FINISHED:
    if (!found_help)
        fprintf(stderr, "No help found for program <%s>\n", program_name);
    else if (found_help != 2)
        fprintf(stderr, "Warning: no matching tag found\n");

    close_file(fp);
}

double get_rowpar(long row_num)
/* read in single parameter from "misc_3dna.par" based on row #:
 * #7 for deciding if a helix is strongly curved to output helical axis
 * #8 for H-bond distance for hydration
 * #9 for maximum O3'--P distance for linkage (rebuild) */
{
    char BDIR[BUF512], str[BUF512];
    double rowpar;
    long i, j, k;
    FILE *fp;

    get_BDIR(BDIR, PAR_FILE, 0);
    strcat(BDIR, PAR_FILE);

    fp = open_file(BDIR, "r");
    for (i = 1; i <= row_num; i++)
        if (fgets(str, sizeof str, fp) == NULL)
            fatal("error reading parameter file\n");
    if (row_num == 8) {  /* for hydration */
        k = sscanf(str, "%lf %ld", &rowpar, &j);
        if (k != 1 && k != 2)
            fatal("error reading hydration distance criterion\n");
        if (k == 1 || (k == 2 && j != 1 && j != 9))
            j = 0;  /* 0 for N/O; 1 for N/O/C; 9 for all atoms */
        rowpar = MFACTOR * rowpar + j;
    } else {
        if (sscanf(str, "%lf", &rowpar) != 1)
            fatal("error reading file <misc_3dna.par>\n");
    }
    close_file(fp);

    return rowpar;
}

void delH_pdbfile(char *inpfile, char *outfile)
/* Delete all H-atoms and write the rest coordinates in PDB format */
{
    char *ChainID, *p, **AtomName, **ResName, **Miscs;
    double **xyz;
    long i, k = 0, num, *ResSeq, *idx;
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

    idx = lvector(1, num);
    atom_idx(num, AtomName, Miscs, idx);  /* atom index to mark H */

    fp = open_file(outfile, "w");
    print_pdb_title(inpfile, "*", fp);
    for (i = 1; i <= num; i++) {
        if (idx[i] == 3)  /* H-atom */
            continue;
        p = Miscs[i];
        fprintf(fp, "%s%5ld %4s%c%3s %c%4ld%c   %8.3f%8.3f%8.3f%s\n",
                (p[0] == 'A') ? "ATOM  " : "HETATM", ++k, AtomName[i],
                p[1], ResName[i], ChainID[i], ResSeq[i], p[2], xyz[i][1],
                xyz[i][2], xyz[i][3], p + 3);
    }
    fprintf(fp, "END\n");
    close_file(fp);

    free_pdb(num, NULL, AtomName, ResName, ChainID, ResSeq, xyz, Miscs);
    free_lvector(idx, 1, num);
}

/* to avoid waring with ignoring return value of fgets() */
static void my_fgets(char *s, int n, FILE * fp)
{
    char *p;

    p = fgets(s, n, fp);
    (void)(p);
}

void contact_msg(long prt_msg)
{
    if (prt_msg) {
        char str[BUF512];
        fprintf(stderr, "%s\n", X3DNA_VER);
        fprintf(stderr, "Press <enter> to continue ... ");
        my_fgets(str, sizeof str, stdin);  /* just pause to view previous output */
        help3dna("contact_info");
    }
    exit(1);
}
