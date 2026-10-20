
/* alc2img.c */
void alc2img_usage(void);
void alc2img_cmdline(int argc, char *argv[], char *inpfile, char *outfile,
                     double *scale_factor, long *opts);

/* ana_fncs.c */
long **read_input(char *inpfile, char *pdbfile, char *outfile, long *ds, long *num_bp,
                  long *ip, long *hetatm);
void print_header(long num_bp, long num, char *pdbfile, FILE * fp);
void atom_list(long ds, long num_bp, long **pair_num, long **seidx, long *RY,
               char **bp_seq, char **AtomName, char **ResName, char *ChainID,
               long *ResSeq, char **Miscs, long **phos, long **c6_c8, long **sugar,
               long **chi);
void get_nt_bb_torsion(double **nt_bb_torsion, long num_residue, long **seidx,
                       long *RY, char **AtomName, char **ResName, char *ChainID,
                       long *ResSeq, char **Miscs, double **xyz);
double idx2torsion(long *idx, double **xyz);
void backbone_torsion(long ds, long num_bp, long **pair_num, char **bp_seq, long **sugar,
                      long **chi, double **xyz, double **nt_bb_torsion, FILE * fp);
void p_c1_dist(long ds, long num_bp, char **bp_seq, long **phos, long **chi,
               double **xyz, long *bphlx, FILE * fp);
void lambda_d3(long num_bp, char **bp_seq, long **chi, long **c6_c8, double **xyz,
               FILE * fp);
void print_axyz(long num_bp, char **bp_seq, long **aidx, char *aname, double **xyz);
double gw_angle(long *idx, double *pvec12, double **xyz);
void groove_width(long parallel, long num_bp, char **bp_seq, long **phos,
                  double **xyz, long *bphlx, FILE * fp);
void ref_frames(long ds, long num_bp, long **pair_num, char **bp_seq, long **seidx,
                long *RY, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                char **Miscs, double **xyz, FILE * fp, double **orien, double **org,
                long *WC_info, long *str_type, long irna, long **o3p_brk);
void check_Watson_Crick(long num_bp, char **bp_seq, double **orien, double **org,
                        long *WC_info);
void helical_par(double **rot1, double *org1, double **rot2, double *org2, double *pars,
                 double **mst_orien, double *mst_org);
void print_par(char **bp_seq, long num_bp, long ich, long ishel, double **param,
               FILE * fp);
void single_helix(long num_bp, char **bp_seq, double **step_par, double **heli_par,
                  double **orien, double **org, FILE * fp);
void prt_stepstr(char **step_str, long num_step, long *bphlx, long ishel, double **param,
                 FILE * fp);
void prt_step_par(char **bp_seq, long num_bp, long *bphlx, long ishel, double **param,
                  FILE * fp);
void double_helix(long num_bp, char **bp_seq, double **step_par, double **heli_par,
                  double **orien, double **org, FILE * fp, double **twist_rise,
                  double *mst_orien, double *mst_org, double *mst_orienH,
                  double *mst_orgH, long *bphlx, long istart, long istep,
                  long bz, long *str_type);
void bz_check(double **r1, double *o1, double **r2, double *o2, long bz,
              long *bz_junction, long *z_step);
double get_mtwist(long nbpm1, long *bphlx, long *WC_info, double **twist_rise);
void get_parameters(long ds, long num_bp, char **bp_seq, double **orien, double **org,
                    FILE * fp, double **twist_rise, double *mst_orien, double *mst_org,
                    double *mst_orienH, double *mst_orgH, long *bphlx, long istart,
                    long istep, long bz, long *str_type);
void parvec2mtx(double *parvec, long num, double **parmtx);
void print_ss_rebuild_pars(double **pars, long num_bp, char *str, char **bp_seq,
                           FILE * fp);
void print_ds_rebuild_pars(double **bp_par, double **step_par, long num_bp, char *str,
                           char **bp_seq, FILE * fp);
void print_ref(char **bp_seq, long num_item, long ich, double *org, double *orien,
               FILE * fp);
void write_mst(long num_bp, long **pair_num, char **bp_seq, double *mst_orien,
               double *mst_org, long **seidx, char **AtomName, char **ResName,
               char *ChainID, long *ResSeq, double **xyz, char **Miscs,
               long **htm_water, char *strfile);
void print_xyzP(long parallel, long nbpm1, char **bp_seq, long **phos, double *mst_orien,
                double *mst_org, double **xyz, FILE * fp, char *title_str, double **aveP,
                long p_offset);
void print_PP(long parallel, double mtwist, double **twist_rise, long num_bp,
              char **bp_seq, long **phos, double *mst_orien, double *mst_org,
              double *mst_orienH, double *mst_orgH, double **xyz, long *WC_info,
              long *bphlx, FILE * fp);
void str_classify(double mtwist, long str_type, long parallel, long num_bp, FILE * fp);
double a_hlxdist(long idx, double **xyz, double *hlx_axis, double *hlx_pos);
void print_radius(char **bp_seq, long nbpm1, long ich, double **p_radius,
                  double **o4_radius, double **c1_radius, long *bphlx, FILE * fp);
void helix_radius(long ds, long num_bp, char **bp_seq, double **orien, double **org,
                  long **phos, long **chi, double **xyz, long *bphlx, FILE * fp);
void print_shlx(char **bp_seq, long nbpm1, long ich, double *shlx_orien,
                double *shlx_org, FILE * fp);
void get_helix_axis(long ds, long num_bp, char **bp_seq, double **orien,
                    double **org, long *bphlx, FILE * fp);
void print_analyze_ref_frames(long ds, long num_bp, char **bp_seq, double *iorg,
                              double *iorien);
void get_axis(long nvec, long **idx, long num, double **xyz, long nb, long *C1b,
              long *C1e, double *std_rise, double *hrise, double *haxis,
              double *hstart, double *hend);
void print_poc_r3d(double *rave, double *hstart, double *hend);
void get_poc_radius(long ds, long num_bp, long **phos, long **chi, double **xyz,
                    double *haxis, double *hstart, double *hend, double *rave,
                    double *rstd);
void global_c1_c1_par(long num_bp, char **bp_seq, long **chi, double **xyz,
                      double *haxis, double *hstart, FILE * fp);
void global_analysis(long ds, long num_bp, long num, char **bp_seq, long **chi,
                     long **phos, double **xyz, FILE * fp);
void base_overlap(long ds, long num_bp, long num, long num_residue, long **pair_num,
                  long *bRY, char **bp_seq, long **seidx, char **AtomName, double **xyz,
                  double **orien, double **org, FILE * fp);
double get_oarea(long r1, long r2, long **ring_atom, double *oave, double *zave,
                 double **xyz, long only_ring);
long ratom_xyz(long *ratom_list, long only_ring, double **xyz, double *oave,
               double **oxyz);
void get_zoave(long istep, long ds, double **orien, double **org, double *oave,
               double *zave);
void get_bp_zoave(long ia, long ib, double **orien, double **org, double *oave,
                  double *zave);
void ring_oidx(long num, long num_residue, long *RY, long **seidx, char **AtomName,
               double **xyz, long **ring_atom);
void get_cntatom(long *ringlist, long **connect, long *idx);
void pia_fit(double minx, double miny, const double mid, double sclx, double scly,
             point * x, long cx, vertex * ix, long fudge);
double pia_area(point a, point p, point q);
void pia_cntrib(double *s, point f, point t, long w);
long pia_ovl(point p, point q);
void pia_cross(double *out_s, vertex * a, vertex * b, vertex * c, vertex * d,
               double a1, double a2, double a3, double a4);
void pia_inness(double *out_s, vertex * P, long cP, vertex * Q, long cQ);
double pia_inter(point * a, long na, point * b, long nb);
void cehs_base_atoms(char **AtomName, long ib, long ie, long *num_batom, long *batom);
void cehs_bppar(double **rot1, double *org1, double **rot2, double *org2, double *pars,
                double **mst_orien, double *mst_org);
void cehs_pars(long num_bp, long istart, long istep, long **pair_num, char **bp_seq,
               long **seidx, long **c6_c8, long *RY, char **AtomName, char **ResName,
               char *ChainID, long *ResSeq, char **Miscs, double **xyz, double *bp_orien,
               double *bp_org, long bz, FILE * fp);
void schnaap_global(long num_bp, long num, char **bp_seq, long **chi, double **xyz,
                    double *bp_orien, double *bp_org, FILE * fp);
void out_cehs(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org,
              FILE * fp);
void compdna(double **rot1, double *org1, double **rot2, double *org2, double *pars,
             double **mst_orien, double *mst_org);
void out_compdna(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org,
                 FILE * fp);
void my_curves(double **rot1, double *org1, double **rot2, double *org2, double *pars);
void curves_mbt(long ibp, double **orien, double **org, double **cvr, double *cvo);
void out_curves(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org,
                FILE * fp);
void freehelix(double **rot1, double *org1, double **rot2, double *org2, double *pars,
               double **mst_orien, double *mst_org);
void out_freehelix(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org,
                   FILE * fp);
void sgl_helix(double **rot1, double **rot2, double *rot_ang, double *rot_hlx);
void ngeom(double **rot1, double *org1, double **rot2, double *org2, double *pars,
           double **mst_orien, double *mst_org);
void out_ngeom(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org,
               FILE * fp);
void nuparm(double **rot1, double *org1, double **rot2, double *org2, double *pars,
            double **mst_orien, double *mst_org, double *hpars, long get_hpar);
void out_nuparm(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org,
                FILE * fp);
void rna(double **rot1, double *org1, double *pvt1, double **rot2, double *org2,
         double *pvt2, double *pars, double **mst_orien, double *mst_org);
void pvt_dxdy(double **rot1, double *org1, double *pvt1, double *pars,
              double **mst_orien, double *mst_org);
void out_rna(long ds, long num_bp, char **bp_seq, long *bphlx, double **orien,
             double **org, FILE * fp);
void other_pars(long num_bp, char **bp_seq, long *bphlx, double **orien, double **org);

/* analyze.c */
void analyze_usage(void);
void process_str(char *inpfile, long icnt, long waters, long istart, long istep, long bz);

/* anyhelix.c */
void anyhelix_usage(void);
long **read_pair(char *inpfile, char *pdbfile, char *outfile, long *num_bp,
                 long *hetatm, long *eqnum);
void cross_pair(long num, long num_residue, long num_bp, long **pair_num,
                char *HB_ATOM, double *HB_UPPER, double *BPRS, char *bseq,
                long **seidx, char **AtomName, char **ResName, char *ChainID,
                long *ResSeq, char **Miscs, double **xyz, double *chi,
                double **orien, double **org, double **NC1xyz, double *bp_orien,
                double *bp_org, long *RY, long *idx, long **htm_water, FILE * fp);
void stepstr(long i, long j, char *bseq, long **pair_num, char *str);
void print_ah_par(long ishel, long num_bp, char *bseq, long **pair_num,
                  double **param, FILE * fp);
void ah_write_mst(long i, long **pair_num, char *bseq, long **seidx, double **mst,
                  double *msto, char **AtomName, char **ResName, char *ChainID,
                  long *ResSeq, char **Miscs, double **xyz, long **htm_water, FILE * fp);
void ah_step_hel(long num_bp, long **pair_num, char *bseq, long **seidx,
                 double *bp_orien, double *bp_org, char **AtomName, char **ResName,
                 char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                 long **htm_water, FILE * fp);
void get_ah_helix(long num_bp, char *bseq, long **pair_num, long num, long **idxCN,
                  long **seidx, char **AtomName, double **xyz, double **NC1xyz,
                  double *bp_org, FILE * fp);
void get_ah_poc_radius(long num_bp, long **pair_num, long **seidx, char **AtomName,
                       double **xyz, double *haxis, double *hstart, double *hend,
                       double *rave, double *rstd);
void print_ah_frame(long num_bp, long **pair_num, char *bseq, double *bp_orien,
                    double *bp_org);
void process_anyhelix(char *inpfile, long waters);
void stnd1_up(double **orien, double **org, double *chi, long num_bp,
              long **pair_num, FILE * fp);
void reverse_xz(long idx, double **orien, long **pair_num, long *s1dir);
void local_c1par(long num_bp, char *bseq, double *bp_orien, double **NC1xyz,
                 long **pair_num, FILE * fp);
void c1_twist_rise(long i, long nb, double *mz, double **NC1xyz, long **pair_num,
                   FILE * fp);
void layer_lsfit(long num_bp, char *bseq, long **pair_num, long **seidx,
                 char **AtomName, double **xyz, FILE * fp);

/* cehs.c */
void cehs_usage(void);
void process_cehs(char *inpfile, long irna, long istart, long istep, long bz);

/* cmn_fncs.c */
void print_used_time(time_t time0);
void parcat(char *str, double par, char *format, char *bstr);
void print_bp_crit(double *BPRS, char *HB_ATOM, FILE * fp);
char *my_getline(FILE * fp);
long csplit(char *str, char *item[], long itemsize, char sepc);
char *trim(char *a);
char *ltrim(char *a);
char *rtrim(char *a);
long itemize(char *str, char *item[], long itemsize);
void refs_right_left(long bnum, double **orien, double **org, double **r1, double *o1,
                     double **r2, double *o2);
void refs_i_j(long b1, long b2, double *bp_orien, double *bp_org, double **r1,
              double *o1, double **r2, double *o2);
void ref_frame_i(long bnum, double *bp_orien, double *bp_org, double **r, double *o);
void mst2orien(double *orien_vec, long ioffset, double **mst);
void orien2mst(double *orien_vec, long ioffset, double **mst);
void x_y_z_2_mtx(double *x, double *y, double *z, double **mtx);
void cehs_average(long inum_base, long *ivec, double **orien, double **org, double **mst,
                  double *morg);
void geom_average(long inum_base, long *ivec, double **orien, double **org, double **mst,
                  double *morg);
void pair2mst(long inum_base, long *ivec, char **AtomName, char **ResName, char *ChainID,
              long *ResSeq, char **Miscs, double **xyz, double **orien, double **org,
              long **seidx, double *mst_orien, double *mst_org, long **htm_water,
              FILE * fp);
void get_chi_angle(long num_residue, long *RY, char *bseq, long **seidx, double **xyz,
                   char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                   char **Miscs, double *chi, long **idxCN);
FILE *open_file(char *filename, char *filemode);
long close_file(FILE * fp);
long exist_file(char *filename);
void remove_file(char *filename);
void rename_file(char *src, char *dst);
void copy_file(char *src, char *dst);
void cat_file(char *src, char *dst);
long upperstr(char *a);
long lowerstr(char *a);
char *my_strdup(char *src);
void print_sep(FILE * fp, char x, long n);
void check_slash(char *BDIR);
void delete_end_slash(char *str);
char *basename(char *str);
long lround(double d);
void del_extension(char *fullname, char *okname);
void fatal(char *fmt, ...);
void print_pdb_title(char *pdbfile, char *chain_list, FILE * fp);
long number_of_atoms(char *pdbfile, long hetatm, char *ALT_LIST);
long read_pdb(char *pdbfile, long *AtomSNum, char **AtomName, char **ResName,
              char *ChainID, long *ResSeq, double **xyz, char **Miscs, long hetatm,
              char *ALT_LIST);
void free_pdb(long num, long *AtomSNum, char **AtomName, char **ResName, char *ChainID,
              long *ResSeq, double **xyz, char **Miscs);
void reset_xyz(long num, double **xyz, char *fmt);
void pdb_record(long ib, long ie, long *inum, long idx, char **AtomName, char **ResName,
                char *ChainID, long *ResSeq, double **xyz, char **Miscs, FILE * fp);
void write_pdb(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               double **xyz, char **Miscs, char *pdbfile);
void write_pdbcnt(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                  double **xyz, long **connect, char *pdbfile);
void move_position(double **d, long nr, long nc, double *mpos);
long **residue_idx(long num, long *ResSeq, char **Miscs, char *ChainID, char **ResName,
                   long *num_residue);
long residue_ident(char **AtomName, double **xyz, char **Miscs, long ib, long ie);
char base_ident(char *rname, char *idmsg, char **baselist, long num_sb);
void get_baselist(char **baselist, long *num_sb);
void get_seq(long num_residue, long **seidx, char **AtomName, char **ResName,
             char *ChainID, long *ResSeq, char **Miscs, double **xyz, char *bseq,
             long *RY);
void get_bpseq(long ds, long num_bp, long **pair_num, long **seidx, char **AtomName,
               char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
               char **bp_seq, long *RY);
long num_strmatch(char *str, char **strmat, long nb, long ne);
long find_1st_atom(char *str, char **strmat, long nb, long ne, char *idmsg);
double torsion(double **d);
double torsion2(double **d);
void get_BDIR(char *BDIR, char *filename, long prt_msg);
void align2zaxis(long num, double *haxis, double **rotmat, double **xyz, double **xyzH);
void cov_matrix(double **a, double **b, long nr, long nc, double **cmtx);
double ls_fitting(double **sxyz, double **exyz, long n, double **fitted_xyz, double **R,
                  double *orgi);
void ls_plane(double **bxyz, long n, double *pnormal, double *ppos, double *odist,
              double *adist);
void arb_rotation(double *va, double ang_deg, double **rot_mtx);
double vec_ang(double *va, double *vb, double *vref);
void get_vector(double *va, double *vref, double deg_ang, double *vo);
void rotate(double **a, long i, long j, long k, long l, double *g, double *h, double s,
            double tau);
void eigsrt(double *d, double **v, long n);
void jacobi(double **a, long n, double *d, double **v);
void dludcmp(double **a, long n, long *indx, double *d);
void dlubksb(double **a, long n, long *indx, double *b);
void dinverse(double **a, long n, double **y);
void rotx(double ang_deg, double **rotmat);
void roty(double ang_deg, double **rotmat);
void rotz(double ang_deg, double **rotmat);
void get_alc_nums(char *alcname, long *num, long *nbond);
void read_alc(char *alcname, long *num, long *nbond, char **AtomName, double **xyz,
              long *ibase, long **linkage);
void write_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
               long **linkage, char *alcfile);
void free_alc(long num, long nbond, char **AtomName, double **xyz, long *ibase,
              long zero_1, long **linkage);
void cnct_org(long num_bp, long ia, long ib, char **tAtomName, double **txyz,
              long *tibase, long **tlinkage, double **org_xyz);
long get_num_bp(char *inpfile);
void dsort(long n, double *a, long *idx);
void lsort(long n, long *a, long *idx);
void lreverse(long ia, long n, long *lvec);
void fig_title(FILE * fp);
void ps_title_cmds(FILE * fp, char *imgfile, long *bbox);
void get_fig_xy(long num, double **xyz, long nO, double **oxyz, long *urxy,
                long frame_box, FILE * fp);
void get_pxy(double *xy1, double *xy2, double r, double *px, double *py);
void alc2fig(long nobj, long *idx, long *depth, long **allobj, double **blkxyz,
             double **oxyz, long *ibase, long faces[][5], long *opts, FILE * fp);
void get_ps_xy(char *imgfile, long *urxy, long frame_box, FILE * fp);
void alc2ps(long nobj, long *idx, long **allobj, double **blkxyz, double **oxyz,
            long *ibase, long faces[][5], long *opts, FILE * fp);
void bring_atoms(long ib, long ie, long ra_num, char **AtomName, long *nmatch,
                 long *batom);
void all_bring_atoms(long num_residue, long *RY, long **seidx, char **AtomName,
                     long *num_ring, long **ring_atom);
void base_idx(long num, char *bseq, long *ibase, long single);
long basepair_idx(char *bpi);
void plane_xyz(long num, double **xyz, double *ppos, double *nml, double **nxyz);
void prj2plane(long num, long ra_num, char **AtomName, double **xyz, double z0,
               double **nxyz);
void adjust_xy(long num, double **xyz, long nO, double **oxyz, double scale_factor,
               long default_size, long *urxy);
void get_depth(long nobj, long *zval, long *depth);
void raster3d_header(long num, double **xyz, double scale_factor, long no_header,
                     long frame_box, FILE * fp);
void get_r3dpars(double **base_col, double *hb_col, double *width3, double **atom_col,
                 char *label_style);
void r3d_rod(long itype, double *xyz1, double *xyz2, double rad, double *rgbv, FILE * fp);
void r3d_dash(double *xyz1, double *xyz2, double hb_width, double *hb_col, FILE * fp);
void r3d_sphere(double *xyz1, double rad, double *rgbv, FILE * fp);
void cpk_model(long num, long *idx, double **xyz, double ballrad, double **colrgb,
               FILE * fp);
void r3d_tripln(long itype, double *xyz1, double *xyz2, double *xyz3, double *rgbv,
                FILE * fp);
void r3d_block_edge(double *rgbv, long ioffset8, double **blkxyz, double w1, FILE * fp);
void base_label(double **rxyz, char *label_style, double *rgbv, char *bname_num,
                FILE * fp);
void fill_base_ring(long num_residue, long num_ring, long **ring_atom, double **xyz,
                    long *ibase, char *bseq, double **base_col, char *label_style,
                    long label_ring, long *ResSeq, FILE * fp);
void process_alc(char *alcfile, char *imgfile, double scale_factor, long *opts);
void alc_3images(long *opts, long nobj, long num_blk, long num_blk8, long nO_lkg,
                 long nO, double **oxyz, double **blkxyz, long *blkibase, long **linkage,
                 double scale_factor, char *imgfile);
void get_alc_objs(long num_blk, double **blkxyz, long nO, double **oxyz, long nO_lkg,
                  long **linkage, long faces[][5], long **allobj);
void get_fig_pars(double *dot_sep, long *dlcol, long *dwidth, long *bp1width,
                  long *bp2width, long bc_idx[][7], double *msat, double *Msat,
                  long *o_sides, long *line_width, long *join_style, long *cap_style,
                  long *mfcol);
void change_xyz(long side_view, double *morg, double **mst, long num, double **xyz);
void get_side_view(long ib, long ie, double **xyz);
void get_CNidx(long ds, long num_bp, long **chi, long **idx, long *nvec, long *C1b,
               long *C1e);
void add_3axes(long *num, char **AtomName, long *ibase, double **xyz, long *nbond,
               long **linkage, long side_view);
void bp_sequence(char *Wbase, char *Cbase, long *num_bp, char **bp_seq);
void combine_pstnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                    char **tAtomName, char **tResName, char *tChainID, long *tResSeq,
                    double **txyz, char **tAtomName2, double **txyz2);
void reverse_stnd2(long num_bp, char **bp_seq, long **s2idx, long *tnum,
                   char **tAtomName, char **tResName, char *tChainID, long *tResSeq,
                   double **txyz, char **tAtomName2, double **txyz2, long basep);
void pair_checking(long ip, long ds, long num_residue, char *pdbfile, long *num_bp,
                   long **pair_num);
void drct_checking(long ds, long num_bp, long **pair_num, long **seidx, char **AtomName,
                   double **xyz, long *parallel, long *bbexist, long **o3p_brk, FILE * fp);
void base_frame(long num_residue, char *bseq, long **seidx, long *RY, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                char *BDIR, double **orien, double **org);
void base_str(char chain_id, long res_seq, char *misc, char *rname, char bcode,
              long stnd, char *idmsg);
void write_lkglist(long nbond, long **linkage, char **AtomName, char **ResName,
                   char *ChainID, long *ResSeq, char **Miscs);
void hbond_info(long **pair_num, char *bseq, long **seidx, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                long *RY, long *num_hbond, long **hb_linkage);
void hbond_pdb(long num, long num_residue, char *bseq, long **seidx, long *RY,
               char **AtomName, char **ResName, char *ChainID, long *ResSeq,
               char **Miscs, double **xyz, long *num_hbond, long **hb_linkage, long pwise);
void hbond_list(long i, long j, char **AtomName, char **ResName, char *ChainID,
                long *ResSeq, double **xyz, char **Miscs, char *bseq, long **seidx,
                long **hb_linkage, double *HB_UPPER, char *HB_ATOM, long **num_list,
                long *num_hbond, long ilayer, FILE * fph);
void hb_numlist(long i, long j, char basei, char basej, long **seidx, char **AtomName,
                double **xyz, double *HB_UPPER, char *HB_ATOM, long *num_hb,
                long **num_list);
void hb_information(long num_bp, long **pair_num, char **bp_seq, long **seidx,
                    char **AtomName, double **xyz, long *WC_info, FILE * fp);
long good_hbatoms(char *HB_ATOM, char *atom1, char *atom2);
void read_lkginfo(char *lkgfile, long num, long *nbond, long **linkage);
void read_hbinfo(char *hbfile, long num, long *num_hbond, long **hb_linkage);
void hb_atompair(long num_hbonds, char **hb_atom1, char **hb_atom2, double *hb_dist,
                 double *HB_UPPER);
void get_hbond_ij(long i, long j, char basei, char basej, double *HB_UPPER, long **seidx,
                  char **AtomName, char *HB_ATOM, double **xyz, char *HB_INFO);
char donor_acceptor(char basei, char basej, char *hb_atom1, char *hb_atom2);
void hb_crt_alt(double *HB_UPPER, char *HB_ATOM, char *ALT_LIST);
void atom_info(long idx, char atoms_list[NELE][3], double *covalence_radii,
               double *vdw_radii);
void aname2asym(char *aname0, char *my_asym);
void atom_idx(long num, char **AtomName, char **Miscs, long *idx);
void get_bonds(long num, char **AtomName, double **xyz, long num_residue, long *RY,
               long **seidx, long **connect);
void atom_linkage(long ib, long ie, long *idx, double **xyz, char **Miscs, char *ChainID,
                  long nbond_estimated, long *nbond, long **linkage);
void lkg2connect(long ib, long ie, long nbond, long **linkage, long **connect);
void init_htm_water(long waters, long num, long num_residue, long *idx, long **htm_water);
void identify_htw(long num_residue, long **seidx, long *RY, char **AtomName,
                  char **ResName, char *ChainID, long *ResSeq, char **Miscs,
                  double **xyz, long **htm_water);
long attached_residues(long inum_base, long *ivec, long *ivec2, long **seidx,
                       double **xyz, long **htm_water);
void print_pairinfo(long i, long j, char basei, char basej, double *rtn_val, double *chi,
                    char *HB_ATOM, double *HB_UPPER, long **seidx, char **AtomName,
                    double **xyz, char *bseq, long detailed, FILE * fp);
void check_pair(long i, long j, char *bseq, long **seidx, double **xyz, double **NC1xyz,
                double **orien, double **org, long *idx, char **AtomName, double *BPRS,
                double *rtn_val, long *bpid, long **ring_atom, long network);
void more_bppars(long i, long j, double dir_x, double dir_y, double dd, double **orien,
                 double **org, char *bseq, double **NC1xyz, double *rtn_val, long *bpid);
void o3_p_xyz(long ib, long ie, char *aname, char **AtomName, double **xyz,
              double *o3_or_p, long idx);
long is_baseatom(char *atomname);
void base_info(long num_residue, char *bseq, long **seidx, long *RY, char **AtomName,
               char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
               double **orien, double **org, double **NC1xyz, double **o3_p, double *BPRS);
void bpstep_par(double **rot1, double *org1, double **rot2, double *org2, double *pars,
                double **mst_orien, double *mst_org);
void help3dna(char *program_name);
double get_rowpar(long row_num);
void delH_pdbfile(char *inpfile, char *outfile);
void contact_msg(long prt_msg);

/* comb_str.c */
void comb_alcfile(char *inpfile1, char *inpfile2, char *outfile);
void comb_pdbfile(long num1, char *inpfile1, char *inpfile2, char *outfile);
long get_pdbconnect(long num, long **connect, char *pdbfile);

/* ex_str.c */
void exstr_usage(void);
void exstr_cmdline(int argc, char *argv[], char *inpfile, char *outfile, long *inum,
                   long *NMR, long *BIOU, long *DELH);
void check_ifNMR(char *pdbfile);
long best_conformer(char *pdbfile);
void extract_NMR(char *inpfile, char *outfile, long inum);
void extract_3DNA(char *inpfile, char *outfile, long inum);
void extract_BIOU(char *inpfile, char *outfile, long inum);
long biounit_cid(char *pdbfile, char **chlist);

/* fiber.c */
void fiber_dir(char *FDIR);
void fiber_usage(void);
void fiber_cmdline(int argc, char *argv[], char *FDIR, long *str_id, char *pdbfile);
void list_str_info(char *FDIR);
void get_twist_rise(long str_id, double *twist, double *rise);
void residue(double twist, double rise, char *b_rep, long s2r, char *FDIR, char *pdbfile);
void multi_bp(long str_id, double twist, double rise, char *FDIR, char *pdbfile);
void triplex(long str_id, double twist, double rise, char *FDIR, char *pdbfile);
void quadruplex(double twist, double rise, char *FDIR, char *pdbfile);
void shelix(double twist, double rise, char *FDIR, char *pdbfile);
void bp_repeat(char *Wbase, char *Cbase, char *crepeat, long *n_rep, char **bp_seq);
long repeat_num(void);

/* find_pair.c */
void fp_usage(void);
void fp_cmdline(int argc, char *argv[], char *inpfile, char *outfile, long *ds,
                long *curves, long *divide, long *hetatm, long *pairs,
                long *detailed, long *waters);
void handle_str(char *pdbfile, char *outfile, long ds, long curves, long divide,
                long hetatm, long pairs, long detailed, long waters);
void multi_bps(char *pdbfile, char *parfile);
void multiplets(long max_ple, long num_residue, long **pair_info, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, char **Miscs, double **xyz,
                double **orien, double **org, long **seidx, char *bseq, long hetatm,
                long **htm_water, FILE * fp);
void base_compatibility(long num_residue, long **pair_info, long *num_match,
                        long **match_list, long *num_partial, long **partial_list,
                        FILE * fp);
void print_list(long num_residue, long **mylist, char *my_str, FILE * fp);
long isequal_list(long *listA, long *listB);
void bp_network(long num_residue, long *RY, long **seidx, char **AtomName,
                char **ResName, char *ChainID, long *ResSeq, char **Miscs, long *idx,
                double **xyz, long **ring_atom, char *bseq, long **pair_info,
                double **NC1xyz, double **orien, double **org, double *BPRS,
                long hetatm, long **htm_water, FILE * fp);
long allbase_cncts(long i, long tnum_base, long *ivec, long **pair_info, FILE * fp);
void bases_elimination(long i, long inum_base, long *ivec, char *bseq, long **seidx,
                       long **ring_atom, double **xyz, double **NC1xyz, double **orien,
                       double **org, char **AtomName, double *BPRS, char *b1, long *idx,
                       long *max_ple, long **pair_info, FILE * fp);
void all_pairs(long num_residue, long *RY, double **NC1xyz, double **orien,
               double **org, double *BPRS, long **seidx, double **xyz, long *idx,
               long **ring_atom, char **AtomName, char **ResName, char *ChainID,
               long *ResSeq, char **Miscs, char *bseq, long hetatm, long **htm_water,
               char *pdbfile, char *outfile);
void best_pair(long i, long num_residue, long *RY, long **seidx, double **xyz, long *idx,
               double **NC1xyz, long *matched_idx, double **orien, double **org,
               long **ring_atom, char **AtomName, char *bseq, double *BPRS,
               long *pair_stat);
long bp_coplanar(long i, double d, double d2, double *txyz, double *txyz2, long n,
                 long *ddidx, double **bp_xyz);
void bp_context(long num_bp, double HELIX_CHG, double **bp_xyz, long **bp_order,
                long **end_list, long *num_ends, FILE * tfp);
void locate_helix(long num_bp, long **helix_idx, long num_ends, long *num_helix,
                  long **end_list, long **bp_order, long *bp_idx, long *helix_marker);
void first_step(long i, long **helix_idx, long *bp_idx, long *swapped, long **base_pairs,
                double **o3_p);
long chain1dir(long m, long n, long *swapped, long **base_pairs, double **o3_p);
double distance_ab(double **o3_p, long ia, long ib, long ipa, long ipb);
void get_ij(long m, long *swapped, long **base_pairs, long *n1, long *n2);
void get_bidx(long m, long *swapped, long *idx1, long *idx2);
double wcbp_zdir(long m, long n, long idxm1, long idxm2, long idxn1, long idxn2,
                 double **bp_xyz);
double wcbp_xang(long m, long n, double **bp_xyz);
long wc_bporien(long m, long n, long *swapped, long **base_pairs, double **bp_xyz,
                double **o3_p);
long check_o3dist(long m, long n, long *swapped, long **base_pairs, double **o3_p);
long is_linked(long i, long j, double **o3_p);
long check_schain(long m, long n, long *swapped, long **base_pairs, double **o3_p);
long check_others(long m, long n, long *swapped, long **base_pairs, double **o3_p,
                  double **bp_xyz);
void check_direction(long i, long **helix_idx, long *bp_idx, long *swapped,
                     long **base_pairs, double **o3_p, long *direction);
void check_strand2(long i, long **helix_idx, long *bp_idx, double **bp_xyz,
                   long *swapped, long **base_pairs, double **o3_p, long *direction,
                   FILE * tfp);
void check_rise(long i, long **helix_idx, long *bp_idx, long *swapped, long **base_pairs,
                double **bp_xyz, double **o3_p);
void five2three(long num_bp, long *num_helix, long **helix_idx, long *bp_idx,
                double **bp_xyz, long **base_pairs, double **o3_p, FILE * tfp);
void check_zdna(long *num_helix, long **helix_idx, long *bp_idx, double **bp_xyz,
                long **base_pairs, FILE * tfp);
void wc3(long *pair_k, char *wc);
void re_ordering(long num_bp, long **base_pairs, long *bp_idx, long *helix_marker,
                 long **helix_idx, double *BPRS, long *num_helix, double **o3_p,
                 char *bseq, long **seidx, char **ResName, char *ChainID,
                 long *ResSeq, char **Miscs);
void helix_info(long **helix_idx, long idx, FILE * fp);
void write_bestpairs(long num_bp, long **base_pairs, long *bp_idx, char *bseq,
                     long **seidx, char **AtomName, char **ResName, char *ChainID,
                     long *ResSeq, char **Miscs, double **xyz, double **orien,
                     double **org, long **htm_water);
void write_fpmst(double *morg, double *morien, FILE * rframe);
void write_helix(long num_helix, long **helix_idx, long *bp_idx, long **seidx,
                 char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                 char **Miscs, double **xyz, long **base_pairs, long **htm_water);
void no_basepairs(char *pdbfile, char *outfile, char *parfile);
void duplex(long num, long num_residue, char *bseq, long **seidx, long *RY,
            char **AtomName, char **ResName, char *ChainID, long *ResSeq, char **Miscs,
            double **xyz, long curves, long divide, long hetatm, long pairs,
            long detailed, long waters, char *pdbfile, char *outfile, char *parfile);
long find_bestpair(long nout, long **base_pairs, long num_residue, char *bseq,
                   long **seidx, long *RY, char **AtomName, double **xyz, long *idx,
                   double **orien, double **org, double **NC1xyz, long **ring_atom,
                   double *BPRS);
void col_helices(long num_helix, long **helix_idx, long *bp_idx, long **base_pairs,
                 long **seidx, char *pdbfile, char *ChainID, long *ResSeq);
void x3dna_input(long idx, long start_num, long end_num, long nbp, char *pdbfile,
                 char *outfile, char *parfile, long hetatm, long *bp_idx,
                 long *helix_marker, long **helix_idx, long **base_pairs, long **seidx,
                 char **ResName, char *ChainID, long *ResSeq, char **Miscs, char *bseq,
                 double *BPRS, long detailed);
void curves_input(long idx, long start_num, long end_num, long nbp, char *pdbfile,
                  char *outfile, char *parfile, long *bp_idx, long **base_pairs,
                  long zdna, long parallel);

/* frame_mol.c */
void frame_mol_usage(void);
void frame_mol_cmdline(int argc, char *argv[], char *parfile, char *inpfile,
                       char *outfile, long *side_view, long *helix_axis, long *ref_axis,
                       long *is_helical, long *n1, long *n2, long *bp_str, long *frame);
void read_frames(char *parfile, long *num_bp, char **bp_seq, double **org,
                 double **orien, long *ds);
void attach_helix(long nbpm1, double **orien, double **org, long *num, char **AtomName,
                  long *ibase, double **xyz, long *nbond, long **linkage);
void ref_duplex_block(long num_bp, char **bp_seq, double **org, double **orien,
                      char *outfile);
void pdbalc_frame(long ref_axis, long num_bp, double **org, double **orien,
                  char *inpfile, char *outfile);
void get_mstorg(long n1, long n2, long is_helical, double **org, double **orien,
                double *morg, double **mst);
void frame_pdbfile(long num, long side_view, char *inpfile, char *outfile, double *morg,
                   double **mst);
void frame_alcfile(long num_bp, long side_view, long helix_axis, long ref_axis,
                   char *inpfile, char *outfile, double **org, double **orien,
                   double *morg, double **mst);

/* get_part.c */
void getpart_usage(void);
void getpart_cmdline(int argc, char *argv[], char *inpfile, char *outfile, long *na,
                     long *protein, long *ligand, long *dwater, long *dhatom, long *bb,
                     long *msc);
long is_ringatom(char *aname, long resid);
void msc_setting(char *inpfile, char *outfile, long *res_type, char *ChainID,
                 long num_residue, long **seidx);
long chain_colors(char *chains, char **ccolors, char **pnl, double *NA_coilradius);
void check_color_name(char *chnname, char *colname, long mscol_num, char **mscol_name);
long molscript_colors(char **mscol);
long expt_chainID(long num_residue, long **seidx, long *res_type, char *ChainID,
                  char *Echains, char *chain_type);
void residue_wtype(long num_residue, long **seidx, char **ResName, char **AtomName,
                   double **xyz, char **Miscs, long *res_type);
long contain_metal(long ib, long ie, char **AtomName);
void molauto_protein(char *chains, char **ccolors, char *Echains, char *chain_type,
                     char **pnl, FILE * fp);
void molscript_na(char *chains, char **ccolors, long num_Echains, char *Echains,
                  char *chain_type, char **pnl, double NA_coilradius, FILE * fp);
void write_selected(char *inpfile, char *outfile, long num_residue, long **seidx,
                    long *res_type, char **AtomName, char **ResName, char *ChainID,
                    long *ResSeq, double **xyz, char **Miscs, long bb, long dwater,
                    long na, long protein, long ligand, long dhatom);

/* o1p_o2p.c */

/* pdb2img.c */
void pdb2img_usage(void);
void pdb2img_cmdline(int argc, char *argv[], char *inpfile, char *outfile,
                     double *scale_factor, long *opts);
void pdb2alc(char *inpfile, char *alcfile, long peptide);
void base_blks(long num_residue, long *RY, double **orien, double **org, char *bseq,
               char *BDIR, char *alcfile);
void peptide_info(long num_residue, long **seidx, char **AtomName, char **ResName,
                  char *ChainID, long *ResSeq, char **Miscs, double **xyz, long *resid,
                  long *cidx, long *mchain);
void peptide_frame(long num_residue, char *BDIR, long *resid, long *mchain, double **xyz,
                   double **orien, double **org);
void peptide_blks(long num_residue, char *BDIR, long *cidx, double **orien, double **org,
                  char *alcfile);

/* r3d_atom.c */
void r3d_atom_usage(void);
void r3d_atom_cmdline(int argc, char *argv[], char *hbfile, char *lkgfile,
                      char *inpfile, char *outfile, double *rodrad, double *ballrad,
                      double *scale_factor, long *opts);
void stick_model(long nobase, char *lkgfile, long num, long *RY, long *idx,
                 char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                 char **Miscs, double **xyz, double rodrad, long *colidx,
                 double **colrgb, FILE * fp);
void r3d_input(char *hbfile, char *lkgfile, char *inpfile, char *outfile, double rodrad,
               double ballrad, double scale_factor, long *opts);
void gray_image(long is_gray, double **atom_col, double **base_col, double *hb_col);
void col_residue(long num, long num_residue, long *RY, long **seidx, long *idx,
                 double **atom_col, double **base_col, long *ibase, long by_residue,
                 long *colidx, double **colrgb);
void render_hbonds(long num, long num_residue, char *bseq, long **seidx, long *RY,
                   char **AtomName, char **ResName, char *ChainID, long *ResSeq,
                   char **Miscs, double **xyz, long pwise, double *width3,
                   double rodrad, char *hbfile, double *hb_col, FILE * fp);

/* reb_fncs.c */
void rebuild_usage(void);
void rebuild_cmdline(int argc, char *argv[], char *options[], char *user_input,
                     char *inpfile, char *outfile, long *neg_xdir);
void read_par(char *inpfile, long *num_bp, long *is_single, long *is_helical, long *xdir,
              long *parallel, char **bp_seq, double **bp_par, double **step_par,
              long *pidx);
void link_o3_p(long num_residue, long **seidx, char **AtomName, double **xyz,
               char *ChainID, long **connect);
void atom_lkg(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
              double **xyz, char *outfile);
void atomic_pdb1(long num_bp, long num_atoms, long is_helical, long xdir, char **bp_seq,
                 double **step_par, char *BDIR, char *outfile);
void atomic_pdb2(long parallel, long num_bp, long num_atoms, long is_helical, long xdir,
                 char **bp_seq, double **bp_par, double **step_par, char *BDIR,
                 char *outfile);
void base_c1_atoms(char **AtomName, long ib, long ie, long *num_batom, long *batom);
void extract_base_atoms(long num, char **AtomName, double **xyz, long *bnum,
                        char **bAtomName, double **bxyz);
void atomic_base_p(long parallel, long num_bp, long num_atoms, long is_helical,
                   long xdir, char **bp_seq, double **bp_par, double **step_par,
                   char *BDIR, long *pidx, char *outfile);
long num_PDB_atoms(long num_bp, long is_single, char **bp_seq, char *BDIR);
void set_bp_pdb(char *fname1, char *fname2, double *param, long xdir, long *num1,
                long *num2, char **AtomName1, char **AtomName2, double **xyz1,
                double **xyz2, char ap);
void base_fname(char bname, char *BDIR, char *fname);
void block_alc1(long num_bp, long is_single, long is_helical, long xdir, char **bp_seq,
                double **step_par, char *BDIR, char *outfile);
void block_alc2(long num_bp, long is_helical, long xdir, char **bp_seq, double **bp_par,
                double **step_par, char *BDIR, char *outfile);
void set_bp_alc(char *fname1, char *fname2, double *param, long xdir, long num,
                long nbond, double **bp_xyz, char ap);
void xbpfunc(double *param, double **orien, double **mst, double *pos, double *mpos);
void xhelfunc(double *param, double **orien, double **mst, double *pos, double *mpos);
void print_ref_frames(FILE * fp, double *pos_next, double **orien_next);

/* rebuild.c */

/* regular_dna.c */
void regular_dna_usage(void);
void regular_dna_cmdline(int argc, char *argv[], char *options[2], char *user_input,
                         char *outfile);
void get_six_pars(long is_bp, double *pars);

/* rotate_mol.c */
void rotmol_usage(void);
void rotate_mol_cmdline(long argc, char *argv[], char *rotfile, char *inpfile,
                        char *outfile, long *fnum, long *to_center, long *pmi,
                        long *to_frame);
void molscr_rotmat(char *rotfile, double **rotmat);
void read_rotmat(char *rotfile, double **rotmat);
void read_bp1frame(char *rotfile, double **rotmat);
void rot_pdbfile(long num, char *inpfile, char *outfile, long to_center, double **rotmat,
                 long pmi, long to_frame);
void pdb_pmi(long num, char **AtomName, char **ResName, char *ChainID, long *ResSeq,
             double **xyz, char **Miscs, long pmi, double **nxyz, double **refmat,
             double **rotmat, long *nnum);
void rot_alcfile(char *inpfile, char *outfile, long to_center, double **rotmat, long pmi,
                 long to_frame);
void cmb_transformation(double **refmat, double **rotmat);
void write_rotmat(double **rotmat);
void base1_check(long num_residue, long *RY, double **org, double **rotmat);
void pmi_rotmat(long num, double **xyz, double **rotmat);
void pmi_view(double **rotmat, double **tmprot, long c1, long c2, long c3);
void reorient_xyz(long num, double **xyz, double *orien1, double *org1, double **refmat);
void align_str(double **rotmat, long num, double **xyz, double **nxyz);

/* stack2img.c */
void stack2img_usage(void);
void stack2img_cmdline(int argc, char *argv[], char *hbfile, char *lkgfile, char *inpfile,
                       char *outfile, double *scale_factor, double *ballrad, long *opts);
void get_pair_num(char *inpfile, long **pair_num);
void process_stack(char *hbfile, char *lkgfile, char *inpfile, char *outfile,
                   double scale_factor, double ballrad, long *opts);
void stack2fig(long nobj, long *idx, long *aidx, long *depth, long **allobj,
               double **xyz, long **linkage, long **hb_linkage, long **ring_atom,
               long *ibase, long is_color, long **pair_num, long label_ring,
               long *ResSeq, double ballrad, FILE * fp);
void stack2ps(long nobj, long *idx, long *aidx, long **allobj, double **xyz,
              long **linkage, long **hb_linkage, long **ring_atom, long *ibase,
              long is_color, long **pair_num, long label_ring, long *ResSeq,
              double ballrad, FILE * fp);
void stack2r3d(long *aidx, long **allobj, long num, double **xyz, long *ResSeq,
               long nbond, long **linkage, long num_ring, long **ring_atom,
               long num_hbond, long **hb_linkage, long **pair_num, long num_ball,
               long num_residue, long *ibase, char *bseq, double scale_factor,
               double ballrad, long *opts, FILE * fp);
void get_stack_objs(long num_ball, double **xyz, long **seidx, long nbond,
                    long **linkage, long num_residue, long **ring_atom, long num_hbond,
                    long **hb_linkage, long **allobj);

/* std_base.c */
void fit2standard_base(long num, char **AtomName, double **xyz, char bname, long ry,
                       double **nxyz);
void std_base_usage(void);
void std_base_cmdline(int argc, char *argv[], char *inpfile, char *outfile, char *bname,
                      double *c1_n_org, double *dist_n2org, double *n_org_x,
                      long *fit_std);

/* step_hel.c */
