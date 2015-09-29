#inspired by sage/src/sage/numerical/backends/glpk_backend.pxd

cdef extern from "glpk.h":
    ctypedef struct glp_prob "glp_prob":
        pass
    ctypedef struct glp_iocp "glp_iocp":
        int msg_lev
        int br_tech
        int bt_tech
        int pp_tech
        int fp_heur
        int gmi_cuts
        int mir_cuts
        int cov_cuts
        int clq_cuts
        double tol_int
        double tol_obj
        double mip_gap
        int tm_lim
        int out_frq
        int out_dly
        int presolve
        int binarize
    ctypedef struct glp_smcp "glp_smcp":
        int msg_lev
        int meth
        int pricing
        int r_test
        double tol_bnd
        double tol_dj
        double tol_piv
        double obj_ll
        double obj_ul
        int it_lim
        int tm_lim
        int out_frq
        int out_dly
        int presolve
    glp_iocp * new_glp_iocp "new glp_iocp" ()
    void glp_init_iocp(glp_iocp *)
    void glp_init_smcp(glp_smcp *)
    glp_prob * glp_create_prob()
    void glp_set_prob_name(glp_prob *, char *)
    void glp_set_obj_dir(glp_prob *, int)
    void glp_add_rows(glp_prob *, int)
    void glp_add_cols(glp_prob *, int)
    void glp_del_rows(glp_prob *, int, int *)
    void glp_set_row_name(glp_prob *, int, char *)
    void glp_set_col_name(glp_prob *, int, char *)
    void glp_set_row_bnds(glp_prob *, int, int, double, double)
    void glp_set_col_bnds(glp_prob *, int, int, double, double)
    void glp_set_obj_coef(glp_prob *, int, double)
    void glp_load_matrix(glp_prob *, int, int *, int *, double *)
    int glp_simplex(glp_prob *, glp_smcp *)
    int glp_exact(glp_prob *, glp_smcp *)  # requires gmp
    int glp_intopt(glp_prob *, glp_iocp *)
    void glp_std_basis(glp_prob *)
    void glp_delete_prob(glp_prob *)
    double glp_get_col_prim(glp_prob *, int)
    double glp_get_obj_val(glp_prob *)
    double glp_get_col_dual(glp_prob *, int)
    double glp_get_row_dual(glp_prob *, int)
    int glp_print_ranges(glp_prob *lp, int,int, int, char *fname)
    int glp_get_num_rows(glp_prob *)
    int glp_get_num_cols(glp_prob *)
    int glp_get_num_int(glp_prob *)
    double glp_mip_col_val(glp_prob *, int)
    double glp_mip_obj_val(glp_prob *)
    void glp_set_col_kind(glp_prob *, int, int)
    int glp_write_mps(glp_prob *lp, int fmt, void *parm, char *fname)
    int glp_write_lp(glp_prob *lp, void *parm, char *fname)
    int glp_write_prob(glp_prob *P, int flags, char *fname)
    int glp_read_prob(glp_prob *P, int flags, char *fname)

    void glp_set_prob_name(glp_prob *lp, char *name)
    void glp_set_obj_name(glp_prob *lp, char *name)
    void glp_set_row_name(glp_prob *lp, int i, char *name)
    void glp_set_col_name(glp_prob *lp, int i, char *name)

    double glp_get_row_ub(glp_prob *lp, int i)
    double glp_get_row_lb(glp_prob *lp, int i)

    double glp_get_col_ub(glp_prob *lp, int i)
    double glp_get_col_lb(glp_prob *lp, int i)
    void glp_set_col_ub(glp_prob *lp, int i, double value)
    void glp_set_col_lb(glp_prob *lp, int i, double value)


    void glp_create_index(glp_prob *P)
    int glp_find_row(glp_prob *P, const char *name)
    int glp_find_col(glp_prob *P, const char *name)
    void glp_delete_index(glp_prob *P)

    double glp_get_col_lb(glp_prob *lp, int i)
    double glp_get_col_ub(glp_prob *lp, int i)

    void glp_scale_prob(glp_prob *lp, int flags)
    void glp_unscale_prob(glp_prob *lp)

    int glp_get_prim_stat(glp_prob *lp)
    int glp_get_status(glp_prob *lp)
    int glp_mip_status(glp_prob *lp)
    int glp_get_num_nz(glp_prob *lp)
    int glp_set_mat_row(glp_prob *lp, int, int, int *, double * )
    int glp_set_mat_col(glp_prob *lp, int, int, int *, double * )
    int glp_get_mat_row(glp_prob *lp, int, int *, double * )
    int glp_get_mat_col(glp_prob *lp, int, int *, double * )
    double glp_get_row_ub(glp_prob *lp, int)
    double glp_get_row_lb(glp_prob *lp, int)
    int glp_get_col_kind(glp_prob *lp, int)
    double glp_get_obj_coef(glp_prob *lp, int)
    int glp_get_obj_dir(glp_prob *lp)
    void glp_copy_prob(glp_prob *dst, glp_prob *src, int names)

    const char *glp_version()

    # output redirection
    int glp_term_out(int flag)
    void glp_term_hook(int (*func)(void *info, const char *s), void *info)

    int glp_warm_up(glp_prob *P)
    void glp_adv_basis(glp_prob *P, int flags)

    # constants

    # constants for smcp control

    int GLP_MSG_OFF
    int GLP_MSG_ERR
    int GLP_MSG_ON
    int GLP_MSG_ALL

    int GLP_PRIMAL
    int GLP_DUALP
    int GLP_DUAL

    int GLP_PT_STD
    int GLP_PT_PSE

    int GLP_RT_STD
    int GLP_RT_HAR

    double DBL_MAX

    int INT_MAX

    int GLP_ON
    int GLP_OFF

    # constants for scaling the problem
    int GLP_SF_AUTO
    int GLP_SF_GM
    int GLP_SF_EQ
    int GLP_SF_2N
    int GLP_SF_SKIP

    # constants for iocp control, not already in simplex

    int GLP_BR_FFV
    int GLP_BR_LFV
    int GLP_BR_MFV
    int GLP_BR_DTH
    int GLP_BR_PCH

    int GLP_BT_DFS
    int GLP_BT_BFS
    int GLP_BT_BLB
    int GLP_BT_BPH

    int GLP_PP_NONE
    int GLP_PP_ROOT
    int GLP_PP_ALL

    # error codes
    int GLP_EBADB
    int GLP_ESING
    int GLP_ECOND
    int GLP_EBOUND
    int GLP_EFAIL
    int GLP_EOBJLL
    int GLP_EOBJUL
    int GLP_EITLIM
    int GLP_ETMLIM
    int GLP_ENOPFS
    int GLP_ENODFS
    int GLP_EROOT
    int GLP_ESTOP
    int GLP_EMIPGAP
    int GLP_ENOFEAS
    int GLP_ENOCVG
    int GLP_EINSTAB
    int GLP_EDATA
    int GLP_ERANGE


    int GLP_UNDEF
    int GLP_OPT
    int GLP_FEAS
    int GLP_NOFEAS
    int GLP_INFEAS
    int GLP_UNBND

    # other constants

    int GLP_MAX
    int GLP_MIN
    int GLP_UP
    int GLP_FR
    int GLP_DB
    int GLP_FX
    int GLP_LO
    int GLP_CV
    int GLP_IV
    int GLP_BV
    int GLP_MPS_DECK
    int GLP_MPS_FILE

    int GLP_MSG_DBG
