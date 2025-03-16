# cython: language_level=3

cdef extern from "bwtaln.h":
    const int BWA_MODE_GAPE
    const int BWA_MODE_COMPREAD
    const int BWA_MODE_LOGGAP
    const int BWA_MODE_NONSTOP

    ctypedef struct gap_opt_t:
        int trim_qual
        int s_mm
        int s_gapo
        int s_gape
        int mode # bit 24-31 are the barcode length
        int indel_end_skip
        int max_del_occ
        int max_entries
        float fnr
        int max_diff
        int max_gapo
        int max_gape
        int max_seed_diff
        int seed_len
        int n_threads
        int max_top2
        int trim_qual
        int sam
        char *rg_line
        int n_occ
        int interactive_mode
        int with_md

    gap_opt_t *gap_init_opt()


cdef class BwaAlnOptions:
    cdef gap_opt_t * gap_opt(self)
    cdef gap_opt_t * _delegate
    cdef public object _max_hits
    cdef public object _with_md
