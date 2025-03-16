# cython: language_level=3

from libc.stdint cimport uint64_t, int8_t

cdef extern from "bwamem.h":
    const int MEM_F_PE
    const int MEM_F_NOPAIRING
    const int MEM_F_ALL
    const int MEM_F_NO_MULTI
    const int MEM_F_NO_RESCUE
    const int MEM_F_REF_HDR
    const int MEM_F_SOFTCLIP
    const int MEM_F_SMARTPE
    const int MEM_F_PRIMARY5
    const int MEM_F_KEEP_SUPP_MAPQ
    const int MEM_F_XB

    ctypedef struct mem_opt_t:
        int a, b                # match score and mismatch penalty
        int o_del, e_del
        int o_ins, e_ins
        int pen_unpaired        # phred-scaled penalty for unpaired reads
        int pen_clip5,pen_clip3 # clipping penalty. This score is not deducted from the DP score.
        int w                   # band width
        int zdrop               # Z-dropoff
        uint64_t max_mem_intv
        int T                  # output score threshold only affecting output
        int flag               # see MEM_F_* macros
        int min_seed_len       # minimum seed length
        int min_chain_weight
        int max_chain_extend
        float split_factor     # split into a seed if MEM is longer than min_seed_len*split_factor
        int split_width        # split into a seed if its occurence is smaller than this value
        int max_occ            # skip a seed if its occurence is larger than this value
        int max_chain_gap      # do not chain seed if it is max_chain_gap-bp away from the closest seed
        int n_threads          # number of threads
        int chunk_size         # process chunk_size-bp sequences in a batch
        float mask_level       # regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
        float drop_ratio       # drop a chain if its seed coverage is below drop_ratio times the seed coverage of a better chain overlapping with the small chain
        float XA_drop_ratio    # when counting hits for the XA tag, ignore alignments with score < XA_drop_ratio * max_score only effective for the XA tag
        float mask_level_redun
        float mapQ_coef_len
        int mapQ_coef_fac
        int max_ins            # when estimating insert size distribution, skip pairs with insert longer than this value
        int max_matesw         # perform maximally max_matesw rounds of mate-SW for each end
        int max_XA_hits, max_XA_hits_alt # if there are max_hits or fewer, output them all
        int8_t mat[25]         # scoring matrix mat[0] == 0 if unset

    mem_opt_t *mem_opt_init()

    void bwa_fill_scmat(int a, int b, int8_t mat[25])


cdef class BwaMemOptions:
    cdef public object _finalized
    cdef public object _ignore_alt
    cdef public object _mode
    cdef mem_opt_t * _options
    cdef mem_opt_t * _options0
    cdef mem_opt_t* mem_opt(self)
