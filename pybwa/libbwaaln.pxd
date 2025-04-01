# cython: language_level=3

from libc.stdint cimport uint8_t, uint64_t, uint16_t, uint32_t, int64_t, int32_t
from libc.stdio cimport FILE

cdef extern from "libbwaaln_utils.h":
    bam1_t **bwa_aln_and_samse(const bntseq_t *bns, bwt_t *const bwt, uint8_t *pac, sam_hdr_t *h, int n_seqs, bwa_seq_t *seqs,
                     const gap_opt_t *opt, int max_hits, int with_md)

cdef extern from "bntseq.h":
    unsigned char nst_nt4_table[256]
    int bns_cnt_ambi(const bntseq_t *bns, int64_t pos_f, int len, int *ref_id)

cdef extern from "bwt.h":
    ctypedef struct bwt_t:
        int sa_intv

cdef extern from "bwtaln.h":
    const int BWA_TYPE_NO_MATCH
    const int BWA_TYPE_MATESW

    const int BWA_MODE_GAPE
    const int BWA_MODE_COMPREAD
    const int BWA_MODE_LOGGAP
    const int BWA_MODE_NONSTOP

    const int FROM_M
    const int FROM_I
    const int FROM_D
    const int FROM_S

    int __cigar_op(uint16_t __cigar)
    int __cigar_len(uint16_t __cigar)

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
    void gap_print_opt(const gap_opt_t *opt)

    void seq_reverse(int len, unsigned char *seq, int is_comp)

    ctypedef struct bwt_aln1_t:
        pass

cdef extern from "bntseq.h":
    ctypedef  struct bntann1_t:
        int64_t offset
        int32_t len
        char *name

    ctypedef struct bntseq_t:
        int64_t l_pac
        bntann1_t *anns
        FILE * fp_pac

    bntseq_t * bns_restore(const char * prefix)
    void bns_destroy(bntseq_t *bns)

cdef extern from "htslib/kstring.h":
    ctypedef struct kstring_t:
        size_t l, m
        char *s

cdef extern void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi)

cdef extern from "bwase.h":
    int64_t pos_end(const bwa_seq_t *p)
    void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, unsigned char *_pacseq, int with_md)
    char *bwa_cal_md1(int n_cigar, uint16_t *cigar, int len, uint64_t pos, unsigned char *seq, uint64_t l_pac, unsigned char *pacseq, kstring_t *str, int *_nm)
    void bwase_initialize()

cdef extern from "bwtaln.h":
    ctypedef struct bwt_multi1_t:
        uint32_t n_cigar
        uint32_t gap
        uint32_t mm
        uint32_t strand
        int ref_shift
        uint64_t pos
        uint16_t *cigar
    ctypedef struct bwa_seq_t:
        char *name
        uint8_t *seq
        uint8_t *rseq
        uint8_t *qual
        uint32_t len
        uint32_t strand
        uint32_t type
        uint32_t n_mm
        uint32_t n_gapo
        uint32_t n_gape
        uint32_t mapQ
        int score
        int clip_len
        int n_aln
        bwt_aln1_t *aln
        int n_multi;
        bwt_multi1_t *multi;
        int n_occ;
        uint64_t pos
        uint64_t c1
        uint64_t c2
        uint16_t *cigar
        int n_cigar
        int tid
        uint32_t full_len
        uint32_t nm
        char *md


    void bwa_free_read_seq(int n_seqs, bwa_seq_t *seqs)

    void bwa_cal_sa_reg_gap(int tid, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int free_seqs)


# from bwase.c
cdef extern int64_t pos_end_multi(const bwt_multi1_t *p, int len)
cdef extern bam1_t* bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2, kstring_t *str, sam_hdr_t *h)

cpdef bint _set_bwa_aln_verbosity(int level)


cdef extern from "bwa.h":
    void bwa_format_sam_hdr(const bntseq_t *bns, const char *rg_line, kstring_t *str)


cdef extern from "htslib/sam.h":
    ctypedef struct bam1_t:
        pass

    void bam_destroy1(bam1_t *b)

    ctypedef struct sam_hdr_t:
        size_t l_text
        char *text

    sam_hdr_t *sam_hdr_parse(size_t l_text, const char *text)


cdef class BwaAlnOptions:
    cdef gap_opt_t* gap_opt(self)
    cdef gap_opt_t * _delegate
    cdef public object _max_hits
    cdef public object _with_md