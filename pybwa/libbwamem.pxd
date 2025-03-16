# cython: language_level=3

from libc.stdint cimport uint8_t, int64_t, int32_t, int8_t, uint32_t
from libc.stdio cimport FILE
from pybwa.libbwamemopt cimport mem_opt_t


cdef extern from "bwa.h":
    ctypedef struct bseq1_t:
        int l_seq, id
        char *name
        char *comment
        char *seq
        char *qual
        char *sam

cdef extern from "libbwamem_utils.h":
    ctypedef  struct mem_alns_t:
        size_t n, m
        mem_aln_t*a
    mem_alns_t * mem_process_seqs_alt(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns,
                                  const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs,
                                  const mem_pestat_t *pes0)

cdef extern from "limits.h":
    cdef int INT_MAX

cdef extern from "bwt.h":
    ctypedef struct bwt_t:
        int sa_intv

cdef extern from "bntseq.h":
    ctypedef  struct bntann1_t:
        int64_t offset
        int32_t len
        char *name
        char *anno

    ctypedef struct bntseq_t:
        int64_t l_pac
        bntann1_t *anns
        FILE * fp_pac

    unsigned char nst_nt4_table[256]

cdef extern from "kstring.h":
    ctypedef struct kstring_t:
        size_t l, m
        char *s

cdef extern from "bwamem.h":
    ctypedef struct mem_pestat_t:
        int low, high   # lower and upper bounds within which a read pair is considered to be properly paired
        int failed      # non-zero if the orientation is not supported by sufficient data
        double avg, std # mean and stddev of the insert size distribution

    ctypedef struct mem_alnreg_t:
        int score  # best local SW score
        int secondary  # index of the parent hit shadowing the current hit; <0 if primary
        int n_comp  # number of sub-alignments chained together
        int is_alt

    ctypedef struct mem_alnreg_v:
        size_t n, m
        mem_alnreg_t *a

    ctypedef struct mem_aln_t:
        int64_t pos     # forward strand 5'-end mapping position
        int rid         # reference sequence index in bntseq_t; <0 for unmapped
        int flag        # extra flag
        uint32_t is_rev  # is_rev: whether on the reverse strand;
        uint32_t is_alt
        uint32_t mapq   # mapq: mapping quality;
        uint32_t NM  # NM: edit distance
        int n_cigar;     # number of CIGAR operations
        uint32_t *cigar; # CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
        char *XA
        int score, sub, alt_sc;

    mem_opt_t *mem_opt_init()
    mem_alnreg_v mem_align1(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq)
    void bwa_fill_scmat(int a, int b, int8_t mat[25])
    mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq, const mem_alnreg_t *ar)

# from bwamem.c
cdef extern void mem_reorder_primary5(int T, mem_alnreg_v *a)
cdef extern void add_cigar(const mem_opt_t *opt, mem_aln_t *p, kstring_t *str, int which)

# from bwamem_extra.c
cdef extern char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, mem_alnreg_v *a, int l_query, const char *query);
