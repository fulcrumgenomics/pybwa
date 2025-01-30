# cython: language_level=3
from pathlib import Path

from libc.stdint cimport int64_t, int32_t, uint8_t
from libc.stdio cimport FILE

cdef extern from "bwa.h":
    char * bwa_idx_infer_prefix(const char * hint)
    ctypedef struct bwaidx_t:
        bwt_t    *bwt # FM-index
        bntseq_t *bns # information on the reference sequences
        uint8_t  *pac # the actual 2-bit encoded reference sequences with 'N' converted to a random base

    bwaidx_t * bwa_idx_load(const char * hint, int which)
    void bwa_idx_destroy(bwaidx_t *idx)
    int BWA_IDX_BWT
    int BWA_IDX_BNS
    int BWA_IDX_PAC
    int bwa_idx_build(const char *fa, const char *prefix, int algo_type, int block_size)


cdef extern from "bwt.h":
    ctypedef struct bwt_t:
        int sa_intv

    bwt_t *bwt_restore_bwt(const char *fn)
    void bwt_restore_sa(const char *fn, bwt_t *bwt);
    void bwt_destroy(bwt_t *bwt)

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

cdef bytes force_bytes(object s)
cdef bytes force_bytes_with(object s, encoding: str | None = *, errors: str | None = *)

cdef class BwaIndex:
    """Contains the index and nucleotide sequence for Bwa"""
    cdef bwaidx_t *_delegate
    cdef public object header
    cdef bwt_t *bwt(self)
    cdef bntseq_t *bns(self)
    cdef uint8_t *pac(self)
    cdef _load_index(self, prefix, mode)
