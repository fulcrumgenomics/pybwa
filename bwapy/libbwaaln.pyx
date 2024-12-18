# cython: language_level=3

from pathlib import Path
from typing import List

from libc.stdint cimport uint8_t
from libc.stdio cimport SEEK_SET
from libc.stdlib cimport calloc, free
from libc.string cimport strncpy
from pysam import FastxRecord, AlignedSegment
from bwapy.libbwaindex cimport force_bytes
from bwapy.libbwaindex cimport BwaIndex


__all__ = [
    "BwaAlnOptions",
    "BwaAln",
]

cdef class BwaAlnOptions:
    """The container for options for [`BwaAln`][bwapy.BwaAln].

    Use [`BwaAlnOptionsBuilder`][bwapy.BwaAlnOptionsBuilder] to use and set custom options.
    """
    cdef gap_opt_t * _delegate

    _max_hits: int

    def __init__(self):
        self._max_hits = 3
        self._cinit()

    cdef _cinit(self):
        self._delegate = gap_init_opt()

    def __dealloc__(self):
        free(self._delegate)

    cdef gap_opt_t* gap_opt(self):
        return self._delegate

    property max_mismatches:
        def __get__(self):
            return self._delegate.s_mm
        def __set__(self, value: int):
            self._delegate.s_mm = value

    property max_gap_opens:
        """bwa aln -o <int>"""
        def __get__(self):
            return self._delegate.max_gapo
        def __set__(self, value: int):
           self._delegate.max_gapo = value

    property max_gap_extensions:
        """bwa aln -e <int>"""
        def __get__(self):
            return self._delegate.max_gapo
        def __set__(self, value: int):
            self._delegate.max_gape = value
            if self._delegate.max_gape > 0:
                self._delegate.mode &= ~BWA_MODE_GAPE

    property min_indel_to_end_distance:
        """bwa aln -i <int>"""
        def __get__(self):
            return self._delegate.indel_end_skip
        def __set__(self, value: int):
           self._delegate.indel_end_skip = value

    property max_occurrences_for_extending_long_deletion:
        """bwa aln -d <int>"""
        def __get__(self):
            return self._delegate.max_del_occ
        def __set__(self, value: int):
           self._delegate.max_del_occ = value

    property seed_length:
        """bwa aln -l <int>"""
        def __get__(self):
            return self._delegate.seed_len
        def __set__(self, value: int):
           self._delegate.seed_len = value

    property max_mismatches_in_seed:
        """bwa aln -k <int>"""
        def __get__(self):
            return self._delegate.max_seed_diff
        def __set__(self, value: int):
           self._delegate.max_seed_diff = value

    property mismatch_penalty:
        """bwa aln -M <int>"""
        def __get__(self):
            return self._delegate.s_mm
        def __set__(self, value: int):
           self._delegate.s_mm = value

    property gap_open_penalty:
        """bwa aln -O <int>"""
        def __get__(self):
            return self._delegate.s_gapo
        def __set__(self, value: int):
           self._delegate.s_gapo = value

    property gap_extension_penalty:
        """bwa aln -E <int>"""
        def __get__(self):
            return self._delegate.s_gape
        def __set__(self, value: int):
           self._delegate.s_gape = value

    property stop_at_max_best_hits:
        """bwa aln -R <int>"""
        def __get__(self):
            return self._delegate.max_top2
        def __set__(self, value: int):
           self._delegate.max_top2 = value

    property max_hits:
        """bwa samse -n <int>"""
        def __get__(self):
            return self._max_hits
        def __set__(self, value: int):
           self._max_hits = value

    property log_scaled_gap_penalty:
        """bwa aln -L"""
        def __get__(self):
            return self._delegate.mode & BWA_MODE_LOGGAP > 0
        def __set__(self, value: int):
            if value:
                self._delegate.mode |= BWA_MODE_LOGGAP
            else:
                self._delegate.mode &= ~BWA_MODE_LOGGAP


cdef class BwaAln:
    """The class to align reads with `bwa aln`."""

    cdef BwaIndex _index

    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None):
        if prefix is not None:
            assert Path(prefix).exists()
            self._index = BwaIndex(prefix=prefix)
        elif index is not None:
            self._index = index
        else:
            raise Exception("Either prefix or index must be given")

        bwase_initialize()

    # TODO: a list of records...
    def align(self, opt: BwaAlnOptions, queries: List[FastxRecord]) -> List[AlignedSegment]:
        """Align one or more queries with `bwa aln`.

        Args:
            opt: the alignment options
            queries: the queries to align

        Returns:
            one alignment per query
        """
        return self._calign(opt,  queries)

    cdef _copy_seq(self, q: FastxRecord, bwa_seq_t* s):
        seq_len = len(q.sequence)
        s.len = seq_len
        s.clip_len = seq_len
        s.full_len = seq_len
        s.seq = <uint8_t *> calloc(sizeof(uint8_t), seq_len + 1)
        s.rseq = <uint8_t *> calloc(sizeof(uint8_t), seq_len + 1)
        s.qual = <uint8_t *> calloc(sizeof(uint8_t), seq_len + 1)

        # use seq_reverse from bwaseqio.c
        for i, base in enumerate(q.sequence):
            s.seq[i] = nst_nt4_table[ord(base)]
            s.qual[i] = 40  # FIXME: parameterize
        s.seq[seq_len] = b'\0'
        s.qual[seq_len] = b'\0'
        seq_reverse(seq_len, s.seq,
                    0)  #  // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
        seq_reverse(seq_len, s.rseq, 0)

        s.name = <char *> calloc(sizeof(char), len(q.name) + 1)
        strncpy(s.name, force_bytes(q.name), len(q.name))
        s.name[len(q.name)] = b'\0'

    cdef _build_alignment(self, query: FastxRecord, bwa_seq_t *seq, kstring_t *kstr):
        cdef int reference_id
        cdef int nm

        # make a default, unmapped, empty record
        rec = AlignedSegment(header=self._index.header)
        rec.query_name = query.name
        rec.reference_id = -1
        rec.reference_start = -1
        rec.mapping_quality = 0
        rec.is_paired = False
        rec.is_read1 = True
        rec.is_read2 = False
        rec.is_qcfail = False
        rec.is_duplicate = False
        rec.is_secondary = False
        rec.is_supplementary = False
        rec.is_unmapped = True
        rec.query_sequence = query.sequence
        rec.query_qualities = query.quality

        if seq.type == BWA_TYPE_NO_MATCH:  # unmapped read
            # TODO: custom bwa tags: RG, BC, XC
            return rec

        ref_len_in_alignment = pos_end(seq) - seq.pos

        # if on the reverse strand, reverse the query sequence and qualities
        rec.query_sequence = query.sequence[::-1]
        if query.quality is not None:
            rec.query_qualities = query.quality[::-1]

        # reference id
        nn = bns_cnt_ambi(self._index.bns(), seq.pos, ref_len_in_alignment, &reference_id)
        rec.reference_id = reference_id

        # make this unmapped if we map off the end of the contig
        rec.is_unmapped = seq.pos + ref_len_in_alignment - self._index.bns().anns[
            rec.reference_id].offset > self._index.bns().anns[rec.reference_id].len

        # strand, reference start, and mapping quality
        if seq.strand:
            rec.is_reverse = True
        rec.reference_start = seq.pos - self._index.bns().anns[rec.reference_id].offset  # 0-based
        rec.mapping_quality = seq.mapQ

        # cigar
        cigar = ""
        if seq.cigar:
            for j in range(seq.n_cigar):
                cigar_len = __cigar_len(seq.cigar[j])
                cigar_op = "MIDS"[__cigar_op(seq.cigar[j])]
                cigar = f"{cigar}{cigar_len}{cigar_op}"
        elif seq.type != BWA_TYPE_NO_MATCH:
            cigar = f"{seq.len}M"
        rec.cigarstring = cigar

        # # tags
        if seq.type != BWA_TYPE_NO_MATCH:
            attrs = dict()
            attrs["MD"] = f"{seq.md}"
            attrs["NM"] = f"{seq.nm}"
            rec.set_tags(list(attrs.items()))
        # # TODO:the custom bwa tags: XT, NM, XN, SM, AM, X0, X1, XM, XO, XG, XA, HN

        return rec

    cdef _calign(self, opt: BwaAlnOptions, queries: List[FastxRecord]):
        cdef bwa_seq_t* seqs
        cdef bwa_seq_t* s
        cdef char* s_char
        cdef kstring_t* kstr
        cdef gap_opt_t* gap_opt

        gap_opt = opt.gap_opt()

        kstr = <kstring_t*>calloc(sizeof(kstring_t), 1)

        # copy FastqProxy into bwa_seq_t
        num_seqs = len(queries)
        seqs = <bwa_seq_t*>calloc(sizeof(bwa_seq_t), num_seqs)
        for i in range(num_seqs):
            self._copy_seq(queries[i], &seqs[i])
            seqs[i].tid = -1

        # this is `bwa aln`, and the rest is `bwa samse`
        bwa_cal_sa_reg_gap(0, self._index.bwt(), num_seqs, seqs, gap_opt)

        # create the full alignment
        for i in range(num_seqs):
            s = &seqs[i]
            # bwa_cal_sa_reg_gap frees name, seq, rseq, and qual, so add them back in again
            self._copy_seq(queries[i], s)
            bwa_aln2seq_core(s.n_aln, s.aln, s, 1, opt.max_hits)

        # # calculate the genomic position given the suffix array offsite
        bwa_cal_pac_pos_with_bwt(self._index.bns(), num_seqs, seqs, gap_opt.max_diff, gap_opt.fnr, self._index.bwt())

        # refine gapped alignment
        bwa_refine_gapped(self._index.bns(), num_seqs, seqs, self._index.pac())

        # create the AlignedSegment from FastxRecord and bwa_seq_t.
        recs = [
            self._build_alignment(query=queries[i], seq=&seqs[i], kstr=kstr)
            for i in range(num_seqs)
        ]

        bwa_free_read_seq(num_seqs, seqs)
        free(kstr.s)
        free(kstr)

        return recs