# cython: language_level=3

from pathlib import Path
from typing import List
from fgpyo.sequence import reverse_complement

from libc.stdint cimport uint8_t
from libc.stdlib cimport calloc, free
from posix.stdlib cimport srand48
from libc.string cimport strncpy
from pysam import FastxRecord, AlignedSegment, qualitystring_to_array, CMATCH, CINS, CDEL, CSOFT_CLIP
from pybwa.libbwaindex cimport force_bytes
from pybwa.libbwaindex cimport BwaIndex
from pybwa.libbwa cimport bwa_verbose

from pysam.libcalignedsegment cimport makeAlignedSegment


__all__ = [
    "BwaAlnOptions",
    "BwaAln",
]

cpdef bint _set_bwa_aln_verbosity(int level):
    """Set the BWA C-API verbosity, returning True if changed, false otherwise."""
    global bwa_verbose
    retval = level != bwa_verbose
    bwa_verbose = level
    return retval


cdef class BwaAlnOptions:
    """The container for options for :class:`pybwa.BwaAln`.
    
    Args:
        max_mismatches (int | None): :code:`-n <int>`
        max_gap_opens (int | None): :code:`-o <int>`
        max_gap_extensions (int | None): :code:`-e <int>`
        min_indel_to_end_distance (int | None): :code:`-i <int>`
        max_occurrences_for_extending_long_deletion (int | None): :code:`-d <int>`
        seed_length (int | None): :code:`-l <int>`
        max_mismatches_in_seed (int | None): :code:`-k <int>`
        mismatch_penalty (int | None): :code:`-M <int>`
        gap_open_penalty (int | None): :code:`-O <int>`
        gap_extension_penalty (int | None): :code:`-E <int>`
        stop_at_max_best_hits (int | None): :code:`-R <int>`
        max_hits (int | None): :code:`bwa samse -n <int>`
        log_scaled_gap_penalty (bool | None): :code:`-L`
        find_all_hits (bool | None): :code:`-N`
        with_md (bool): output the MD to each alignment in the XA tag, otherwise use :code:`"."`
        threads (int): the number of threads to use
    """

    def __init__(self,
                 max_mismatches: int = -1,
                 max_gap_opens: int = 1,
                 max_gap_extensions: int = 6,
                 min_indel_to_end_distance: int = 5,
                 max_occurrences_for_extending_long_deletion: int = 10,
                 seed_length: int = 32,
                 max_mismatches_in_seed: int = 2,
                 mismatch_penalty: int = 3,
                 gap_open_penalty: int = 11,
                 gap_extension_penalty: int = 4,
                 stop_at_max_best_hits: int = 30,
                 max_hits: int = 3,
                 log_scaled_gap_penalty: bool = False,
                 find_all_hits: bool = False,
                 with_md: bool = False,
                 threads: int = 1
                 ):
        if max_mismatches is not None:
            self.max_mismatches = max_mismatches
        if max_gap_opens is not None:
            self.max_gap_opens = max_gap_opens
        if max_gap_extensions is not None:
            self.max_gap_extensions = max_gap_extensions
            # this mimics the behavior in bwtaln.c, where the default # of gap extensions is 6
            # and BWA_MODE_GAPE is set.  But when a value is given on the command line, then
            # BWA_MODE_GAPE is unset.  So in this constructor if the default value is given, then
            # we keep BWA_MODE_GAPE set, and if this parameter is set via the setter, we unset it.
            if max_gap_extensions == 6:
                self._delegate.mode |= BWA_MODE_GAPE
        if min_indel_to_end_distance is not None:
            self.min_indel_to_end_distance = min_indel_to_end_distance
        if max_occurrences_for_extending_long_deletion is not None:
            self.max_occurrences_for_extending_long_deletion = max_occurrences_for_extending_long_deletion
        if seed_length is not None:
            self.seed_length = seed_length
        if max_mismatches_in_seed is not None:
            self.max_mismatches_in_seed = max_mismatches_in_seed
        if mismatch_penalty is not None:
            self.mismatch_penalty = mismatch_penalty
        if gap_open_penalty is not None:
            self.gap_open_penalty = gap_open_penalty
        if gap_extension_penalty is not None:
            self.gap_extension_penalty = gap_extension_penalty
        if stop_at_max_best_hits is not None:
            self.stop_at_max_best_hits = stop_at_max_best_hits
        if max_hits is not None:
            self.max_hits = max_hits
        if log_scaled_gap_penalty is not None:
            self.log_scaled_gap_penalty = log_scaled_gap_penalty
        if find_all_hits is not None:
            self.find_all_hits = find_all_hits
        if with_md is not None:
            self.with_md = with_md
        if threads is not None:
            self.threads = threads

    def __cinit__(self):
        self._delegate = gap_init_opt()

    def __dealloc__(self):
        free(self._delegate)

    cdef gap_opt_t* gap_opt(self):
        """Returns the options struct to use with the bwa C library methods"""
        return self._delegate

    property max_mismatches:
        """:code:`bwa aln -n <int>`"""
        def __get__(self) -> int:
            return self._delegate.max_diff
        def __set__(self, value: int):
            if value >= 0:
                self._delegate.fnr = -1.0
                self._delegate.max_diff = value
            else:
                self._delegate.max_diff = -1
                self._delegate.fnr = 0.04

    property max_gap_opens:
        """:code:`bwa aln -o <int>`"""
        def __get__(self) -> int:
            return self._delegate.max_gapo
        def __set__(self, value: int):
           self._delegate.max_gapo = value

    property max_gap_extensions:
        """:code:`bwa aln -e <int>`"""
        def __get__(self) -> int:
            return self._delegate.max_gape
        def __set__(self, value: int):
            self._delegate.max_gape = value
            # the BWA_MODE_GAPE mode indicates that gap extensions
            # should count against the maximum # of mismatches
            if self._delegate.max_gape > 0:
                self._delegate.mode &= ~BWA_MODE_GAPE
            else:
                self._delegate.mode |= BWA_MODE_GAPE

    property min_indel_to_end_distance:
        """:code:`bwa aln -i <int>`"""
        def __get__(self) -> int:
            return self._delegate.indel_end_skip
        def __set__(self, value: int):
           self._delegate.indel_end_skip = value

    property max_occurrences_for_extending_long_deletion:
        """:code:`bwa aln -d <int>`"""
        def __get__(self) -> int:
            return self._delegate.max_del_occ
        def __set__(self, value: int):
           self._delegate.max_del_occ = value

    property seed_length:
        """:code:`bwa aln -l <int>`"""
        def __get__(self) -> int:
            return self._delegate.seed_len
        def __set__(self, value: int):
           self._delegate.seed_len = value

    property max_mismatches_in_seed:
        """:code:`bwa aln -k <int>`"""
        def __get__(self) -> int:
            return self._delegate.max_seed_diff
        def __set__(self, value: int):
           self._delegate.max_seed_diff = value

    property mismatch_penalty:
        """:code:`bwa aln -M <int>`"""
        def __get__(self) -> int:
            return self._delegate.s_mm
        def __set__(self, value: int):
           self._delegate.s_mm = value

    property gap_open_penalty:
        """:code:`bwa aln -O <int>`"""
        def __get__(self) -> int:
            return self._delegate.s_gapo
        def __set__(self, value: int):
           self._delegate.s_gapo = value

    property gap_extension_penalty:
        """:code:`bwa aln -E <int>`"""
        def __get__(self) -> int:
            return self._delegate.s_gape
        def __set__(self, value: int):
           self._delegate.s_gape = value

    property stop_at_max_best_hits:
        """:code:`bwa aln -R <int>`"""
        def __get__(self) -> int:
            return self._delegate.max_top2
        def __set__(self, value: int):
           self._delegate.max_top2 = value

    property max_hits:
        """:code:`bwa samse -n <int>`"""
        def __get__(self) -> int:
            return self._max_hits
        def __set__(self, value: int):
           self._max_hits = value

    property log_scaled_gap_penalty:
        """:code:`bwa aln -L`"""
        def __get__(self) -> bool:
            return self._delegate.mode & BWA_MODE_LOGGAP > 0
        def __set__(self, value: bool):
            if value:
                self._delegate.mode |= BWA_MODE_LOGGAP
            else:
                self._delegate.mode &= ~BWA_MODE_LOGGAP

    property find_all_hits:
        """:code:`bwa aln -N`"""
        def __get__(self) -> bool:
            return self._delegate.mode & BWA_MODE_NONSTOP > 0
        def __set__(self, value: bool):
            if value:
                self._delegate.mode |= BWA_MODE_NONSTOP
                self.stop_at_max_best_hits = 0x7fffffff
            else:
                self._delegate.mode &= ~BWA_MODE_NONSTOP

    property with_md:
        """:code:`bwa samse -d`
        
        Output the MD to each alignment in the XA tag, otherwise use :code:`"."`.
        """
        def __get__(self) -> bool:
            return self._with_md
        def __set__(self, value: bool):
           self._with_md = value

    property threads:
        """:code:`bwa aln -t`"""
        def __get__(self) -> int:
            return self._delegate.n_threads
        def __set__(self, value: int):
            self._delegate.n_threads = value


cdef class BwaAln:
    """The class to align reads with :code:`bwa aln`."""

    cdef BwaIndex _index

    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None):
        """Constructs the :code:`bwa aln` aligner.

        One of `prefix` or `index` must be specified.

        Args:
            prefix: the path prefix for the BWA index (typically a FASTA)
            index: the index to use
        """
        if prefix is not None:
            assert Path(prefix).exists()
            self._index = BwaIndex(prefix=prefix)
        elif index is not None:
            self._index = index
        else:
            raise ValueError("Either prefix or index must be given")

        bwase_initialize()

    def align(self, queries: List[FastxRecord] | List[str] | None = None, opt: BwaAlnOptions | None = None) -> List[AlignedSegment]:
        """Align one or more queries with `bwa aln`.

        Args:
            queries: the queries to align
            opt: the alignment options | None = None, or None to use the default options

        Returns:
            one alignment (:class:`~pysam.AlignedSegment`) per query
            :code:`List[List[AlignedSegment]]`.
        """
        if len(queries) == 0:
            return []
        elif isinstance(queries[0], str):
            queries = [
                FastxRecord(name=f"read.{i}", sequence=sequence)
                for i, sequence in enumerate(queries)
            ]
        opt = BwaAlnOptions() if opt is None else opt

        return self._calign(opt,  queries)

    cdef _copy_seq(self, q: FastxRecord, bwa_seq_t* s, is_comp: bool):
        seq_len = len(q.sequence)
        s.len = seq_len
        s.clip_len = seq_len
        s.full_len = seq_len
        s.seq = <uint8_t *> calloc(sizeof(uint8_t), seq_len + 1)
        s.rseq = <uint8_t *> calloc(sizeof(uint8_t), seq_len + 1)

        # convert char into int
        for i, base in enumerate(q.sequence):
            s.seq[i] = nst_nt4_table[ord(base)]
            s.rseq[i] = nst_nt4_table[ord(base)]
        s.seq[seq_len] = b'\0'

        # qualities
        if q.quality is None:
            s.qual = NULL
        else:
            s.qual = <uint8_t *> calloc(sizeof(uint8_t), seq_len + 1)
            qual_str = force_bytes(q.quality)
            for i in range(seq_len):
                s.qual[i] = qual_str[i]
            s.qual[seq_len] = b'\0'

        # use seq_reverse from bwaseqio.c
        seq_reverse(seq_len, s.seq,
                    0)  #  // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
        seq_reverse(seq_len, s.rseq, 1 if is_comp else 0)

        s.name = <char *> calloc(sizeof(char), len(q.name) + 1)
        strncpy(s.name, force_bytes(q.name), len(q.name))
        s.name[len(q.name)] = b'\0'

    cdef _calign(self, opt: BwaAlnOptions, queries: List[FastxRecord]):
        cdef bwa_seq_t* seqs
        cdef bwa_seq_t* s
        cdef char* s_char
        cdef gap_opt_t* gap_opt
        cdef sam_hdr_t *h
        cdef kstring_t hdr_str
        cdef bam1_t **bams

        hdr_str.l = hdr_str.m = 0
        hdr_str.s = NULL
        bwa_format_sam_hdr(self._index.bns(), NULL, &hdr_str)
        h = sam_hdr_parse(hdr_str.l, hdr_str.s)
        h.l_text = hdr_str.l
        h.text = hdr_str.s

        gap_opt = opt.gap_opt()

        # copy FastqProxy into bwa_seq_t
        num_seqs = len(queries)
        seqs = <bwa_seq_t*>calloc(sizeof(bwa_seq_t), num_seqs)
        for i in range(num_seqs):
            self._copy_seq(queries[i], &seqs[i], (opt._delegate.mode & BWA_MODE_COMPREAD) != 0)
            seqs[i].tid = -1

        # this is `bwa aln` and `bwa se` combined to be multi-threaded.
        bams = bwa_aln_and_samse(self._index.bns(), self._index.bwt(),  self._index.pac(), h, num_seqs, seqs, gap_opt, opt.max_hits, opt.with_md)

        # create the AlignedSegment
        recs = []
        for i in range(num_seqs):
            rec = makeAlignedSegment(bams[i], self._index.header)
            bam_destroy1(bams[i])
            recs.append(rec)
        free(bams)

        bwa_free_read_seq(num_seqs, seqs)
        free(hdr_str.s)

        return recs

    cpdef reinitialize_seed(self):
        """Re-initializes the random seed."""
        srand48(self._index.bns().seed)


###########################################################
# Below is code to facilitate testing
###########################################################


cdef _assert_gap_opt_are_the_same_c():
    """Tests that the defaults are synced between bwa and pybwa."""
    cdef gap_opt_t *bwa_gap_opt
    cdef gap_opt_t *pybwa_gap_opt
    bwa_gap_opt = gap_init_opt()
    options = BwaAlnOptions()
    pybwa_gap_opt = options.gap_opt()

    # Check member by member for the gap_opt_t struct
    assert pybwa_gap_opt.s_mm == bwa_gap_opt.s_mm
    assert pybwa_gap_opt.s_gape == bwa_gap_opt.s_gape
    assert pybwa_gap_opt.mode == bwa_gap_opt.mode
    assert pybwa_gap_opt.indel_end_skip == bwa_gap_opt.indel_end_skip
    assert pybwa_gap_opt.max_del_occ == bwa_gap_opt.max_del_occ
    assert pybwa_gap_opt.max_entries == bwa_gap_opt.max_entries
    assert pybwa_gap_opt.fnr == bwa_gap_opt.fnr
    assert pybwa_gap_opt.max_diff == bwa_gap_opt.max_diff
    assert pybwa_gap_opt.max_gapo == bwa_gap_opt.max_gapo
    assert pybwa_gap_opt.max_gape == bwa_gap_opt.max_gape
    assert pybwa_gap_opt.max_seed_diff == bwa_gap_opt.max_seed_diff
    assert pybwa_gap_opt.seed_len == bwa_gap_opt.seed_len
    assert pybwa_gap_opt.n_threads == bwa_gap_opt.n_threads
    assert pybwa_gap_opt.max_top2 == bwa_gap_opt.max_top2
    assert pybwa_gap_opt.trim_qual == bwa_gap_opt.trim_qual

    free(bwa_gap_opt)


def _assert_gap_opt_are_the_same() -> None:
    """Tests that the defaults are synced between bwa and pybwa."""
    _assert_gap_opt_are_the_same_c()
