# cython: language_level=3

from libc.stdlib cimport free


__all__ = [
    "BwaAlnOptions",
]


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

    _max_hits: int
    _with_md: bool

    def __init__(self,
                 max_mismatches: int | None = None,
                 max_gap_opens: int | None = None,
                 max_gap_extensions: int | None = None,
                 min_indel_to_end_distance: int | None = None,
                 max_occurrences_for_extending_long_deletion: int | None = None,
                 seed_length: int | None = None,
                 max_mismatches_in_seed: int | None = None,
                 mismatch_penalty: int | None = None,
                 gap_open_penalty: int | None = None,
                 gap_extension_penalty: int | None = None,
                 stop_at_max_best_hits: int | None = None,
                 max_hits: int | None = 3,
                 log_scaled_gap_penalty: bool | None = None,
                 find_all_hits: bool | None = None,
                 with_md: bool | None = False,
                 threads: int | None = None
                 ):
        if max_mismatches is not None:
            self.max_mismatches = max_mismatches
        if max_gap_opens is not None:
            self.max_gap_opens = max_gap_opens
        if max_gap_extensions is not None:
            self.max_gap_extensions = max_gap_extensions
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
            self._delegate.fnr = -1.0
            self._delegate.max_diff = value

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
        """:code:`bwa samse -d
        
        Output the MD to each alignment in the XA tag, otherwise use :code:`"."`.
        """
        def __get__(self) -> bool:
            return self._with_md
        def __set__(self, value: bool):
           self._with_md = value

    property threads:
        """:code:`bwa aln -t"""
        def __get__(self) -> int:
            return self._delegate.n_threads
        def __set__(self, value: int):
            self._delegate.n_threads = value

