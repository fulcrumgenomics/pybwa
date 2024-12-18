# cython: language_level=3
from ctypes import memset

from libc.stdlib cimport calloc, free
import enum

cdef class BwaMemOptions:
    """The container for options for [`BwaMem`][bwapy.BwaMem].

    Use [`BwaMemptionsBuilder`][bwapy.BwaMemOptionsBuilder] to use and set custom options.
    """
    cdef mem_opt_t * _delegate
    _ignore_alt: bool

    def __init__(self):
        self._ignore_alt = False
        self._cinit()

    cdef _cinit(self):
        self._delegate = mem_opt_init()

    def __dealloc__(self):
        free(self._delegate)


# class syntax
@enum.unique
class BwaMemMode(enum.Enum):
    """The read type for overriding multiple options"""
    PACBIO = enum.auto()
    ONT2D = enum.auto()
    INTRACTG = enum.auto()


cdef class BwaMemOptionsBuilder:
    """Builder for options for [`BwaAln`][bwapy.BwaAln]"""
    _options: BwaMemOptions
    _options0: BwaMemOptions  # attributes are set to `1` when options are set in _options
    _mode: BwaMemMode | None

    def __init__(self, options: BwaMemOptions | None = None) -> None:
        self._options: BwaMemOptions = BwaMemOptions() if options is None else options
        self._options0: BwaMemOptions = BwaMemOptions()
        self._mode = None
        memset(self._options0._delegate, 0, sizeof(mem_opt_t))

    def build(self) -> BwaMemOptions:
        if self._mode is None:
            # matching score is changed so scale the rest of the penalties/scores
            if self._options0._delegate.a == 1:
                if self._options0._delegate.b != 1:
                    self._options._delegate.b *= self._options._delegate.a
                if self._options0._delegate.T != 1:
                    self._options._delegate.T *= self._options._delegate.a
                if self._options0._delegate.o_del != 1:
                    self._options._delegate.o_del *= self._options._delegate.a
                if self._options0._delegate.e_del != 1:
                    self._options._delegate.e_del *= self._options._delegate.a
                if self._options0._delegate.o_ins != 1:
                    self._options._delegate.o_ins *= self._options._delegate.a
                if self._options0._delegate.e_ins != 1:
                    self._options._delegate.e_ins *= self._options._delegate.a
                if self._options0._delegate.zdrop != 1:
                    self._options._delegate.zdrop *= self._options._delegate.a
                if self._options0._delegate.pen_clip5 != 1:
                    self._options._delegate.pen_clip5 *= self._options._delegate.a
                if self._options0._delegate.pen_clip3 != 1:
                    self._options._delegate.pen_clip3 *= self._options._delegate.a
                if self._options0._delegate.pen_unpaired != 1:
                    self._options._delegate.pen_unpaired *= self._options._delegate.a
        elif self._mode == BwaMemMode.INTRACTG:
            if self._options0._delegate.o_del != 1:
                self._options._delegate.o_del = 16
            if self._options0._delegate.o_ins != 1:
                self._options._delegate.o_ins = 16
            if self._options0._delegate.b != 1:
                self._options._delegate.b = 9
            if self._options0._delegate.pen_clip5 != 1:
                self._options._delegate.pen_clip5 = 5
            if self._options0._delegate.pen_clip3 != 1:
                self._options._delegate.pen_clip3 = 5
        else:
            if self._options0._delegate.o_del != 1:
                self._options._delegate.o_del = 1
            if self._options0._delegate.o_ins != 1:
                self._options._delegate.o_ins = 1
            if self._options0._delegate.e_del != 1:
                self._options._delegate.e_del = 1
            if self._options0._delegate.e_ins != 1:
                self._options._delegate.e_ins = 1
            if self._options0._delegate.b != 1:
                self._options._delegate.b = 1
            if self._options0._delegate.split_factor == 0.0:
                self._options._delegate.split_factor = 10.0
            if self._options0._delegate.pen_clip5 != 1:
                self._options._delegate.pen_clip5 = 0
            if self._options0._delegate.pen_clip3 != 1:
                self._options._delegate.pen_clip3 = 0
            # ONT2D vs PACBIO options
            if self._options0._delegate.min_chain_weight != 1:
                self._options._delegate.min_chain_weight = 20 if self._mode == BwaMemMode.ONT2D else 40
            if self._options0._delegate.min_seed_len != 1:
                self._options._delegate.min_seed_len = 14 if self._mode == BwaMemMode.ONT2D else 17

        bwa_fill_scmat(
            self._options._delegate.a, self._options._delegate.b, self._options._delegate.mat
        )

        return self._options

    cdef min_seed_len(self, value: int):
        """bwa mem -k <int>"""
        self._options._delegate.min_seed_len = value
        self._options0._delegate.min_seed_len = 1
        return self

    # FIXME: convert to an enum enum
    cdef mode(self, value: BwaMemMode):
        """bwa mem -x <str>"""
        self._mode = value
        return self

    cdef band_width(self, value: int):
        """bwa mem -w <int>"""
        self._options._delegate.w = value
        self._options0._delegate.w = 1
        return self

    cdef band_width(self, value: int):
        """bwa mem -w <int>"""
        self._options._delegate.w = value
        self._options0._delegate.w = 1
        return self

    cdef match_score(self, value: int):
        """bwa mem -A <int>"""
        self._options._delegate.a = value
        self._options0._delegate.a = 1
        return self

    cdef mismatch_penalty(self, value: int):
        """bwa mem -A <int>"""
        self._options._delegate.b = value
        self._options0._delegate.b = 1
        return self

    cdef minimum_score(self, value: int):
        """bwa mem -T <int>"""
        self._options._delegate.T = value
        self._options0._delegate.T = 1
        return self

    cdef unpaired_penalty(self, value: int):
        """bwa mem -U <int>"""
        self._options._delegate.pen_unpaired = value
        self._options0._delegate.pen_unpaired = 1
        return self

    cdef unpaired_penalty(self, value: int):
        """bwa mem -U <int>"""
        self._options._delegate.pen_unpaired = value
        self._options0._delegate.pen_unpaired = 1
        return self

    cdef num_threads(self, value: int):
        """bwa mem -t <int>"""
        self._options._delegate.num_threads = value if value > 1 else 1
        return self

    def _set_flag(self, value: bool, flag: int):
        if value:
            self._options._delegate.flag |= flag
        else:
            self._options._delegate.flag |= ~flag
        return self

    cdef skip_pairing(self, value: bool):
        """bwa mem -P"""
        return self._set_flag(value, MEM_F_NOPAIRING)

    cdef output_all_for_fragments(self, value: bool):
        """bwa mem -a"""
        return self._set_flag(value, MEM_F_ALL)

    cdef interleaved_paired_end(self, value: bool):
        """bwa mem -p"""
        return self._set_flag(value, MEM_F_PE | MEM_F_SMARTPE)

    cdef skip_mate_rescue(self, value: bool):
        """bwa mem -S"""
        return self._set_flag(value, MEM_F_NO_MULTI)

    cdef soft_clip_supplementary(self, value: bool):
        """bwa mem -Y"""
        return self._set_flag(value, MEM_F_SOFTCLIP)

    cdef with_xr_tag(self, value: bool):
        """bwa mem -V"""
        return self._set_flag(value, MEM_F_REF_HDR)

    cdef query_coord_as_primary(self, value: bool):
        """bwa mem -5"""
        return self._set_flag(value, MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ) #/ always apply MEM_F_KEEP_SUPP_MAPQ with -5

    cdef keep_mapq_for_supplementary(self, value: bool):
        """bwa mem -q"""
        return self._set_flag(value, MEM_F_KEEP_SUPP_MAPQ)

    cdef with_xb_tag(self, value: bool):
        """bwa mem -u"""
        return self._set_flag(value, MEM_F_XB)

    cdef max_occurrences(self, value: int):
        """bwa mem -c <int>"""
        self._options._delegate.max_occ = value
        self._options0._delegate.max_occ = 1
        return self

    cdef off_diagonal_x_dropoff(self, value: int):
        """bwa mem -d <float>"""
        self._options._delegate.XA_drop_ratio = value
        self._options0._delegate.XA_drop_ratio = 1
        return self

    cdef ignore_alternate_contigs(self, value: bool):
        """bwa mem -j"""
        self._options._delegate._ignore_alt = value
        return self

    cdef internal_seed_split_factor(self, value: float):
        """bwa mem -r <float>"""
        self._options._delegate.split_factor = value
        self._options0._delegate.split_factor = 1
        return self

    cdef drop_chain_fraction(self, value: float):
        """bwa mem -D <float>"""
        self._options._delegate.drop_ratio = value
        self._options0._delegate.drop_ratio = 1
        return self

    cdef max_mate_rescue_rounds(self, value: int):
        """bwa mem -m <float>"""
        self._options._delegate.max_matesw = value
        self._options0._delegate.max_matesw = 1
        return self

    cdef min_seeded_bases_in_chain(self, value: int):
        """bwa mem -W <float>"""
        self._options._delegate.min_chain_weight = value
        self._options0._delegate.min_chain_weight = 1
        return self

    cdef seed_occurrence_in_3rd_round(self, value: int):
        """bwa mem -y <float>"""
        self._options._delegate.max_mem_intv = value
        self._options0._delegate.max_mem_intv = 1
        return self

    cdef xa_max_hits(self, value: int, alt_value: int | None):
        """bwa mem -h <int<,int>>"""
        self._options0._delegate.max_XA_hits = 1
        self._options0._delegate.max_XA_hits_alt = 1
        self._options._delegate.max_XA_hits = value
        self._options._delegate.max_XA_hits_alt = value if alt_value is None else alt_value

    cdef xa_drop_ration(self, value: float):
        """bwa mem -y <float>"""
        self._options._delegate.XA_drop_ratio = value
        return self

    cdef gap_open_penalty(self, deletions: int, insertions: int):
        """bwa mem -O <int<,int>>"""
        self._options0._delegate.o_del = 1
        self._options0._delegate.o_ins = 1
        self._options._delegate.o_del = deletions
        self._options._delegate.o_ins = insertions
        return self

    cdef gap_extension_penalty(self, deletions: int, insertions: int):
        """bwa mem -E <int<,int>>"""
        self._options0._delegate.e_del = 1
        self._options0._delegate.e_ins = 1
        self._options._delegate.e_del = deletions
        self._options._delegate.e_ins = insertions
        return self

    cdef clipping_penalty(self, five_prime: int, three_prime: int):
        """bwa mem -L <int<,int>>"""
        self._options0._delegate.pen_clip5 = 1
        self._options0._delegate.pen_clip3 = 1
        self._options._delegate.pen_clip5 = five_prime
        self._options._delegate.pen_clip3 = three_prime


    cdef insert_size_distribution(self,
                                  mean: float,
                                  stddev: float | None = None,
                                  max: float | None = None,
                                  min: float | None = None):
        """bwa mem -I float<,float<,intL<,int>>>"""
        raise NotImplementedError


#