# cython: language_level=3
from pathlib import Path
from typing import List

from fgpyo.sequence import reverse_complement
from libc.string cimport memcpy
from libc.stdlib cimport calloc, free
import enum
from pybwa.libbwaindex cimport BwaIndex
from pysam import FastxRecord, AlignedSegment, qualitystring_to_array
from libc.string cimport strncpy
from pybwa.libbwaindex cimport force_bytes


__all__ = [
    "BwaMemMode",
    "BwaMemOptions",
    "BwaMem",
]


@enum.unique
class BwaMemMode(enum.Enum):
    """The read type for overriding multiple options"""
    PACBIO = enum.auto()
    """PacBio reads to ref"""
    ONT2D = enum.auto()
    """Oxford Nanopore 2D-reads to ref"""
    INTRACTG = enum.auto()
    """intra-species contigs to ref"""


cdef class BwaMemOptions:
    """The container for options for :class:`~pybwa.BwaMem`.

    Args:
        min_seed_len (int): :code:`bwa mem -k <int>`
        mode (BwaMemMode): :code:`bwa mem -x <str>`
        band_width (int): :code:`bwa mem -w <int>`
        match_score (int): :code:`bwa mem -A <int>`
        mismatch_penalty (int): :code:`bwa mem -B <int>`
        minimum_score (int): :code:`bwa mem -T <int>`
        unpaired_penalty (int): :code:`bwa mem -U <int>`
        n_threads (int): :code:`bwa mem -t <int>`
        skip_pairing (bool): :code:`bwa mem -P`
        output_all_for_fragments (bool): :code:`bwa mem -a`
        interleaved_paired_end (bool): :code:`bwa mem -p`
        short_split_as_secondary (bool): :code:`bwa mem -M`
        skip_mate_rescue (bool): :code:`bwa mem -S`
        soft_clip_supplementary (bool): :code:`bwa mem -Y`
        with_xr_tag (bool): :code:`bwa mem -V`
        query_coord_as_primary (bool): :code:`bwa mem -5`
        keep_mapq_for_supplementary (bool): :code:`bwa mem -q`
        with_xb_tag (bool): :code:`bwa mem -u`
        max_occurrences (int): :code:`bwa mem -c <int>`
        off_diagonal_x_dropoff (int): :code:`bwa mem -d <float>`
        ignore_alternate_contigs (bool): :code:`bwa mem -j`
        internal_seed_split_factor (float): :code:`bwa mem -r <float>`
        drop_chain_fraction (float): :code:`bwa mem -D <float>`
        max_mate_rescue_rounds (int): :code:`bwa mem -m <int>`
        min_seeded_bases_in_chain (int): :code:`bwa mem -W <int>`
        seed_occurrence_in_3rd_round (int): :code:`bwa mem -y <int>`
        xa_max_hits (int | tuple[int, int]): :code:`bwa mem -h <int<,int>>`
        xa_drop_ratio (float): :code:`bwa mem -z <float>`
        gap_open_penalty (int | tuple[int, int]): :code:`bwa mem -O <int<,int>>`
        gap_extension_penalty (int | tuple[int, int]): :code:`bwa mem -E <int<,int>>`
        clipping_penalty (int | tuple[int, int]): :code:`bwa mem -L <int<,int>>`
    """
    _finalized: bool
    _ignore_alt: bool
    _mode: BwaMemMode | None
    cdef mem_opt_t* _options
    cdef mem_opt_t* _options0

    def _assert_not_finalized(self, attr_name: str) -> None:
        """Raises an AttributeError if the options have been finalized.

        This is used in each setter below to enforce that setters are not used after
         :meth:`~pybwa.BwaMemOptions.finalize` is called.  This could be a property decorator,
         but cython does not support decorating properties.
        """
        if self._finalized:
            raise AttributeError(f"can't set attribute: {attr_name}")

    def __init__(self,
                 min_seed_len: int | None = None,
                 mode: BwaMemMode | None = None,
                 band_width: int | None = None,
                 match_score: int | None = None,
                 mismatch_penalty: int | None = None,
                 minimum_score: int | None = None,
                 unpaired_penalty: int | None = None,
                 n_threads: int | None = None,
                 skip_pairing: bool | None = None,
                 output_all_for_fragments: bool | None = None,
                 interleaved_paired_end: bool | None = None,
                 short_split_as_secondary: bool | None = None,
                 skip_mate_rescue: bool | None = None,
                 soft_clip_supplementary: bool | None = None,
                 with_xr_tag: bool | None = None,
                 query_coord_as_primary: bool | None = None,
                 keep_mapq_for_supplementary: bool | None = None,
                 with_xb_tag: bool | None = None,
                 max_occurrences: int | None = None,
                 off_diagonal_x_dropoff: int | None = None,
                 ignore_alternate_contigs: bool | None = None,
                 internal_seed_split_factor: float | None = None,
                 drop_chain_fraction: float | None = None,
                 max_mate_rescue_rounds: int | None = None,
                 min_seeded_bases_in_chain: int | None = None,
                 seed_occurrence_in_3rd_round: int | None = None,
                 xa_max_hits: int | tuple[int, int] | None = None,
                 xa_drop_ratio: float | None = None,
                 gap_open_penalty: int | tuple[int, int] | None = None,
                 gap_extension_penalty: int | tuple[int, int] | None = None,
                 clipping_penalty: int | tuple[int, int] | None = None,
                 threads: int | None = None,
                 chunk_size: int | None = None) -> None:
        self._finalized = False
        self._ignore_alt = False
        self._mode = None

        # for options that scale by match score or are set by mode, zero out the `options0`, which
        # tracks if they have been set or not.
        self._options0.b = 0
        self._options0.T = 0
        self._options0.o_del = 0
        self._options0.e_del = 0
        self._options0.o_ins = 0
        self._options0.e_ins = 0
        self._options0.zdrop = 0
        self._options0.pen_clip5 = 0
        self._options0.pen_clip3 = 0
        self._options0.pen_unpaired = 0
        self._options0.split_factor = 0
        self._options0.min_chain_weight = 0
        self._options0.min_seed_len = 0

        if min_seed_len is not None:
            self.min_seed_len = min_seed_len
        if mode is not None:
            self.mode = mode
        if band_width is not None:
            self.band_width = band_width
        if match_score is not None:
            self.match_score = match_score
        if mismatch_penalty is not None:
            self.mismatch_penalty = mismatch_penalty
        if minimum_score is not None:
            self.minimum_score = minimum_score
        if unpaired_penalty is not None:
            self.unpaired_penalty = unpaired_penalty
        if n_threads is not None:
            self.n_threads = n_threads
        if skip_pairing is not None:
            self.skip_pairing = skip_pairing
        if output_all_for_fragments is not None:
            self.output_all_for_fragments = output_all_for_fragments
        if interleaved_paired_end is not None:
            self.interleaved_paired_end = interleaved_paired_end
        if short_split_as_secondary is not None:
            self.short_split_as_secondary = short_split_as_secondary
        if skip_mate_rescue is not None:
            self.skip_mate_rescue = skip_mate_rescue
        if soft_clip_supplementary is not None:
            self.soft_clip_supplementary = soft_clip_supplementary
        if with_xr_tag is not None:
            self.with_xr_tag = with_xr_tag
        if query_coord_as_primary is not None:
            self.query_coord_as_primary = query_coord_as_primary
        if keep_mapq_for_supplementary is not None:
            self.keep_mapq_for_supplementary = keep_mapq_for_supplementary
        if with_xb_tag is not None:
            self.with_xb_tag = with_xb_tag
        if max_occurrences is not None:
            self.max_occurrences = max_occurrences
        if off_diagonal_x_dropoff is not None:
            self.off_diagonal_x_dropoff = off_diagonal_x_dropoff
        if ignore_alternate_contigs is not None:
            self.ignore_alternate_contigs = ignore_alternate_contigs
        if internal_seed_split_factor is not None:
            self.internal_seed_split_factor = internal_seed_split_factor
        if drop_chain_fraction is not None:
            self.drop_chain_fraction = drop_chain_fraction
        if max_mate_rescue_rounds is not None:
            self.max_mate_rescue_rounds = max_mate_rescue_rounds
        if min_seeded_bases_in_chain is not None:
            self.min_seeded_bases_in_chain = min_seeded_bases_in_chain
        if seed_occurrence_in_3rd_round is not None:
            self.seed_occurrence_in_3rd_round = seed_occurrence_in_3rd_round
        if xa_max_hits is not None:
            self.xa_max_hits = xa_max_hits
        if xa_drop_ratio is not None:
            self.xa_drop_ratio = xa_drop_ratio
        if gap_open_penalty is not None:
            self.gap_open_penalty = gap_open_penalty
        if gap_extension_penalty is not None:
            self.gap_extension_penalty = gap_extension_penalty
        if clipping_penalty is not None:
            self.clipping_penalty = clipping_penalty
        if threads is not None:
            self.threads = threads
        if chunk_size is not None:
            self.chunk_size = chunk_size

    def __cinit__(self):
        self._options = mem_opt_init()
        self._options0 = mem_opt_init()

    def __dealloc__(self):
        free(self._options)
        self._options = NULL
        free(self._options0)
        self._options0 = NULL

    cdef mem_opt_t* mem_opt(self):
        """Returns the options struct to use with the bwa C library methods"""
        if not self._finalized:
            raise Exception("Cannot call `mem_opt` until `finalize()` is called")
        return self._options

    @property
    def finalized(self) -> bool:
        """True if the options have been finalized with :meth:`~pybwa.BwaMemOptions.finalize`."""
        return self._finalized

    def finalize(self, copy: bool = False) -> BwaMemOptions:
        """Performs final initialization of these options.  The object returned may not be
        modified further.

        If the mode is given, then the presets are applied to options that have not been explicitly
        set.  Otherwise, if the match score has been set, the match score scales the (-TdBOELU)
        options if they have not been explicitly set.

        Args:
            copy: true to return a finalized copy of this object, false to finalize this object
        """
        opt: BwaMemOptions
        if copy:
            opt = BwaMemOptions()
            memcpy(opt._options, self._options, sizeof(mem_opt_t))
        else:
            opt = self

        if opt._mode is None:
            # matching score is changed so scale the rest of the penalties/scores
            if opt._options0.a == 1:
                if opt._options0.b != 1:
                    opt._options.b *= opt._options.a
                if opt._options0.T != 1:
                    opt._options.T *= opt._options.a
                if opt._options0.o_del != 1:
                    opt._options.o_del *= opt._options.a
                if opt._options0.e_del != 1:
                    opt._options.e_del *= opt._options.a
                if opt._options0.o_ins != 1:
                    opt._options.o_ins *= opt._options.a
                if opt._options0.e_ins != 1:
                    opt._options.e_ins *= opt._options.a
                if opt._options0.zdrop != 1:
                    opt._options.zdrop *= opt._options.a
                if opt._options0.pen_clip5 != 1:
                    opt._options.pen_clip5 *= opt._options.a
                if opt._options0.pen_clip3 != 1:
                    opt._options.pen_clip3 *= opt._options.a
                if opt._options0.pen_unpaired != 1:
                    opt._options.pen_unpaired *= opt._options.a
        elif opt._mode == BwaMemMode.INTRACTG:
            if opt._options0.o_del != 1:
                opt._options.o_del = 16
            if opt._options0.o_ins != 1:
                opt._options.o_ins = 16
            if opt._options0.b != 1:
                opt._options.b = 9
            if opt._options0.pen_clip5 != 1:
                opt._options.pen_clip5 = 5
            if opt._options0.pen_clip3 != 1:
                opt._options.pen_clip3 = 5
        else:
            if opt._options0.o_del != 1:
                opt._options.o_del = 1
            if opt._options0.o_ins != 1:
                opt._options.o_ins = 1
            if opt._options0.e_del != 1:
                opt._options.e_del = 1
            if opt._options0.e_ins != 1:
                opt._options.e_ins = 1
            if opt._options0.b != 1:
                opt._options.b = 1
            if opt._options0.split_factor == 0.0:
                opt._options.split_factor = 10.0
            if opt._options0.pen_clip5 != 1:
                opt._options.pen_clip5 = 0
            if opt._options0.pen_clip3 != 1:
                opt._options.pen_clip3 = 0
            # ONT2D vs PACBIO options
            if opt._options0.min_chain_weight != 1:
                opt._options.min_chain_weight = 20 if opt._mode == BwaMemMode.ONT2D else 40
            if opt._options0.min_seed_len != 1:
                opt._options.min_seed_len = 14 if opt._mode == BwaMemMode.ONT2D else 17

        # the **actual** chunk size by default scales by the # of threads, otherwise use a fixed
        # number
        if opt._options0.chunk_size != 1 and opt._options.chunk_size > 0:
            opt._options.chunk_size *= opt._options.n_threads

        bwa_fill_scmat(
            opt._options.a, opt._options.b, opt._options.mat
        )

        opt._finalized = True
        return opt

    @property
    def min_seed_len(self) -> int:
        """:code:`bwa mem -k <int>`"""
        return self._options.min_seed_len

    @min_seed_len.setter
    def min_seed_len(self, value: int):
        self._assert_not_finalized(attr_name=BwaMemOptions.min_seed_len.__name__)
        self._options.min_seed_len = value
        self._options0.min_seed_len = 1

    @property
    def mode(self) -> BwaMemMode | None:
        """:code:`bwa mem -x <str>`"""
        return self._mode

    @mode.setter
    def mode(self, value: BwaMemMode) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.mode.__name__)
        self._mode = value

    @property
    def band_width(self) -> int:
        """:code:`bwa mem -w <int>`"""
        return self._options.w

    @band_width.setter
    def band_width(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.band_width.__name__)
        self._options.w = value
        self._options0.w = 1

    @property
    def match_score(self) -> int:
        """:code:`bwa mem -A <int>`"""
        return self._options.a


    @match_score.setter
    def match_score(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.match_score.__name__)
        self._options.a = value
        self._options0.a = 1

    @property
    def mismatch_penalty(self) -> int:
        """:code:`bwa mem -B <int>`"""
        return self._options.b

    @mismatch_penalty.setter
    def mismatch_penalty(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.mismatch_penalty.__name__)
        self._options.b = value
        self._options0.b = 1

    @property
    def minimum_score(self) -> int:
        """:code:`bwa mem -T <int>`"""
        return self._options.T

    @minimum_score.setter
    def minimum_score(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.minimum_score.__name__)
        self._options.T = value
        self._options0.T = 1

    @property
    def unpaired_penalty(self) -> int:
        """:code:`bwa mem -U <int>`"""
        return self._options.pen_unpaired

    @unpaired_penalty.setter
    def unpaired_penalty(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.unpaired_penalty.__name__)
        self._options.pen_unpaired = value
        self._options0.pen_unpaired = 1

    @property
    def n_threads(self) -> int:
        """:code:`bwa mem -t <int>`"""
        return self._options.n_threads

    @n_threads.setter
    def n_threads(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.n_threads.__name__)
        self._options.n_threads = value if value > 1 else 1

    def _set_flag(self, value: bool, flag: int):
        if value:
            self._options.flag |= flag
        else:
            self._options.flag &= ~flag
        return self

    @property
    def skip_pairing(self) -> bool:
        """:code:`bwa mem -P`"""
        return (self._options.flag & MEM_F_NOPAIRING) != 0

    @skip_pairing.setter
    def skip_pairing(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.skip_pairing.__name__)
        self._set_flag(value, MEM_F_NOPAIRING)

    @property
    def output_all_for_fragments(self) -> bool:
        """:code:`bwa mem -a`"""
        return (self._options.flag & MEM_F_ALL) != 0

    @output_all_for_fragments.setter
    def output_all_for_fragments(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.output_all_for_fragments.__name__)
        self._set_flag(value, MEM_F_ALL)

    @property
    def interleaved_paired_end(self) -> bool:
        """:code:`bwa mem -p`"""
        return (self._options.flag & (MEM_F_PE | MEM_F_SMARTPE)) != 0

    @interleaved_paired_end.setter
    def interleaved_paired_end(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.interleaved_paired_end.__name__)
        self._set_flag(value, MEM_F_PE | MEM_F_SMARTPE)
    
    @property
    def short_split_as_secondary(self) -> bool:
        """:code:`bwa mem -M`"""
        return (self._options.flag & MEM_F_NO_MULTI) != 0

    @short_split_as_secondary.setter
    def short_split_as_secondary(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.short_split_as_secondary.__name__)
        self._set_flag(value, MEM_F_NO_MULTI)

    @property
    def skip_mate_rescue(self) -> bool:
        """:code:`bwa mem -S`"""
        return (self._options.flag & MEM_F_NO_RESCUE) != 0

    @skip_mate_rescue.setter
    def skip_mate_rescue(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.skip_mate_rescue.__name__)
        self._set_flag(value, MEM_F_NO_RESCUE)

    @property
    def soft_clip_supplementary(self) -> bool:
        """:code:`bwa mem -Y`"""
        return (self._options.flag & MEM_F_SOFTCLIP) != 0

    @soft_clip_supplementary.setter
    def soft_clip_supplementary(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.soft_clip_supplementary.__name__)
        self._set_flag(value, MEM_F_SOFTCLIP)

    @property
    def with_xr_tag(self) -> bool:
        """:code:`bwa mem -V`"""
        return (self._options.flag & MEM_F_REF_HDR) != 0

    @with_xr_tag.setter
    def with_xr_tag(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.with_xr_tag.__name__)
        self._set_flag(value, MEM_F_REF_HDR)

    @property
    def query_coord_as_primary(self) -> bool:
        """:code:`bwa mem -5`"""
        return (self._options.flag & (MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ)) != 0

    @query_coord_as_primary.setter
    def query_coord_as_primary(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.query_coord_as_primary.__name__)
        self._set_flag(value, MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ) # always apply MEM_F_KEEP_SUPP_MAPQ with -5

    @property
    def keep_mapq_for_supplementary(self) -> bool:
        """:code:`bwa mem -q`"""
        return (self._options.flag & MEM_F_KEEP_SUPP_MAPQ) != 0

    @keep_mapq_for_supplementary.setter
    def keep_mapq_for_supplementary(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.keep_mapq_for_supplementary.__name__)
        self._set_flag(value, MEM_F_KEEP_SUPP_MAPQ)

    @property
    def with_xb_tag(self) -> bool:
        """:code:`bwa mem -u`"""
        return (self._options.flag & MEM_F_XB) != 0

    @with_xb_tag.setter
    def with_xb_tag(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.with_xb_tag.__name__)
        self._set_flag(value, MEM_F_XB)

    @property
    def max_occurrences(self) -> int:
        """:code:`bwa mem -c <int>`"""
        return self._options.max_occ

    @max_occurrences.setter
    def max_occurrences(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.max_occurrences.__name__)
        self._options.max_occ = value
        self._options0.max_occ = 1

    @property
    def off_diagonal_x_dropoff(self) -> int:
        """:code:`bwa mem -d <float>`"""
        return self._options.zdrop

    @off_diagonal_x_dropoff.setter
    def off_diagonal_x_dropoff(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.off_diagonal_x_dropoff.__name__)
        self._options.zdrop = value
        self._options0.zdrop = 1

    @property
    def ignore_alternate_contigs(self) -> bool:
        """:code:`bwa mem -j`"""
        return self._ignore_alt

    @ignore_alternate_contigs.setter
    def ignore_alternate_contigs(self, value: bool) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.ignore_alternate_contigs.__name__)
        self._ignore_alt = value

    @property
    def internal_seed_split_factor(self) -> float:
        """:code:`bwa mem -r <float>`"""
        return self._options.split_factor

    @internal_seed_split_factor.setter
    def internal_seed_split_factor(self, value: float) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.internal_seed_split_factor.__name__)
        self._options.split_factor = value
        self._options0.split_factor = 1

    @property
    def drop_chain_fraction(self) -> float:
        """:code:`bwa mem -D <float>`"""
        return self._options.drop_ratio

    @drop_chain_fraction.setter
    def drop_chain_fraction(self, value: float) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.drop_chain_fraction.__name__)
        self._options.drop_ratio = value
        self._options0.drop_ratio = 1

    @property
    def max_mate_rescue_rounds(self) -> int:
        """:code:`bwa mem -m <int>`"""
        return self._options.max_matesw

    @max_mate_rescue_rounds.setter
    def max_mate_rescue_rounds(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.max_mate_rescue_rounds.__name__)
        self._options.max_matesw = value
        self._options0.max_matesw = 1

    @property
    def min_seeded_bases_in_chain(self) -> int:
        """:code:`bwa mem -W <int>`"""
        return self._options.min_chain_weight

    @min_seeded_bases_in_chain.setter
    def min_seeded_bases_in_chain(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.min_seeded_bases_in_chain.__name__)
        self._options.min_chain_weight = value
        self._options0.min_chain_weight = 1

    @property
    def seed_occurrence_in_3rd_round(self) -> int:
        """:code:`bwa mem -y <int>`"""
        return self._options.max_mem_intv

    @seed_occurrence_in_3rd_round.setter
    def seed_occurrence_in_3rd_round(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.seed_occurrence_in_3rd_round.__name__)
        self._options.max_mem_intv = value
        self._options0.max_mem_intv = 1

    @property
    def xa_max_hits(self) -> int | tuple[int, int]:
        """:code:`bwa mem -h <int<,int>>`"""
        if self._options.max_XA_hits == self._options.max_XA_hits_alt:
            return self._options.max_XA_hits
        else:
            return self._options.max_XA_hits, self._options.max_XA_hits_alt

    @xa_max_hits.setter
    def xa_max_hits(self, value: int | tuple[int, int]) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.xa_max_hits.__name__)
        self._options0.max_XA_hits = 1
        self._options0.max_XA_hits_alt = 1
        if isinstance(value, int):
            self._options.max_XA_hits = value
            self._options.max_XA_hits_alt = value
        else:
            left, right = value
            self._options.max_XA_hits = left
            self._options.max_XA_hits_alt = right

    @property
    def xa_drop_ratio(self) -> float:
        """:code:`bwa mem -z <float>`"""
        return self._options.XA_drop_ratio

    @xa_drop_ratio.setter
    def xa_drop_ratio(self, value: float) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.xa_drop_ratio.__name__)
        self._options.XA_drop_ratio = value

    @property
    def gap_open_penalty(self) -> int | tuple[int, int]:
        """:code:`bwa mem -O <int<,int>>`"""
        if self._options.o_del == self._options.o_ins:
            return self._options.o_del
        else:
            return self._options.o_del, self._options.o_ins

    @gap_open_penalty.setter
    def gap_open_penalty(self, value: int | tuple[int, int]) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.gap_open_penalty.__name__)
        self._options0.o_del = 1
        self._options0.o_ins = 1
        if isinstance(value, int):
            self._options.o_del = value
            self._options.o_ins = value
        else:
            deletions, insertions = value
            self._options.o_del = deletions
            self._options.o_ins = insertions

    @property
    def gap_extension_penalty(self) -> int | tuple[int, int]:
        """:code:`bwa mem -E <int<,int>>`"""
        if self._options.e_del == self._options.e_ins:
            return self._options.e_del
        else:
            return self._options.e_del, self._options.e_ins

    @gap_extension_penalty.setter
    def gap_extension_penalty(self, value: int | tuple[int, int]) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.gap_extension_penalty.__name__)
        self._options0.e_del = 1
        self._options0.e_ins = 1
        if isinstance(value, int):
            self._options.e_del = value
            self._options.e_ins = value
        else:
            deletions, insertions = value
            self._options.e_del = deletions
            self._options.e_ins = insertions

    @property
    def clipping_penalty(self) -> int | tuple[int, int]:
        """:code:`bwa mem -L <int<,int>>`"""
        if self._options.pen_clip5 == self._options.pen_clip3:
            return self._options.pen_clip5
        else:
            return self._options.pen_clip5, self._options.pen_clip3

    @clipping_penalty.setter
    def clipping_penalty(self, value: int | tuple[int, int]) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.clipping_penalty.__name__)
        self._options0.pen_clip5 = 1
        self._options0.pen_clip3 = 1
        if isinstance(value, int):
            self._options.pen_clip5 = value
            self._options.pen_clip3 = value
        else:
            five_prime, three_prime = value
            self._options.pen_clip5 = five_prime
            self._options.pen_clip3 = three_prime

    @property
    def threads(self) -> int:
        """:code:`bwa mem -t <int>`"""
        return self._options.n_threads

    @threads.setter
    def threads(self, value: int) -> None:
        self._options.n_threads = value

    @property
    def chunk_size(self) -> int:
        """:code:`bwa mem -K <int>`"""
        return self._options.chunk_size

    @chunk_size.setter
    def chunk_size(self, value: int) -> None:
        self._options0.chunk_size = 1
        self._options.chunk_size = value


cdef class BwaMem:
    """The class to align reads with :code:`bwa mem`."""

    cdef BwaIndex _index

    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None):
        """Constructs the :code:`bwa mem` aligner.

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

    # TODO: support paired end
    def align(self, queries: List[FastxRecord] | List[str], opt: BwaMemOptions | None = None) -> List[List[AlignedSegment]]:
        """Align one or more queries with `bwa aln`.

        Args:
            queries: the queries to align
            opt: the alignment options, or None to use the default options

        Returns:
            a list of alignments (:class:`~pysam.AlignedSegment`) per query
            :code:`List[List[AlignedSegment]]`.
        """
        if opt is None:
            opt = BwaMemOptions().finalize()
        elif not opt.finalized:
            opt = opt.finalize(copy=True)

        if len(queries) == 0:
            return []
        elif isinstance(queries[0], str):
            queries = [
                FastxRecord(name=f"read.{i}", sequence=sequence)
                for i, sequence in enumerate(queries)
            ]

        # This mimics how the `bwa mem -K` option works, where we process reads in chunks based on
        # the total number of bases in the reads in the chunk, and making sure we have an _even_
        # number of reads in the chunk.
        start = 0
        results: List[List[AlignedSegment]] = []
        while start < len(queries):
            num_bases = 0
            end = start
            while end < len(queries) and (num_bases < opt.chunk_size or (end-start)&1 == 1):
                num_bases += len(queries[end].sequence)
                end += 1
            assert start < end
            results.extend(self._calign(opt, queries[start:end], n_processed=start))
            start = end

        return results
    
    @staticmethod
    def __to_str(_bytes: bytes) -> str:
        return _bytes.decode('utf-8')

    cdef _copy_seq(self, q: FastxRecord, bseq1_t *s):

        # name
        s.name = <char *> calloc(sizeof(char), len(q.name) + 1)
        strncpy(s.name, force_bytes(q.name), len(q.name))
        s.name[len(q.name)] = b'\0'

        # comment
        # NB: bwa mem supports appending the comment to the SAM tags verbatim! We do not.
        s.comment = NULL

        # sequence
        s.l_seq = len(q.sequence)
        s.seq = <char *> calloc(sizeof(char), s.l_seq + 1)
        for i, base in enumerate(q.sequence):
            s.seq[i] = nst_nt4_table[ord(base)]
        s.seq[s.l_seq] = b'\0'

        # qualities (always ignore)
        s.qual = NULL

    def _unmapped(self, query: FastxRecord) -> AlignedSegment:
        # make a default, unmapped, empty record
        rec = AlignedSegment(header=self._index.header)
        rec.query_name = query.name
        rec.reference_id = -1
        rec.reference_start = -1
        rec.mapping_quality = 0
        rec.is_paired = False
        rec.is_read1 = False
        rec.is_read2 = False
        rec.is_qcfail = False
        rec.is_duplicate = False
        rec.is_secondary = False
        rec.is_supplementary = False
        rec.is_unmapped = True
        rec.query_sequence = query.sequence
        rec.query_qualities = qualitystring_to_array(query.quality)
        return rec

    def _add_sa_tag(self, records: list[AlignedSegment]) -> None:
        num_non_secondary = sum(1 for record in records if not record.is_secondary)
        if num_non_secondary <= 1:
            return
        for i, record in enumerate(records):
            if record.is_secondary:
                continue
            # TODO: bwa mem uses the original cigar, not the one modified 
            # when not using opt.soft_clip_supplementary.  Therefore, we
            # change any hard-clip back to soft-clip
            SA = ""
            for j, other in enumerate(records):
                if i == j or other.is_secondary:
                    continue
                SA += f"{other.reference_name},{other.reference_start+1},"
                SA += "+" if other.is_forward else "-"
                SA += f",{other.cigarstring.replace('H', 'S')}"
                SA += f",{other.mapq},{other.get_tag('NM')};"
            record.set_tag("SA", SA)

    cdef _calign(self, opt: BwaMemOptions, queries: List[FastxRecord], n_processed: int = 0):
        # TODO: ignore_alt
        # TODO: refactor to make this more readable
        cdef bseq1_t* seqs
        cdef bseq1_t* seq
        cdef char* s_char
        cdef kstring_t* kstr
        cdef int take_all
        cdef size_t j
        cdef char **XA
        cdef mem_alnreg_t *mem_alnreg
        cdef mem_aln_t mem_aln
        cdef char *md
        cdef mem_opt_t *mem_opt
        cdef mem_alns_t *mem_alns_vec
        cdef mem_alns_t *mem_alns
        cdef bntann1_t *anno

        recs_to_return: List[List[AlignedSegment]] = []

        # copy FastxRecord into bwa_seq_t
        num_seqs = len(queries)
        mem_opt = opt.mem_opt()
        seqs = <bseq1_t*>calloc(sizeof(bseq1_t), num_seqs)
        for i in range(num_seqs):
            self._copy_seq(queries[i], &seqs[i])

        # process the sequences (ignores the paired end stats)
        mem_alns_vec = mem_process_seqs_alt(mem_opt, self._index.bwt(), self._index.bns(), self._index.pac(),
                         n_processed, num_seqs, seqs, NULL)

        for i in range(num_seqs):
            seq = &seqs[i]
            query = queries[i]
            mem_alns = &mem_alns_vec[i]

            mapped_recs = []
            for j in range(mem_alns.n):
                rec = self._unmapped(query=query)

                mem_aln = mem_alns.a[j]

                # set the flags
                mem_aln.flag |= 0x4 if mem_aln.rid < 0 else 0
                mem_aln.flag |= 0x10 if mem_aln.is_rev > 0 else 0
                rec.flag = (mem_aln.flag & 0xffff) | (0x100 if (mem_aln.flag & 0x10000) != 0 else 0)
                if rec.is_unmapped:
                    continue

                # for secondary alignments, do not add sequence and qualities
                if rec.is_secondary:
                    rec.query_sequence = None
                    rec.query_qualities = None
                elif rec.is_reverse:
                    rec.query_sequence = reverse_complement(query.sequence)
                    if query.quality is not None:
                        # NB: cannot use rec.query_qualities as it is invalidated by setting
                        # query_sequence above
                        rec.query_qualities = qualitystring_to_array(query.quality[::-1])

                # reference id, position, mapq, and cigar
                rec.reference_id = mem_aln.rid
                rec.reference_start = mem_aln.pos
                rec.mapping_quality = mem_aln.mapq
                cigar = ""
                cigar_len_sum = 0
                for k in range(mem_aln.n_cigar):
                    cigar_op = mem_aln.cigar[k] & 0xf
                    if not opt.soft_clip_supplementary and mem_aln.is_alt == 0 and (
                            cigar_op == 3 or cigar_op == 4):
                        cigar_op = 4 if j > 0 else 3  # // use hard clipping for supplementary alignments
                    cigar_len = mem_aln.cigar[k] >> 4
                    cigar += f"{cigar_len}" + "MIDSH"[cigar_op]
                    if cigar_op < 4:
                        cigar_len_sum += cigar_len
                rec.cigarstring = cigar

                # remove leading and trailing soft-clipped bases for non-primary etc.
                if mem_aln.n_cigar > 0 and j > 0 and not opt.soft_clip_supplementary and mem_aln.is_alt == 0:
                    qb = 0
                    qe = len(query.sequence)
                    leading_op = mem_aln.cigar[0] & 0xf
                    trailing_op = mem_aln.cigar[mem_aln.n_cigar - 1] & 0xf
                    if leading_op == 3 or leading_op == 4:
                        qb += mem_aln.cigar[0] >> 4
                    if trailing_op == 3 or trailing_op == 4:
                        qe -= mem_aln.cigar[mem_aln.n_cigar - 1] >> 4
                    query_qualities = rec.query_qualities # as setting query_sequence may invalidate this
                    if rec.query_sequence is not None:
                        rec.query_sequence = rec.query_sequence[qb:qe]
                    if query_qualities is not None:
                        rec.query_qualities = query_qualities[qb:qe]

                # Optional tags
                attrs = dict()
                if mem_aln.n_cigar > 0:
                    attrs["NM"] = mem_aln.NM
                    md = <char *> (mem_aln.cigar + mem_aln.n_cigar)
                    attrs["MD"] = self.__to_str(md)
                # NB: mate tags are not output: MC, MQ
                if mem_aln.score >= 0:
                    attrs["AS"] = mem_aln.score
                if mem_aln.sub >= 0:
                    attrs["XS"] = mem_aln.sub
                # NB: SA is added after all the records have been created
                if mem_aln.XA != NULL:
                    attrs["XB" if opt.with_xb_tag else "XA"] = mem_aln.XA
                anno = &self._index.bns().anns[rec.reference_id]
                if opt.with_xr_tag and anno.anno != NULL and anno.anno[0] != b'\0':
                    attrs["XR"] = str(anno.anno)
                rec.set_tags(list(attrs.items()))

                mapped_recs.append(rec)

            for j in range(mem_alns.n):
                mem_aln = mem_alns.a[j]
                free(mem_aln.cigar)
                free(mem_aln.XA)
            if len(mapped_recs) == 0:
                recs_to_return.append([self._unmapped(query=query)])
            else:
                self._add_sa_tag(mapped_recs)
                recs_to_return.append(mapped_recs)

        for i in range(num_seqs):
            free(seqs[i].name)
            free(seqs[i].comment)
            free(seqs[i].seq)
            free(seqs[i].qual)
            free(mem_alns_vec[i].a)
        free(seqs)
        free(mem_alns_vec)

        return recs_to_return
