# cython: language_level=3
from pathlib import Path
from typing import List

from fgpyo.sequence import reverse_complement
from libc.string cimport memcpy
from libc.stdlib cimport calloc, free
import enum
from pybwa.libbwaindex cimport BwaIndex
from pysam import FastxRecord, AlignedSegment, qualitystring_to_array, CMATCH, CINS, CDEL, CSOFT_CLIP, CHARD_CLIP
from libc.string cimport strncpy
from pybwa.libbwaindex cimport force_bytes
from pybwa.libbwa cimport bwa_verbose

from pysam.libcalignedsegment cimport makeAlignedSegment


__all__ = [
    "BwaMemMode",
    "BwaMemOptions",
    "BwaMem",
]


cpdef bint _set_bwa_mem_verbosity(int level):
    """Set the BWA C-API verbosity, returning True if changed, false otherwise."""
    global bwa_verbose
    retval = level != bwa_verbose
    bwa_verbose = level
    return retval

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
        threads (int): :code:`bwa mem -t <int>`
        chunk_size (int): :code:`bwa mem -K <int>`
    """

    def _assert_not_finalized(self, attr_name: str) -> None:
        """Raises an AttributeError if the options have been finalized.

        This is used in each setter below to enforce that setters are not used after
         :meth:`~pybwa.BwaMemOptions.finalize` is called.  This could be a property decorator,
         but cython does not support decorating properties.
        """
        if self._finalized:
            raise AttributeError(f"can't set attribute: {attr_name}")

    def __init__(self,
                 min_seed_len: int = 19,
                 mode: BwaMemMode = None,
                 band_width: int = 100,
                 match_score: int = 1,
                 mismatch_penalty: int = 4,
                 minimum_score: int = 30,
                 unpaired_penalty: int = 17,
                 skip_pairing: bool = False,
                 output_all_for_fragments: bool = False,
                 interleaved_paired_end: bool = False,
                 short_split_as_secondary: bool = False,
                 skip_mate_rescue: bool = False,
                 soft_clip_supplementary: bool = False,
                 with_xr_tag: bool = False,
                 query_coord_as_primary: bool = False,
                 keep_mapq_for_supplementary: bool = False,
                 with_xb_tag: bool = False,
                 max_occurrences: int = 500,
                 off_diagonal_x_dropoff: int = 100,
                 ignore_alternate_contigs: bool = False,
                 internal_seed_split_factor: float = 1.5,
                 drop_chain_fraction: float = 0.50,
                 max_mate_rescue_rounds: int = 50,
                 min_seeded_bases_in_chain: int = 0,
                 seed_occurrence_in_3rd_round: int = 20,
                 xa_max_hits: int | tuple[int, int] = (5, 200),
                 xa_drop_ratio: float = 0.80,
                 gap_open_penalty: int | tuple[int, int] = 6,
                 gap_extension_penalty: int | tuple[int, int] = 1,
                 clipping_penalty: int | tuple[int, int] = 5,
                 threads: int = 1,
                 chunk_size: int = 10_000_000) -> None:
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
        self._options0.chunk_size = 0

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
        if not self._finalized:  # pragma: no cover
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
        if self._options.min_seed_len != value:
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
        if self._options.w != value:
            self._options.w = value
            self._options0.w = 1

    @property
    def match_score(self) -> int:
        """:code:`bwa mem -A <int>`"""
        return self._options.a


    @match_score.setter
    def match_score(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.match_score.__name__)
        if self._options.a != value:
            self._options.a = value
            self._options0.a = 1

    @property
    def mismatch_penalty(self) -> int:
        """:code:`bwa mem -B <int>`"""
        return self._options.b

    @mismatch_penalty.setter
    def mismatch_penalty(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.mismatch_penalty.__name__)
        if self._options.b != value:
            self._options.b = value
            self._options0.b = 1

    @property
    def minimum_score(self) -> int:
        """:code:`bwa mem -T <int>`"""
        return self._options.T

    @minimum_score.setter
    def minimum_score(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.minimum_score.__name__)
        if self._options.T != value:
            self._options.T = value
            self._options0.T = 1

    @property
    def unpaired_penalty(self) -> int:
        """:code:`bwa mem -U <int>`"""
        return self._options.pen_unpaired

    @unpaired_penalty.setter
    def unpaired_penalty(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.unpaired_penalty.__name__)
        if self._options.pen_unpaired != value:
            self._options.pen_unpaired = value
            self._options0.pen_unpaired = 1

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
        if self._options.max_occ != value:
            self._options.max_occ = value
            self._options0.max_occ = 1

    @property
    def off_diagonal_x_dropoff(self) -> int:
        """:code:`bwa mem -d <float>`"""
        return self._options.zdrop

    @off_diagonal_x_dropoff.setter
    def off_diagonal_x_dropoff(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.off_diagonal_x_dropoff.__name__)
        if self._options.zdrop != value:
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
        if self._options.split_factor != value:
            self._options.split_factor = value
            self._options0.split_factor = 1

    @property
    def drop_chain_fraction(self) -> float:
        """:code:`bwa mem -D <float>`"""
        return self._options.drop_ratio

    @drop_chain_fraction.setter
    def drop_chain_fraction(self, value: float) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.drop_chain_fraction.__name__)
        if self._options.drop_ratio != value:
            self._options.drop_ratio = value
            self._options0.drop_ratio = 1

    @property
    def max_mate_rescue_rounds(self) -> int:
        """:code:`bwa mem -m <int>`"""
        return self._options.max_matesw

    @max_mate_rescue_rounds.setter
    def max_mate_rescue_rounds(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.max_mate_rescue_rounds.__name__)
        if self._options.max_matesw != value:
            self._options.max_matesw = value
            self._options0.max_matesw = 1

    @property
    def min_seeded_bases_in_chain(self) -> int:
        """:code:`bwa mem -W <int>`"""
        return self._options.min_chain_weight

    @min_seeded_bases_in_chain.setter
    def min_seeded_bases_in_chain(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.min_seeded_bases_in_chain.__name__)
        if self._options.min_chain_weight != value:
            self._options.min_chain_weight = value
            self._options0.min_chain_weight = 1

    @property
    def seed_occurrence_in_3rd_round(self) -> int:
        """:code:`bwa mem -y <int>`"""
        return self._options.max_mem_intv

    @seed_occurrence_in_3rd_round.setter
    def seed_occurrence_in_3rd_round(self, value: int) -> None:
        self._assert_not_finalized(attr_name=BwaMemOptions.seed_occurrence_in_3rd_round.__name__)
        if self._options.max_mem_intv != value:
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
        if isinstance(value, int):
            max_XA_hits = max_XA_hits_alt = value
        else:
            max_XA_hits, max_XA_hits_alt = value
        if self._options.max_XA_hits != max_XA_hits:
            self._options.max_XA_hits = max_XA_hits
            self._options0.max_XA_hits = 1
        if self._options.max_XA_hits_alt != max_XA_hits_alt:
            self._options0.max_XA_hits_alt = 1
            self._options.max_XA_hits_alt = max_XA_hits_alt

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
        if isinstance(value, int):
            deletions = insertions = value
        else:
            deletions, insertions = value
        if self._options.o_del != deletions:
            self._options.o_del = deletions
            self._options0.o_del = 1
        if self._options.o_ins != insertions:
            self._options0.o_ins = 1
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
        if isinstance(value, int):
            deletions = insertions = value
        else:
            deletions, insertions = value
        if self._options.e_del != deletions:
            self._options.e_del = deletions
            self._options0.e_del = 1
        if self._options.e_ins != insertions:
            self._options0.e_ins = 1
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
        if isinstance(value, int):
            five_prime = three_prime = value
        else:
            five_prime, three_prime = value
        if self._options.pen_clip5 != five_prime:
            self._options.pen_clip5 = five_prime
            self._options0.pen_clip5 = 1
        if self._options.pen_clip3 != three_prime:
            self._options0.pen_clip3 = 1
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
        if self._options.chunk_size != value:
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
            end = start + 1
            while end < len(queries) and (num_bases < opt.chunk_size or (end-start)&1 == 1):
                num_bases += len(queries[end].sequence)
                end += 1
            assert start < end
            results.extend(self._calign(opt, queries[start:end], n_processed=start))
            start = end

        return results
    
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

        # qualities
        if q.quality is None:
            s.qual = NULL
        else:
            s.qual = <char *> calloc(sizeof(char), s.l_seq + 1)
            strncpy(s.qual, force_bytes(q.quality), s.l_seq)
            s.qual[s.l_seq] = b'\0'

    cdef _calign(self, opt: BwaMemOptions, queries: List[FastxRecord], n_processed: int = 0):
        # TODO: ignore_alt
        # TODO: refactor to make this more readable
        cdef bseq1_t* seqs
        cdef bseq1_t* seq
        cdef char* s_char
        cdef kstring_t kstr
        cdef int take_all
        cdef size_t j
        cdef mem_alnreg_t *mem_alnreg
        cdef mem_aln_t mem_aln
        cdef char *md
        cdef mem_opt_t *mem_opt
        cdef bntann1_t *anno
        cdef sam_hdr_t *h

        kstr.l = kstr.m = 0
        kstr.s = NULL
        bwa_format_sam_hdr(self._index.bns(), NULL, &kstr)
        h = sam_hdr_parse(kstr.l, kstr.s)
        h.l_text = kstr.l
        h.text = kstr.s

        recs_to_return: List[List[AlignedSegment]] = []

        # copy FastxRecord into bwa_seq_t
        num_seqs = len(queries)
        mem_opt = opt.mem_opt()
        seqs = <bseq1_t*>calloc(sizeof(bseq1_t), num_seqs)
        for i in range(num_seqs):
            self._copy_seq(queries[i], &seqs[i])

        # process the sequences (ignores the paired end stats)
        mem_process_seqs(mem_opt, self._index.bwt(), self._index.bns(), self._index.pac(),
                         n_processed, num_seqs, seqs, NULL, h)

        for i in range(num_seqs):
            seq = &seqs[i]
            recs = []
            for j in range(seq.bams.l):
                rec = makeAlignedSegment(seq.bams.bams[j], self._index.header)
                recs.append(rec)
            recs_to_return.append(recs)
            bams_destroy(seq.bams)
            seq.bams = NULL

        for i in range(num_seqs):
            free(seqs[i].name)
            free(seqs[i].comment)
            free(seqs[i].seq)
            free(seqs[i].qual)
        free(seqs)
        free(kstr.s)

        return recs_to_return


###########################################################
# Below is code to facilitate testing
###########################################################


from pybwa.libbwamem cimport BwaMemOptions
from pybwa.libbwamem cimport mem_opt_t
from pybwa.libbwamem cimport mem_opt_init
from libc.stdlib cimport free

cdef _assert_mem_opt_are_the_same_c():
    cdef mem_opt_t *bwa_mem_opt
    cdef mem_opt_t *pybwa_mem_opt
    bwa_mem_opt = mem_opt_init()
    options = BwaMemOptions()
    pybwa_mem_opt = options._options

    # Check member by member for the mem_opt_t struct
    assert pybwa_mem_opt.a == bwa_mem_opt.a
    assert pybwa_mem_opt.b == bwa_mem_opt.b
    assert pybwa_mem_opt.o_del == bwa_mem_opt.o_del
    assert pybwa_mem_opt.e_del == bwa_mem_opt.e_del
    assert pybwa_mem_opt.o_ins == bwa_mem_opt.o_ins
    assert pybwa_mem_opt.e_ins == bwa_mem_opt.e_ins
    assert pybwa_mem_opt.pen_unpaired == bwa_mem_opt.pen_unpaired
    assert pybwa_mem_opt.pen_clip5 == bwa_mem_opt.pen_clip5
    assert pybwa_mem_opt.pen_clip3 == bwa_mem_opt.pen_clip3
    assert pybwa_mem_opt.w == bwa_mem_opt.w
    assert pybwa_mem_opt.zdrop == bwa_mem_opt.zdrop
    assert pybwa_mem_opt.max_mem_intv == bwa_mem_opt.max_mem_intv
    assert pybwa_mem_opt.T == bwa_mem_opt.T
    assert pybwa_mem_opt.flag == bwa_mem_opt.flag
    assert pybwa_mem_opt.min_seed_len == bwa_mem_opt.min_seed_len
    assert pybwa_mem_opt.min_chain_weight == bwa_mem_opt.min_chain_weight
    assert pybwa_mem_opt.max_chain_extend == bwa_mem_opt.max_chain_extend
    assert pybwa_mem_opt.split_factor == bwa_mem_opt.split_factor
    assert pybwa_mem_opt.split_width == bwa_mem_opt.split_width
    assert pybwa_mem_opt.max_occ == bwa_mem_opt.max_occ
    assert pybwa_mem_opt.max_chain_gap == bwa_mem_opt.max_chain_gap
    assert pybwa_mem_opt.n_threads == bwa_mem_opt.n_threads
    assert pybwa_mem_opt.chunk_size == bwa_mem_opt.chunk_size
    assert pybwa_mem_opt.mask_level == bwa_mem_opt.mask_level
    assert pybwa_mem_opt.drop_ratio == bwa_mem_opt.drop_ratio
    assert pybwa_mem_opt.XA_drop_ratio == bwa_mem_opt.XA_drop_ratio
    assert pybwa_mem_opt.mask_level_redun == bwa_mem_opt.mask_level_redun
    assert pybwa_mem_opt.mapQ_coef_len == bwa_mem_opt.mapQ_coef_len
    assert pybwa_mem_opt.max_ins == bwa_mem_opt.max_ins
    assert pybwa_mem_opt.max_matesw == bwa_mem_opt.max_matesw
    assert pybwa_mem_opt.max_XA_hits == bwa_mem_opt.max_XA_hits
    assert pybwa_mem_opt.max_XA_hits_alt == bwa_mem_opt.max_XA_hits_alt

    free(bwa_mem_opt)


def _assert_mem_opt_are_the_same() -> None:
    """Tests that the defaults are synced between bwa and pybwa."""
    _assert_mem_opt_are_the_same_c()


cdef _call_mem_opt_when_not_finalized_c():
    options = BwaMemOptions()
    assert options.mem_opt() != NULL  # this should except since it is not finalized


def _call_mem_opt_when_not_finalized() -> None:
    _call_mem_opt_when_not_finalized_c()
