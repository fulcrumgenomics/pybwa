# cython: language_level=3
from pathlib import Path
from typing import List
from libc.string cimport memset, memcpy
from libc.stdlib cimport calloc, free
import enum
from bwapy.libbwaindex cimport BwaIndex
from pysam import FastxRecord, AlignedSegment

__all__ = [
    "BwaMemMode",
    "BwaMemOptions",
    "BwaMemOptionsBuilder",
    "BwaMem",
]


# class syntax
@enum.unique
class BwaMemMode(enum.Enum):
    """The read type for overriding multiple options"""
    PACBIO = enum.auto()
    ONT2D = enum.auto()
    INTRACTG = enum.auto()


cdef class BwaMemOptions:
    """The container for options for [`BwaMem`][bwapy.BwaMem]."""
    _ignore_alt: bool
    _mode: BwaMemMode | None
    cdef mem_opt_t* _options

    def __init__(self):
        self._ignore_alt = False
        self._mode = None

    def __cinit__(self):
        self._options = mem_opt_init()

    def __dealloc__(self):
        free(self._options)
        self._options = NULL

    cdef mem_opt_t* mem_opt(self):
        return self._options

    property min_seed_len:
        """bwa mem -k <int>"""
        def __get__(self):
            return self._options.min_seed_len

    property mode:
        """bwa mem -x <str>"""
        def __get__(self) -> BwaMemMode:
            return self._mode

    property band_width:
        """bwa mem -w <int>"""
        def __get__(self):
            return self._options.w

    property match_score:
        """bwa mem -A <int>"""
        def __get__(self):
            return self._options.a

    property mismatch_penalty:
        """bwa mem -A <int>"""
        def __get__(self):
            return self._options.b

    property minimum_score:
        """bwa mem -T <int>"""
        def __get__(self):
            return self._options.T

    property unpaired_penalty:
        """bwa mem -U <int>"""
        def __get__(self):
            return self._options.pen_unpaired

    property n_threads:
        """bwa mem -t <int>"""
        def __get__(self):
            return self._options.n_threads

    property skip_pairing:
        def __get__(self):
            return (self._options.flag & MEM_F_NOPAIRING) != 0

    property output_all_for_fragments:
        def __get__(self):
            return (self._options.flag & MEM_F_ALL) != 0

    property interleaved_paired_end:
        def __get__(self):
            return (self._options.flag & (MEM_F_PE | MEM_F_SMARTPE)) != 0

    property skip_mate_rescue:
        def __get__(self):
            return (self._options.flag & MEM_F_NO_MULTI) != 0

    property soft_clip_supplementary:
        def __get__(self):
            return (self._options.flag & MEM_F_SOFTCLIP) != 0

    property with_xr_tag:
        def __get__(self):
            return (self._options.flag & MEM_F_REF_HDR) != 0

    property query_coord_as_primary:
        def __get__(self):
            return (self._options.flag & (MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ)) != 0

    property keep_mapq_for_supplementary:
        def __get__(self):
            return (self._options.flag & MEM_F_KEEP_SUPP_MAPQ) != 0

    property with_xb_tag:
        def __get__(self):
            return (self._options.flag & MEM_F_XB) != 0

    property max_occurrences:
        """bwa mem -c <int>"""
        def __get__(self):
            return self._options.max_occ

    property off_diagonal_x_dropoff:
        """bwa mem -d <float>"""
        def __get__(self):
            return self._options.XA_drop_ratio

    property ignore_alternate_contigs:
        """bwa mem -j"""
        def __get__(self):
            return self._ignore_alt

    property internal_seed_split_factor:
        """bwa mem -r <float>"""
        def __get__(self):
            return self._options.split_factor

    property drop_chain_fraction:
        """bwa mem -D <float>"""
        def __get__(self):
            return self._options.drop_ratio
        def __set__(self, value: float):
            self._options.drop_ratio = value
            self._options0.drop_ratio = 1

    property max_mate_rescue_rounds:
        """bwa mem -m <int>"""
        def __get__(self):
            return self._options.max_matesw

    property min_seeded_bases_in_chain:
        """bwa mem -W <int>"""
        def __get__(self):
            return self._options.min_chain_weight

    property seed_occurrence_in_3rd_round:
        """bwa mem -y <int>"""
        def __get__(self):
            return self._options.max_mem_intv

    property xa_max_hits:
        """bwa mem -h <int<,int>>"""
        def __get__(self):
            return self._options.max_XA_hits, self._options.max_XA_hits_alt

    property xa_drop_ration:
        """bwa mem -y <float>"""
        def __get__(self):
            return self._options.XA_drop_ratio

    property gap_open_penalty:
        """bwa mem -O <int<,int>>"""
        def __get__(self):
            if self._options.o_del == self._options.o_ins:
                return self._options.o_del
            else:
                return self._options.o_del, self._options.o_ins

    property gap_extension_penalty:
        """bwa mem -E <int<,int>>"""
        def __get__(self):
            if self._options.e_del == self._options.e_ins:
                return self._options.e_del
            else:
                return self._options.e_del, self._options.e_ins


    property clipping_penalty:
        """bwa mem -L <int<,int>>"""
        def __get__(self):
            if self._options.pen_clip5 == self._options.pen_clip3:
                return self._options.pen_clip5
            else:
                return self._options.pen_clip5, self._options.pen_clip3



cdef class BwaMemOptionsBuilder(BwaMemOptions):
    """The container for options for [`BwaMem`][bwapy.BwaMem]."""
    cdef mem_opt_t* _options0

    def __init__(self, options: BwaMemOptions | None = None):
        super().__init__()
        if options is not None:
            self.ignore_alternate_contigs = options.ignore_alternate_contigs
            self.mode = options.mode
            self._copy_options(options)

    cdef _copy_options(self, options: BwaMemOptions):
        memcpy(self._options, options.mem_opt(), sizeof(mem_opt_t))

    def __cinit__(self, options: BwaMemOptions | None = None):
        self._options0 = mem_opt_init()
        memset(self._options0, 0, sizeof(mem_opt_t))

    def __dealloc__(self):
        free(self._options0)
        self._options0 = NULL

    def build(self) -> BwaMemOptions:
        if self._mode is None:
            # matching score is changed so scale the rest of the penalties/scores
            if self._options0.a == 1:
                if self._options0.b != 1:
                    self._options.b *= self._options.a
                if self._options0.T != 1:
                    self._options.T *= self._options.a
                if self._options0.o_del != 1:
                    self._options.o_del *= self._options.a
                if self._options0.e_del != 1:
                    self._options.e_del *= self._options.a
                if self._options0.o_ins != 1:
                    self._options.o_ins *= self._options.a
                if self._options0.e_ins != 1:
                    self._options.e_ins *= self._options.a
                if self._options0.zdrop != 1:
                    self._options.zdrop *= self._options.a
                if self._options0.pen_clip5 != 1:
                    self._options.pen_clip5 *= self._options.a
                if self._options0.pen_clip3 != 1:
                    self._options.pen_clip3 *= self._options.a
                if self._options0.pen_unpaired != 1:
                    self._options.pen_unpaired *= self._options.a
        elif self._mode == BwaMemMode.INTRACTG:
            if self._options0.o_del != 1:
                self._options.o_del = 16
            if self._options0.o_ins != 1:
                self._options.o_ins = 16
            if self._options0.b != 1:
                self._options.b = 9
            if self._options0.pen_clip5 != 1:
                self._options.pen_clip5 = 5
            if self._options0.pen_clip3 != 1:
                self._options.pen_clip3 = 5
        else:
            if self._options0.o_del != 1:
                self._options.o_del = 1
            if self._options0.o_ins != 1:
                self._options.o_ins = 1
            if self._options0.e_del != 1:
                self._options.e_del = 1
            if self._options0.e_ins != 1:
                self._options.e_ins = 1
            if self._options0.b != 1:
                self._options.b = 1
            if self._options0.split_factor == 0.0:
                self._options.split_factor = 10.0
            if self._options0.pen_clip5 != 1:
                self._options.pen_clip5 = 0
            if self._options0.pen_clip3 != 1:
                self._options.pen_clip3 = 0
            # ONT2D vs PACBIO options
            if self._options0.min_chain_weight != 1:
                self._options.min_chain_weight = 20 if self._mode == BwaMemMode.ONT2D else 40
            if self._options0.min_seed_len != 1:
                self._options.min_seed_len = 14 if self._mode == BwaMemMode.ONT2D else 17

        bwa_fill_scmat(
            self._options.a, self._options.b, self._options.mat
        )
        return self

    property min_seed_len:
        """bwa mem -k <int>"""
        def __get__(self):
            return super().min_seed_len
        def __set__(self, value: int):
            self._options.min_seed_len = value
            self._options0.min_seed_len = 1

    property mode:
        """bwa mem -x <str>"""
        def __get__(self) -> BwaMemMode:
            return super().mode
        def __set__(self, value: BwaMemMode):
            self._mode = value

    property band_width:
        """bwa mem -w <int>"""
        def __get__(self):
            return super().band_width
        def __set__(self, value: int):
            self._options.w = value
            self._options0.w = 1

    property match_score:
        """bwa mem -A <int>"""
        def __get__(self):
            return super().match_score
        def __set__(self, value: int):
            self._options.a = value
            self._options0.a = 1

    property mismatch_penalty:
        """bwa mem -A <int>"""
        def __get__(self):
            return super().mismatch_penalty
        def __set__(self, value: int):
            self._options.b = value
            self._options0.b = 1

    property minimum_score:
        """bwa mem -T <int>"""
        def __get__(self):
            return super().minimum_score
        def __set__(self, value: int):
            self._options.T = value
            self._options0.T = 1

    property unpaired_penalty:
        """bwa mem -U <int>"""
        def __get__(self):
            return super().unpaired_penalty
        def __set__(self, value: int):
            self._options.pen_unpaired = value
            self._options0.pen_unpaired = 1

    property n_threads:
        """bwa mem -t <int>"""
        def __get__(self):
            return super().n_threads
        def __set__(self, value: int):
            self._options.n_threads = value if value > 1 else 1

    def _set_flag(self, value: bool, flag: int):
        if value:
            self._options.flag |= flag
        else:
            self._options.flag |= ~flag
        return self

    property skip_pairing:
        def __get__(self):
            return super().skip_pairing
        def __set__(self, value: bool):
            self._set_flag(value, MEM_F_NOPAIRING)

    property output_all_for_fragments:
        def __get__(self):
            return super().output_all_for_fragments
        def __set__(self, value: bool):
            self._set_flag(value, MEM_F_ALL)

    property interleaved_paired_end:
        def __get__(self):
            return super().interleaved_paired_end
        def __set__(self, value: bool):
            self._set_flag(value, MEM_F_PE | MEM_F_SMARTPE)

    property skip_mate_rescue:
        def __get__(self):
            return super().skip_pairing
        def __set__(self, value: bool):
            self._set_flag(value, MEM_F_NO_MULTI)

    property soft_clip_supplementary:
        def __get__(self):
            return super().soft_clip_supplementary
        def __set__(self, value: bool):
            self._set_flag(value, MEM_F_SOFTCLIP)

    property with_xr_tag:
        def __get__(self):
            return super().with_xr_tag
        def __set__(self, value: bool):
            self._set_flag(value, MEM_F_REF_HDR)

    property query_coord_as_primary:
        def __get__(self):
            return super().query_coord_as_primary
        def __set__(self, value: bool):
            self._set_flag(value, MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ) #/ always apply MEM_F_KEEP_SUPP_MAPQ with -5

    property keep_mapq_for_supplementary:
        def __get__(self):
            return super().keep_mapq_for_supplementary
        def __set__(self, value: bool):
            self._set_flag(value, MEM_F_KEEP_SUPP_MAPQ)

    property with_xb_tag:
        def __get__(self):
            return super().with_xb_tag
        def __set__(self, value: bool):
            self._set_flag(value, MEM_F_XB)

    property max_occurrences:
        """bwa mem -c <int>"""
        def __get__(self):
            return super().max_occurrences
        def __set__(self, value: int):
            self._options.max_occ = value
            self._options0.max_occ = 1

    property off_diagonal_x_dropoff:
        """bwa mem -d <float>"""
        def __get__(self):
            return super().off_diagonal_x_dropoff
        def __set__(self, value: float):
            self._options.XA_drop_ratio = value
            self._options0.XA_drop_ratio = 1

    property ignore_alternate_contigs:
        """bwa mem -j"""
        def __get__(self):
            return super()._ignore_alt
        def __set__(self, value: bool):
           self._ignore_alt = value

    property internal_seed_split_factor:
        """bwa mem -r <float>"""
        def __get__(self):
            return super().internal_seed_split_factor
        def __set__(self, value: float):
            self._options.split_factor = value
            self._options0.split_factor = 1

    property drop_chain_fraction:
        """bwa mem -D <float>"""
        def __get__(self):
            return super().drop_chain_fraction
        def __set__(self, value: float):
            self._options.drop_ratio = value
            self._options0.drop_ratio = 1

    property max_mate_rescue_rounds:
        """bwa mem -m <int>"""
        def __get__(self):
            return super().max_mate_rescue_rounds
        def __set__(self, value: int):
            self._options.max_matesw = value
            self._options0.max_matesw = 1

    property min_seeded_bases_in_chain:
        """bwa mem -W <int>"""
        def __get__(self):
            return super().min_seeded_bases_in_chain
        def __set__(self, value: int):
            self._options.min_chain_weight = value
            self._options0.min_chain_weight = 1

    property seed_occurrence_in_3rd_round:
        """bwa mem -y <int>"""
        def __get__(self):
            return super().seed_occurrence_in_3rd_round
        def __set__(self, value: int):
            self._options.max_mem_intv = value
            self._options0.max_mem_intv = 1

    property xa_max_hits:
        """bwa mem -h <int<,int>>"""
        def __get__(self):
            return super().xa_max_hits
        def __set__(self, value: int | tuple[int, int]):
            self._options0.max_XA_hits = 1
            self._options0.max_XA_hits_alt = 1
            if isinstance(value, int):
                self._options.max_XA_hits = value
                self._options.max_XA_hits_alt = value
            else:
                left, right = value
                self._options.max_XA_hits = left
                self._options.max_XA_hits_alt = right

    property xa_drop_ration:
        """bwa mem -y <float>"""
        def __get__(self):
            return self._options.XA_drop_ratio
        def __set__(self, value: float):
           self._options.XA_drop_ratio = value

    property gap_open_penalty:
        """bwa mem -O <int<,int>>"""
        def __get__(self):
            return super().gap_open_penalty
        def __set__(self, value: int | tuple[int, int]):
            self._options0.o_del = 1
            self._options0.o_ins = 1
            if isinstance(value, int):
                self._options.o_del = value
                self._options.o_ins = value
            else:
                deletions, insertions = value
                self._options.o_del = deletions
                self._options.o_ins = insertions

    property gap_extension_penalty:
        """bwa mem -E <int<,int>>"""
        def __get__(self):
            return super().gap_extension_penalty
        def __set__(self, value: int | tuple[int, int]):
            self._options0.e_del = 1
            self._options0.e_ins = 1
            if isinstance(value, int):
                self._options.e_del = value
                self._options.e_ins = value
            else:
                deletions, insertions = value
                self._options.e_del = deletions
                self._options.e_ins = insertions


    property clipping_penalty:
        """bwa mem -L <int<,int>>"""
        def __get__(self):
            return super().clipping_penalty
        def __set__(self, value: int | tuple[int, int]):
            self._options0.pen_clip5 = 1
            self._options0.pen_clip3 = 1
            if isinstance(value, int):
                self._options.pen_clip5 = value
                self._options.pen_clip3 = value
            else:
                five_prime, three_prime = value
                self._options.pen_clip5 = five_prime
                self._options.pen_clip3 = three_prime


cdef class BwaMem:
    """The class to align reads with `bwa mem`."""

    cdef BwaIndex _index

    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None):
        if prefix is not None:
            assert Path(prefix).exists()
            self._index = BwaIndex(prefix=prefix)
        elif index is not None:
            self._index = index
        else:
            raise Exception("Either prefix or index must be given")

    # TODO: support paired end
    def align(self, opt: BwaMemOptions, queries: List[FastxRecord]) -> List[List[AlignedSegment]]:
        """Align one or more queries with `bwa aln`.

        Args:
            opt: the alignment options
            queries: the queries to align

        Returns:
            one alignment per query
        """
        return self._calign(opt, queries)

    cdef _copy_seq(self, q: FastxRecord, kstring_t *seq):
        seq_len = len(q.sequence)
        seq.s = <char *> calloc(sizeof(char), seq_len + 1)
        seq.l = seq_len
        seq.m = seq_len + 1
        for i, base in enumerate(q.sequence):
            seq.s[i] = nst_nt4_table[ord(base)]
        seq.s[seq_len] = b'\0'

    def _unmapped(self, query: FastxRecord) -> AlignedSegment:
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
        return rec

    # mimics mem_aln2sam in bwamem.c
    cdef _build_alignment(self, opt: BwaMemOptions, query: FastxRecord, mem_aln_t *mem_aln_v, int mem_aln_n, int which):
        cdef mem_aln_t mem_aln
        cdef mem_aln_t* mem_aln_other

        # create a AlignedSegment record here
        rec = self._unmapped(query=query)

        mem_aln = mem_aln_v[which]

        # set the flags
        mem_aln.flag |= 0x4 if mem_aln.rid < 0 else 0
        mem_aln.flag |= 0x10 if mem_aln.is_rev > 0 else 0
        rec.flag = (mem_aln.flag & 0xffff) | (0x100 if (mem_aln.flag & 0x10000) != 0 else 0)

        # sequence and qualities
        rec.query_sequence = query.sequence if rec.is_forward else query.sequence[::-1]
        if query.quality is not None:
            rec.query_qualities = query.quality if rec.is_forward else query.quality[::-1]

        if rec.is_unmapped:
            return rec

        # reference id, position, mapq, and cigar
        rec.reference_id = mem_aln.rid
        rec.reference_start = mem_aln.pos
        rec.mapping_quality = mem_aln.mapq
        cigar = ""
        for i in range(mem_aln.n_cigar):
            cigar_op = mem_aln.cigar[i] & 0xf
            if opt.soft_clip_supplementary and mem_aln.is_alt == 0 and (
                    cigar_op == 3 or cigar_op == 4):
                cigar_op = 4 if mem_aln_n > 0 else 3  # // use hard clipping for supplementary alignments
            cigar_len = mem_aln.cigar[i] >> 4
            cigar = cigar + f"{cigar_len}" + "MIDS"[cigar_op]
        rec.cigarstring = cigar

        # remove leading and trailing soft-clipped bases for non-primary etc.
        if mem_aln.n_cigar > 0 and mem_aln_n > 0 and not opt.soft_clip_supplementary and not mem_aln.is_alt:
            qb = 0
            qe = len(query.sequence)
            if (mem_aln.cigar[0] & 0xf) == 4 or (mem_aln.cigar[0] & 0xf) == 3:
                qb += mem_aln.cigar[0] >> 4
            if (mem_aln.cigar[mem_aln.n_cigar - 1] & 0xf) == 4 or (
                    mem_aln.cigar[mem_aln.n_cigar - 1] & 0xf) == 3:
                qe -= mem_aln.cigar[mem_aln.n_cigar - 1] >> 4
            rec.query_sequence = rec.query_sequence[qb:qe]
            if query.quality is not None:
                rec.query_qualities = rec.query_qualities[qb:qe]

        # Optional tags
        md = <char *> (mem_aln.cigar + mem_aln.n_cigar)
        attrs = dict()
        attrs["NM"] = f"{mem_aln.NM}"
        attrs["MD"] = f"{md}"
        # NB: mate tags are not output: MC, MQ
        if mem_aln.score >= 0:
            attrs["AS"] = mem_aln.score
        if mem_aln.sub >= 0:
            attrs["XS"] = mem_aln.sub
        # if mem_aln.flag & 0x100 == 0:  # secondary
        #     # find if there are other primary hits, and if so, output them in the SA tag
        #     for i in range(mem_aln_n):
        #         if i == which or mem_aln_v[i].flag & 0x100 == 0:
        #             break
        #     if i < mem_aln_n:
        #         for i in range(mem_aln_n):
        #             mem_aln_other = &mem_aln_v[i]
        #             if i == which or mem_aln_other.flag & 0x100 != 0:
        #                 continue
        #             SA = self._index.bns().anns[mem_aln_other.rid].name + ","
        #             SA += f"{mem_aln_other.pos + 1},"
        #             SA += "+-"[mem_aln_other.is_rev] + ","
        #             for k in range(mem_aln_other.n_cigar):
        #                 cigar_op = mem_aln_other.cigar[i] & 0xf
        #                 cigar_len = mem_aln_other.cigar[i] >> 4
        #                 SA += f"{cigar_len}" + "MIDSH"[cigar_op]
        #             SA += f",{mem_aln_other.mapq}"
        #             SA += f",{mem_aln_other.NM}"
        #             SA += ";"
        #         attrs["SA"] = SA
        #     if mem_aln.alt_sc > 0:
        #         attrs['pa'] = mem_aln.score / float(mem_aln.alt_sc)
        # if mem_aln.XA != NULL:
        #     attrs["XB" if opt.with_xb_tag else "XA"] = mem_aln.XA
        # if opt.with_xr_tag and self._index.bns().anns[rec.reference_id].anno != 0 and \
        #         self._index.bns().anns[rec.reference_id].anno[0] != 0:
        #     attrs["XR"] = self._index.bns().anns[rec.reference_id].anno
        rec.set_tags(list(attrs.items()))

        return rec

    cdef _calign(self, opt: BwaMemOptions, queries: List[FastxRecord]):
        # TODO: ignore_alt
        # TODO: refactor to make this more readable
        cdef kstring_t* seqs
        cdef kstring_t* seq
        cdef char* s_char
        cdef kstring_t* kstr
        cdef mem_alnreg_v mem_alnregs
        cdef int take_all
        cdef size_t j
        cdef char **XA
        cdef mem_alnreg_t *mem_alnreg
        cdef mem_aln_t mem_aln
        cdef mem_aln_t *mem_aln_v
        cdef int mem_aln_n
        cdef char *md
        cdef mem_opt_t *mem_opt
        print(f"bwamem [start] num_seqs: {len(queries)}")

        mem_opt = opt.mem_opt()

        recs_to_return: List[List[AlignedSegment]] = []

        # copy FastqProxy into bwa_seq_t
        num_seqs = len(queries)
        print(f"bwamem [start] num_seqs: {num_seqs}")
        seqs = <kstring_t*>calloc(sizeof(kstring_t), num_seqs)
        for i in range(num_seqs):
            print(f"bwamem [start] {i}")
            seq = &seqs[i]
            query = queries[i]
            self._copy_seq(queries[i], seq)
            print(f"bwamem [mem_align1] i: {i}")
            mem_alnregs = mem_align1(mem_opt, self._index.bwt(), self._index.bns(), self._index.pac(), seq.l, seq.s)
            if opt.query_coord_as_primary:
                print(f"bwamem [mem_reorder_primary5] i: {i}")
                mem_reorder_primary5(opt.minimum_score, &mem_alnregs)

            # mimic mem_reg2sam from bwamem.c
            mem_aln_n = 0
            mem_aln_v = <mem_aln_t*>calloc(sizeof(mem_aln_t), mem_alnregs.n)

            recs = []
            XA = NULL
            keep_all = opt.output_all_for_fragments
            if not keep_all:
                print(f"bwamem [mem_gen_alt] i: {i}")
                XA = mem_gen_alt(mem_opt, self._index.bns(), self._index.pac(), &mem_alnregs, seq.l, seq.s)
            for j in range(mem_alnregs.n):
                mem_alnreg = &mem_alnregs.a[j]

                print(f"bwamem [loop] i: {i} j: {j}")
                if mem_alnreg.score < opt.minimum_score:
                    continue
                if mem_alnreg.secondary >= 0 and (mem_alnreg.is_alt or not keep_all):
                    continue
                if 0 <= mem_alnreg.secondary < INT_MAX and mem_alnreg.score < mem_alnregs.a[mem_alnreg.secondary].score * opt.xa_drop_ration:
                    continue
                print(f"bwamem [mem_reg2aln] i: {i} j: {j}")
                mem_aln = mem_reg2aln(mem_opt, self._index.bns(), self._index.pac(), seq.l, seq.s, mem_alnreg)
                mem_aln.XA = XA[j] if not keep_all else NULL
                if mem_alnreg.secondary >= 0:
                    mem_aln.sub = -1  # don't output sub-optimal score
                if mem_aln_n > 0 > mem_alnreg.secondary:  # if supplementary
                    mem_aln.flag |= 0x10000 if opt.skip_mate_rescue else 0x800
                if opt.keep_mapq_for_supplementary and mem_aln_n > 0 and mem_alnreg.is_alt > 0 and mem_aln.mapq > mem_alnregs.a[0].mapq:
                    mem_aln.mapq = mem_alnregs.a[0].mapq  # lower

                mem_aln_v[mem_aln_n] = mem_aln
                mem_aln_n += 1

            print(f"bwamem [after] i: {i} mem_aln_n: {mem_aln_n}")
            if mem_aln_n == 0:
                recs = [self._unmapped(query=queries[i])]
            else:
                print(f"bwamem [_build_alignment] i: {i} mem_aln_n: {mem_aln_n}")
                recs = [
                    self._build_alignment(opt, query, mem_aln_v, mem_aln_n, j)
                    for j in range(mem_aln_n)
                ]
            print(f"bwamem [freeing] i: {i} mem_aln_n: {mem_aln_n}")
            for j in range(mem_aln_n):
                free(XA[j])
                free(mem_aln_v[j].cigar)
            free(XA)
            free(mem_aln_v)
            free(mem_alnregs.a)
            recs_to_return.append(recs)

        print(f"bwamem [freeing]")
        for i in range(num_seqs):
            free(seqs[i].s)
        free(seqs)

        return recs_to_return