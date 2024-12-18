# cython: language_level=3
from ctypes import memset
from pathlib import Path
from typing import List

from libc.stdlib cimport calloc, free, malloc
import enum
from bwapy.libbwaindex cimport force_bytes
from bwapy.libbwaindex cimport BwaIndex
from pysam import FastxRecord, AlignedSegment

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
    cdef mem_opt_t _options0  # attributes are set to `1` when options are set in _options
    _mode: BwaMemMode | None

    def __init__(self, options: BwaMemOptions | None = None) -> None:
        self._options: BwaMemOptions = BwaMemOptions() if options is None else options
        self._mode = None
        self.__cinit__()

    cdef _cinit(self):
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
                    self._options._delegate.b *= self._options._delegate.a
                if self._options0.T != 1:
                    self._options._delegate.T *= self._options._delegate.a
                if self._options0.o_del != 1:
                    self._options._delegate.o_del *= self._options._delegate.a
                if self._options0.e_del != 1:
                    self._options._delegate.e_del *= self._options._delegate.a
                if self._options0.o_ins != 1:
                    self._options._delegate.o_ins *= self._options._delegate.a
                if self._options0.e_ins != 1:
                    self._options._delegate.e_ins *= self._options._delegate.a
                if self._options0.zdrop != 1:
                    self._options._delegate.zdrop *= self._options._delegate.a
                if self._options0.pen_clip5 != 1:
                    self._options._delegate.pen_clip5 *= self._options._delegate.a
                if self._options0.pen_clip3 != 1:
                    self._options._delegate.pen_clip3 *= self._options._delegate.a
                if self._options0.pen_unpaired != 1:
                    self._options._delegate.pen_unpaired *= self._options._delegate.a
        elif self._mode == BwaMemMode.INTRACTG:
            if self._options0.o_del != 1:
                self._options._delegate.o_del = 16
            if self._options0.o_ins != 1:
                self._options._delegate.o_ins = 16
            if self._options0.b != 1:
                self._options._delegate.b = 9
            if self._options0.pen_clip5 != 1:
                self._options._delegate.pen_clip5 = 5
            if self._options0.pen_clip3 != 1:
                self._options._delegate.pen_clip3 = 5
        else:
            if self._options0.o_del != 1:
                self._options._delegate.o_del = 1
            if self._options0.o_ins != 1:
                self._options._delegate.o_ins = 1
            if self._options0.e_del != 1:
                self._options._delegate.e_del = 1
            if self._options0.e_ins != 1:
                self._options._delegate.e_ins = 1
            if self._options0.b != 1:
                self._options._delegate.b = 1
            if self._options0.split_factor == 0.0:
                self._options._delegate.split_factor = 10.0
            if self._options0.pen_clip5 != 1:
                self._options._delegate.pen_clip5 = 0
            if self._options0.pen_clip3 != 1:
                self._options._delegate.pen_clip3 = 0
            # ONT2D vs PACBIO options
            if self._options0.min_chain_weight != 1:
                self._options._delegate.min_chain_weight = 20 if self._mode == BwaMemMode.ONT2D else 40
            if self._options0.min_seed_len != 1:
                self._options._delegate.min_seed_len = 14 if self._mode == BwaMemMode.ONT2D else 17

        bwa_fill_scmat(
            self._options._delegate.a, self._options._delegate.b, self._options._delegate.mat
        )

        return self._options

    cdef min_seed_len(self, value: int):
        """bwa mem -k <int>"""
        self._options._delegate.min_seed_len = value
        self._options0.min_seed_len = 1
        return self

    # FIXME: convert to an enum enum
    cdef mode(self, value: BwaMemMode):
        """bwa mem -x <str>"""
        self._mode = value
        return self

    cdef band_width(self, value: int):
        """bwa mem -w <int>"""
        self._options._delegate.w = value
        self._options0.w = 1
        return self

    cdef match_score(self, value: int):
        """bwa mem -A <int>"""
        self._options._delegate.a = value
        self._options0.a = 1
        return self

    cdef mismatch_penalty(self, value: int):
        """bwa mem -A <int>"""
        self._options._delegate.b = value
        self._options0.b = 1
        return self

    cdef minimum_score(self, value: int):
        """bwa mem -T <int>"""
        self._options._delegate.T = value
        self._options0.T = 1
        return self

    cdef unpaired_penalty(self, value: int):
        """bwa mem -U <int>"""
        self._options._delegate.pen_unpaired = value
        self._options0.pen_unpaired = 1
        return self

    cdef num_threads(self, value: int):
        """bwa mem -t <int>"""
        self._options._delegate.n_threads = value if value > 1 else 1
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
        self._options0.max_occ = 1
        return self

    cdef off_diagonal_x_dropoff(self, value: int):
        """bwa mem -d <float>"""
        self._options._delegate.XA_drop_ratio = value
        self._options0.XA_drop_ratio = 1
        return self

    cdef ignore_alternate_contigs(self, value: bool):
        """bwa mem -j"""
        self._options._ignore_alt = value
        return self

    cdef internal_seed_split_factor(self, value: float):
        """bwa mem -r <float>"""
        self._options._delegate.split_factor = value
        self._options0.split_factor = 1
        return self

    cdef drop_chain_fraction(self, value: float):
        """bwa mem -D <float>"""
        self._options._delegate.drop_ratio = value
        self._options0.drop_ratio = 1
        return self

    cdef max_mate_rescue_rounds(self, value: int):
        """bwa mem -m <float>"""
        self._options._delegate.max_matesw = value
        self._options0.max_matesw = 1
        return self

    cdef min_seeded_bases_in_chain(self, value: int):
        """bwa mem -W <float>"""
        self._options._delegate.min_chain_weight = value
        self._options0.min_chain_weight = 1
        return self

    cdef seed_occurrence_in_3rd_round(self, value: int):
        """bwa mem -y <float>"""
        self._options._delegate.max_mem_intv = value
        self._options0.max_mem_intv = 1
        return self

    cdef xa_max_hits(self, value: int, alt_value: int | None):
        """bwa mem -h <int<,int>>"""
        self._options0.max_XA_hits = 1
        self._options0.max_XA_hits_alt = 1
        self._options._delegate.max_XA_hits = value
        self._options._delegate.max_XA_hits_alt = value if alt_value is None else alt_value

    cdef xa_drop_ration(self, value: float):
        """bwa mem -y <float>"""
        self._options._delegate.XA_drop_ratio = value
        return self

    cdef gap_open_penalty(self, deletions: int, insertions: int):
        """bwa mem -O <int<,int>>"""
        self._options0.o_del = 1
        self._options0.o_ins = 1
        self._options._delegate.o_del = deletions
        self._options._delegate.o_ins = insertions
        return self

    cdef gap_extension_penalty(self, deletions: int, insertions: int):
        """bwa mem -E <int<,int>>"""
        self._options0.e_del = 1
        self._options0.e_ins = 1
        self._options._delegate.e_del = deletions
        self._options._delegate.e_ins = insertions
        return self

    cdef clipping_penalty(self, five_prime: int, three_prime: int):
        """bwa mem -L <int<,int>>"""
        self._options0.pen_clip5 = 1
        self._options0.pen_clip3 = 1
        self._options._delegate.pen_clip5 = five_prime
        self._options._delegate.pen_clip3 = three_prime


    cdef insert_size_distribution(self,
                                  mean: float,
                                  stddev: float | None = None,
                                  max: float | None = None,
                                  min: float | None = None):
        """bwa mem -I float<,float<,intL<,int>>>"""
        raise NotImplementedError


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
        cdef char *md

        recs_to_return: List[List[AlignedSegment]] = []

        # copy FastqProxy into bwa_seq_t
        num_seqs = len(queries)
        seqs = <kstring_t*>calloc(sizeof(kstring_t), num_seqs)
        for i in range(num_seqs):
            seq = &seqs[i]
            query = queries[i]
            self._copy_seq(queries[i], seq)
            mem_alnregs = mem_align1(opt._delegate, self._index.bwt(), self._index.bns(), self._index.pac(), seq.l, seq.s)
            if opt._delegate.flag & MEM_F_PRIMARY5:
                mem_reorder_primary5(opt._delegate.T, &mem_alnregs)

            # mimic mem_reg2sam from bwamem.c
            recs = []
            XA = NULL
            keep_all = opt._delegate.flag & MEM_F_ALL != 0
            if not keep_all:
                XA = mem_gen_alt(opt._delegate, self._index.bns(), self._index.pac(), &mem_alnregs, seq.l, seq.s)
            num_mem_aln = 0
            for j in range(mem_alnregs.n):
                mem_alnreg = &mem_alnregs.a[j]

                if mem_alnreg.score < opt._delegate.T:
                    continue
                if mem_alnreg.secondary >= 0 and (mem_alnreg.is_alt or not keep_all):
                    continue
                if mem_alnreg.secondary >= 0 and mem_alnreg.secondary < INT_MAX and mem_alnreg.score < mem_alnregs.a[mem_alnreg.secondary].score * opt._delegate.drop_ratio:
                    continue
                mem_aln = mem_reg2aln(opt._delegate, self._index.bns(), self._index.pac(), seq.l, seq.s, mem_alnreg)
                mem_aln.XA = XA[j] if not keep_all else NULL
                if mem_alnreg.secondary >= 0:
                    mem_aln.sub = 1  # don't output sub-optimal score
                if num_mem_aln > 0 and mem_alnreg.secondary < 0:  # if supplementary
                    mem_aln.flag |= 0x10000 if opt._delegate.flag & MEM_F_NO_MULTI else 0x800
                if (opt._delegate.flag & MEM_F_KEEP_SUPP_MAPQ) and num_mem_aln > 0 and mem_alnreg.is_alt > 0 and mem_aln.mapq > mem_alnregs.a[0].mapq:
                    mem_aln.mapq = mem_alnregs.a[0].mapq  # lower
                # create a AlignedSegment record here
                rec = self._unmapped(query=queries[i])

                # set the flags
                mem_aln.flag |= 0x4 if mem_aln.rid < 0 else 0
                mem_aln.flag |= 0x10 if mem_aln.is_rev > 0 else 0
                rec.flag = (mem_aln.flag&0xffff) | (0x100 if mem_aln.flag&0x10000 > 0 else 0)

                # reference id, position, mapq, and cigar
                if rec.is_mapped:
                    rec.reference_id = mem_aln.rid
                    rec.reference_start = mem_aln.pos
                    rec.mapping_quality = mem_aln.mapq
                    cigar = ""
                    for k in range(mem_aln.n_cigar):
                        cigar_opt = mem_aln.cigar[k] & 0xf
                        if opt._delegate.flag & MEM_F_SOFTCLIP == 0 and mem_aln.is_alt == 0 and (cigar_opt == 3 or cigar_opt == 4):
                            cigar_opt = 4 if num_mem_aln > 0 else 3  # // use hard clipping for supplementary alignments
                        cigar = cigar + f"{mem_aln.cigar[k] >> 4}:d" + "MIDS"[cigar_opt]
                    rec.cigarstring = cigar

                # sequence and qualities
                rec.query_sequence = query.sequence if rec.is_forward else query.sequence[::-1]
                if query.quality is not None:
                    rec.query_qualities = query.quality if rec.is_forward else query.quality[::-1]

                # remove leading and trailing soft-clipped bases for non-primary etc.
                if rec.is_mapped and mem_aln.n_cigar > 0 and num_mem_aln > 0 and opt._delegate.flag & MEM_F_SOFTCLIP == 0 and mem_aln.is_alt == 0:
                    qb = 0
                    qe = seq.l
                    if mem_aln.cigar[0] & 0xf == 4 or mem_aln.cigar[0] & 0xf == 3:
                        qb += mem_aln.cigar[0] >> 4
                    if mem_aln.cigar[mem_aln.n_cigar-1] & 0xf == 4 or mem_aln.cigar[mem_aln.n_cigar-1] & 0xf == 3:
                        qe -= mem_aln.cigar[mem_aln.n_cigar-1] >> 4
                    rec.query_sequence = rec.query_sequence [qb:qe]
                    if query.quality is not None:
                        rec.query_qualities = rec.query_qualities[qb:qe]

                if rec.is_mapped:
                    md = <char*>(mem_aln.cigar + mem_aln.n_cigar)
                    attrs = dict()
                    attrs["MD"] = f"{md}"
                    attrs["NM"] = f"{mem_aln.NM}"
                    rec.set_tags(list(attrs.items()))
                    # TODO: other tags: MC, MQ, AS, XS, RG, SA, pa, XA, XB, XR

                num_mem_aln += 1
                recs.append(rec)
            if num_mem_aln == 0:  # unmapped
                recs.append(self._unmapped(query=queries[i]))
            if not keep_all:
                free(XA)
            recs_to_return.append(recs)

        # TODO: free stuff (see mem_reg2sam)
        # TODO: free cigars?
        for i in range(num_seqs):
            free(seqs[i].s)
        free(seqs)

        # TODO: how do we handle retval
        return recs_to_return