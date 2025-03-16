# cython: language_level=3

from pathlib import Path
from typing import List
from fgpyo.sequence import reverse_complement

from libc.stdint cimport uint8_t
from libc.stdlib cimport calloc, free
from libc.string cimport strncpy
from pysam import FastxRecord, AlignedSegment, qualitystring_to_array, CMATCH, CINS, CDEL, CSOFT_CLIP
from pybwa.libbwaindex cimport force_bytes
from pybwa.libbwaindex cimport BwaIndex
from pybwa.libbwaalnopt cimport BwaAlnOptions


__all__ = [
    "BwaAln",
]

cdef int _BWA_ALN_TO_PYSAM_CIGAR_OPERATOR[4]
_BWA_ALN_TO_PYSAM_CIGAR_OPERATOR[FROM_M] = CMATCH
_BWA_ALN_TO_PYSAM_CIGAR_OPERATOR[FROM_I] = CINS
_BWA_ALN_TO_PYSAM_CIGAR_OPERATOR[FROM_D] = CDEL
_BWA_ALN_TO_PYSAM_CIGAR_OPERATOR[FROM_S] = CSOFT_CLIP

cdef inline int _to_pysam_cigar_op(int x):
    return _BWA_ALN_TO_PYSAM_CIGAR_OPERATOR[x]


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
    
    @staticmethod
    def __to_str(_bytes: bytes) -> str:
        return _bytes.decode('utf-8')

    cdef _copy_seq(self, q: FastxRecord, bwa_seq_t* s, is_comp: bool):
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
            s.rseq[i] = nst_nt4_table[ord(base)]
            s.qual[i] = 40  # FIXME: parameterize
        s.seq[seq_len] = b'\0'
        s.qual[seq_len] = b'\0'
        seq_reverse(seq_len, s.seq,
                    0)  #  // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
        seq_reverse(seq_len, s.rseq, 1 if is_comp else 0)

        s.name = <char *> calloc(sizeof(char), len(q.name) + 1)
        strncpy(s.name, force_bytes(q.name), len(q.name))
        s.name[len(q.name)] = b'\0'

    cdef _build_alignment(self, query: FastxRecord, bwa_seq_t *seq, opt: BwaAlnOptions, kstring_t *kstr):
        cdef int reference_id
        cdef int nm
        cdef char XT
        cdef bwt_multi1_t* hit

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

        if seq.type == BWA_TYPE_NO_MATCH:  # unmapped read
            # TODO: custom bwa tags: RG, BC, XC
            return rec

        # get the span of reference bases
        ref_len_in_alignment = pos_end(seq) - seq.pos

        # reference id
        nn = bns_cnt_ambi(self._index.bns(), seq.pos, ref_len_in_alignment, &reference_id)
        rec.reference_id = reference_id

        # make this unmapped if we map off the end of the contig
        rec.is_unmapped = seq.pos + ref_len_in_alignment - self._index.bns().anns[
            rec.reference_id].offset > self._index.bns().anns[rec.reference_id].len

        # strand, reference start, and mapping quality
        if seq.strand:
            rec.is_reverse = True
            rec.query_sequence = reverse_complement(query.sequence)
            if query.quality is not None:
                rec.query_qualities = qualitystring_to_array(query.quality[::-1])
        rec.reference_start = seq.pos - self._index.bns().anns[rec.reference_id].offset  # 0-based
        rec.mapping_quality = seq.mapQ

        # cigar
        cigartuples = []
        if seq.cigar:
            for j in range(seq.n_cigar):
                cigar_len = __cigar_len(seq.cigar[j])
                cigar_op = __cigar_op(seq.cigar[j])
                cigartuples.append((_to_pysam_cigar_op(cigar_op), cigar_len))
        elif seq.type != BWA_TYPE_NO_MATCH:
            cigartuples.append((CMATCH, seq.len))
        rec.cigartuples = cigartuples

        # # tags
        if seq.type != BWA_TYPE_NO_MATCH:
            attrs = dict()
            tag = "NM" if ((opt._delegate.mode & BWA_MODE_COMPREAD) != 0) else "CM"
            attrs[tag] = int(seq.nm)
            if nn > 0:
                attrs["XN"] = nn
            if seq.type != BWA_TYPE_MATESW:  # X0 and X1 are not available for this type of alignment
                attrs["X0"] = seq.c1
                if seq.c1 <= opt.stop_at_max_best_hits:
                    attrs["X1"] = seq.c2
            # SM, AM, X0, and X1 are for paired end reads, so omitted
            attrs["XM"] = seq.n_mm
            attrs["XO"] = seq.n_gapo
            attrs["XG"] = seq.n_gapo + seq.n_gape
            if seq.md != NULL:
                attrs["MD"] = self.__to_str(seq.md)
            # print multiple hits
            if seq.n_multi > 0:
                XA = ""
                for j in range(seq.n_multi):
                    hit = &seq.multi[j]
                    end = pos_end_multi(hit, seq.len - hit.pos)
                    nn = bns_cnt_ambi(self._index.bns(), hit.pos, end, &reference_id)
                    XA += self.__to_str(self._index.bns().anns[reference_id].name)
                    XA += "," + ('-' if hit.strand != 0 else '+')
                    XA += str(hit.pos - self._index.bns().anns[reference_id].offset + 1)
                    XA += ","
                    if hit.cigar == NULL:
                        XA += f"{seq.len}M"
                    else:
                        for k in range(hit.n_cigar):
                            cigar_len = __cigar_len(hit.cigar[k])
                            cigar_op = "MIDS"[__cigar_op(hit.cigar[k])]
                            XA += f"{cigar_len}{cigar_op}"
                    XA += f",{hit.gap + hit.mm};"
                attrs["XA"] = XA
            attrs["HN"] = seq.n_occ
            # NB: the type 'A' (printable character) of the XA tag must be explicitly given
            # below, therefore we extract the attrs dictionary as a list of tuples first.
            tags = [('XT', b'N' if nn > 10 else "NURM"[seq.type], 'A')]
            tags.extend(list(attrs.items()))
            rec.set_tags(tags)

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
            self._copy_seq(queries[i], &seqs[i], (opt._delegate.mode & BWA_MODE_COMPREAD) != 0)
            seqs[i].tid = -1

        # this is `bwa aln`, and the rest is `bwa samse`
        bwa_cal_sa_reg_gap_threaded(0, self._index.bwt(), num_seqs, seqs, gap_opt)

        # create the full alignment
        for i in range(num_seqs):
            s = &seqs[i]
            # bwa_cal_sa_reg_gap frees name, seq, rseq, and qual, so add them back in again
            self._copy_seq(queries[i], s, (opt._delegate.mode & BWA_MODE_COMPREAD) != 0)
            bwa_aln2seq_core(s.n_aln, s.aln, s, 1, opt.max_hits)

        # # calculate the genomic position given the suffix array offsite
        bwa_cal_pac_pos_with_bwt(self._index.bns(), num_seqs, seqs, gap_opt.max_diff, gap_opt.fnr, self._index.bwt())

        # refine gapped alignment
        bwa_refine_gapped(self._index.bns(), num_seqs, seqs, self._index.pac(), opt.with_md)

        # create the AlignedSegment from FastxRecord and bwa_seq_t.
        recs = [
            self._build_alignment(query=queries[i], seq=&seqs[i], opt=opt, kstr=kstr)
            for i in range(num_seqs)
        ]

        bwa_free_read_seq(num_seqs, seqs)
        free(kstr.s)
        free(kstr)

        return recs
