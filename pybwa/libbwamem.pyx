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
from pybwa.libbwamemopt cimport BwaMemOptions

__all__ = [
    "BwaMem",
]

cdef int _BWA_MEM_TO_PYSAM_CIGAR_OPERATOR[5]
_BWA_MEM_TO_PYSAM_CIGAR_OPERATOR[0] = CMATCH
_BWA_MEM_TO_PYSAM_CIGAR_OPERATOR[1] = CINS
_BWA_MEM_TO_PYSAM_CIGAR_OPERATOR[2] = CDEL
_BWA_MEM_TO_PYSAM_CIGAR_OPERATOR[3] = CSOFT_CLIP
_BWA_MEM_TO_PYSAM_CIGAR_OPERATOR[4] = CHARD_CLIP

cdef inline int _to_pysam_cigar_op(int x):
    return _BWA_MEM_TO_PYSAM_CIGAR_OPERATOR[x]


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
                cigartuples = []
                cigar_len_sum = 0
                for k in range(mem_aln.n_cigar):
                    cigar_op = mem_aln.cigar[k] & 0xf
                    if not opt.soft_clip_supplementary and mem_aln.is_alt == 0 and (
                            cigar_op == 3 or cigar_op == 4):
                        cigar_op = 4 if j > 0 else 3  # // use hard clipping for supplementary alignments
                    cigar_len = mem_aln.cigar[k] >> 4
                    cigartuples.append((_to_pysam_cigar_op(cigar_op), cigar_len))
                    if cigar_op < 4:
                        cigar_len_sum += cigar_len
                rec.cigartuples = cigartuples

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
