# cython: language_level=3

from pathlib import Path

from cpython cimport PyBytes_Check, PyUnicode_Check
from pysam import AlignmentFile, FastxFile
import enum


__all__ = [
    "BwaIndex",
    "BwaIndexBuildMethod",
]

@enum.unique
class BwaIndexBuildMethod(enum.Enum):
    """The BWT construction algorithm (:code:`bwa index -a <str>`)"""

    AUTO = enum.auto()
    RB2 = enum.auto()
    BWTSW = enum.auto()
    IS = enum.auto()


cdef str ERROR_HANDLER = 'strict'
cdef str TEXT_ENCODING = 'utf-8'


cdef bytes force_bytes(object s):
    return force_bytes_with(s, None, None)


cdef bytes force_bytes_with(object s, encoding: str | None = None, errors: str | None = None):
    """convert string or unicode object to bytes, assuming
    utf8 encoding.
    """
    if s is None:
        return None
    elif PyBytes_Check(s):
        return s
    elif PyUnicode_Check(s):
        return s.encode(encoding or TEXT_ENCODING, errors or ERROR_HANDLER)
    else:
        raise TypeError("Argument must be string, bytes or unicode.")



cdef class BwaIndex:
    """Contains the index and nucleotide sequence for Bwa.  Use :code:`bwa index` on the command
    line to generate the bwa index.

    Note: the accompanying sequence dictionary must exist (i.e. `.dict` file, which can be generated
    with :code:`samtools dict <fasta>`).

    Args:
        prefix (str | Path): the path prefix for the BWA index (typically a FASTA)
        bwt (bool): load the BWT (FM-index)
        bns (bool): load the BNS (reference sequence metadata)
        pac (bool): load the PAC (the actual 2-bit encoded reference sequences with 'N' converted to a
             random base)
    """

    def __init__(self, prefix: str | Path, bwt: bool = True, bns: bool = True, pac: bool = True) -> None:
        cdef int mode

        mode = 0
        if bwt:
            mode |= BWA_IDX_BWT
        if bns:
            mode |= BWA_IDX_BNS
        if pac:
            mode |= BWA_IDX_PAC
        self._load_index(f"{prefix}", mode)

    @classmethod
    def index(cls,
              fasta: str | Path,
              method: BwaIndexBuildMethod= BwaIndexBuildMethod.AUTO,
              prefix: str | Path | None = None,
              block_size: int = 10000000,
              out_64: bool = False) -> None:
        """Indexes a given FASTA.  Also builds the sequence dictionary (.dict).

        Args:
            fasta (str | Path): the path to the FASTA to index
            method (BwaIndexBuildMethod): the BWT construction algorithm (:code:`bwa index -a <str>`)
            prefix (str | Path | None): the path prefix for the BWA index (typically a FASTA)
            block_size (int): block size for the bwtsw algorithm (effective with -a bwtsw)
            out_64 (bool): index files named as :code:`<in.fasta>.64.*` instead of :code:`<in.fasta>.*`
        """
        if prefix is None:
            prefix = fasta

        # Build the BWA index
        bwa_idx_build(force_bytes(f"{fasta}"), force_bytes(f"{prefix}"), method.value, block_size)

        # Build the sequence dictionary
        dict_fn = Path(prefix.with_suffix(".dict"))
        with FastxFile(filename=f"{fasta}") as reader, dict_fn.open("w") as writer:
            writer.write("@HD\tVN:1.5\tSO:unsorted\n")
            for rec in reader:
                writer.write(f"@SQ\tSN:{rec.name}\tLN:{len(rec.sequence)}\n")


    cdef _load_index(self, prefix, mode):
        prefix = bwa_idx_infer_prefix(force_bytes(prefix))
        if not prefix:
            raise Exception(f"Could not find the index at: {prefix}")

        self._delegate = bwa_idx_load(prefix, mode)

        # the SAM header from the sequence dictionary
        seq_dict = Path(prefix.decode("utf-8")).with_suffix(".dict")
        # TODO: error message when seq_dict is missing?
        with seq_dict.open("r") as fh:
            with AlignmentFile(fh) as reader:
                self.header = reader.header

    cdef bwt_t *bwt(self):
        return self._delegate.bwt

    cdef bntseq_t *bns(self):
        return self._delegate.bns

    cdef uint8_t *pac(self):
        return self._delegate.pac

    def __dealloc__(self):
        bwa_idx_destroy(self._delegate)
        self._delegate = NULL
