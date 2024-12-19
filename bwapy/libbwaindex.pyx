# cython: language_level=3

from pathlib import Path

from cpython cimport PyBytes_Check, PyUnicode_Check
from pysam import AlignmentFile

__all__ = [
    "BwaIndex",
]


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
    """Contains the index and nucleotide sequence for Bwa.

    Args:
        prefix (str | Path): the path prefix for teh BWA index
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

    cdef _load_index(self, prefix, mode):
        prefix = bwa_idx_infer_prefix(force_bytes(prefix))
        if not prefix:
            # FIXME: better error message
            raise Exception("Could not find the index")

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
