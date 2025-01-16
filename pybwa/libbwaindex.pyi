import enum
from pathlib import Path

from pysam import AlignmentHeader

ERROR_HANDLER: str
TEXT_ENCODING: str

@enum.unique
class BwaIndexBuildMethod(enum.Enum):
    """The BWT construction algorithm (:code:`bwa index -a <str>`)"""

    AUTO = enum.auto()
    """Choose the algorithm automatically (based on the reference length)"""
    RB2 = enum.auto()
    """rb2 (does not work for "long" genomes)"""
    BWTSW = enum.auto()
    """bwtsw (does not work for "short" genomes)"""
    IS = enum.auto()
    """is (does not work for "long" genomes)"""

class BwaIndex:
    header: AlignmentHeader
    def __init__(
        self, prefix: str | Path, bwt: bool = ..., bns: bool = ..., pac: bool = ...
    ) -> None: ...
    @classmethod
    def index(
        cls,
        fasta: str | Path,
        method: BwaIndexBuildMethod = BwaIndexBuildMethod.AUTO,
        prefix: str | Path | None = None,
        block_size: int = 10000000,
        out_64: bool = False,
    ) -> None: ...
