import enum
from pathlib import Path

from pysam import AlignmentHeader

ERROR_HANDLER: str
TEXT_ENCODING: str

@enum.unique
class BwaIndexBuildMethod(enum.Enum):
    """The BWT construction algorithm (:code:`bwa index -a <str>`)"""

    AUTO = enum.auto()
    RB2 = enum.auto()
    BWTSW = enum.auto()
    IS = enum.auto()

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
