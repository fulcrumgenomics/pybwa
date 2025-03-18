import enum

from pybwa.libbwaaln import set_bwa_aln_verbosity
from pybwa.libbwaindex import set_bwa_idx_verbosity
from pybwa.libbwamem import set_bwa_mem_verbosity

__all__ = [
    "BwaVerbosity",
    "set_bwa_verbosity",
]


# class syntax
@enum.unique
class BwaVerbosity(enum.IntEnum):
    """The verbosity level for the BWA C-API"""

    QUIET = 0
    ERROR = 1
    WARNING = 2
    INFO = 3
    DEBUG = 4


def set_bwa_verbosity(level: BwaVerbosity) -> bool:
    changed = set_bwa_idx_verbosity(level)
    changed |= set_bwa_mem_verbosity(level)
    changed |= set_bwa_aln_verbosity(level)
    return changed
