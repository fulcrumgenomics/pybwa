import enum

from pybwa.libbwaaln import _set_bwa_aln_verbosity
from pybwa.libbwaindex import _set_bwa_idx_verbosity
from pybwa.libbwamem import _set_bwa_mem_verbosity

__all__ = [
    "BwaVerbosity",
    "set_bwa_verbosity",
]


# class syntax
@enum.unique
class BwaVerbosity(enum.IntEnum):
    """The verbosity level for the BWA C-API."""

    QUIET = 0
    ERROR = 1
    WARNING = 2
    INFO = 3
    DEBUG = 4


def set_bwa_verbosity(level: BwaVerbosity) -> bool:
    """Set the BWA C-API verbosity, returning True if changed, false otherwise."""
    changed = _set_bwa_idx_verbosity(level)
    changed |= _set_bwa_mem_verbosity(level)
    changed |= _set_bwa_aln_verbosity(level)
    return changed
