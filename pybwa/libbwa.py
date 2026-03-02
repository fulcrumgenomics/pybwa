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
    """Suppress all output."""
    ERROR = 1
    """Only output errors."""
    WARNING = 2
    """Output errors and warnings."""
    INFO = 3
    """Output errors, warnings, and informational messages."""
    DEBUG = 4
    """Output all messages including debug details."""


def set_bwa_verbosity(level: BwaVerbosity) -> bool:
    """
    Set the BWA C-API verbosity level.

    By default BWA outputs informational and warning messages to stderr.  Use this to suppress
    or increase the output level.

    Args:
        level: the desired verbosity level

    Returns:
        True if the verbosity level was changed, False otherwise
    """
    changed = _set_bwa_idx_verbosity(level)
    changed |= _set_bwa_mem_verbosity(level)
    changed |= _set_bwa_aln_verbosity(level)
    return changed
