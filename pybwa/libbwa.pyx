# cython: language_level=3
import enum


__all__ = [
    "BwaVerbosity",
    "set_bwa_verbosity",
]


# class syntax
@enum.unique
class BwaVerbosity(enum.IntEnum):
    """The verbosity level for the BWA C-API"""

    ERROR = 1
    WARNING = 2
    INFO = 3
    DEBUG = 4



def set_bwa_verbosity(level: BwaVerbosity) -> bool:
    """Set the BWA C-API verbosity, returning True if changed, false otherwise."""
    return 1 == set_bwa_c_verbosity(int(level))
