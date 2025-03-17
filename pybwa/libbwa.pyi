import enum


# class syntax
@enum.unique
class BwaVerbosity(enum.IntEnum):
    """The verbosity level for the BWA C-API"""

    ERROR = 1
    WARNING = 2
    INFO = 3
    DEBUG = 4

def set_bwa_verbosity(level: BwaVerbosity) -> bool: ...
