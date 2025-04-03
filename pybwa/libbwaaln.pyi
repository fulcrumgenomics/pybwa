from pathlib import Path
from typing import List

from pysam import AlignedSegment
from pysam import FastxRecord

from pybwa.libbwaindex import BwaIndex

class BwaAlnOptions:
    def __init__(
        self,
        max_mismatches: int = -1,
        max_gap_opens: int = 1,
        max_gap_extensions: int = 6,
        min_indel_to_end_distance: int = 5,
        max_occurrences_for_extending_long_deletion: int = 10,
        seed_length: int = 32,
        max_mismatches_in_seed: int = 2,
        mismatch_penalty: int = 3,
        gap_open_penalty: int = 11,
        gap_extension_penalty: int = 4,
        stop_at_max_best_hits: int = 30,
        max_hits: int = 3,
        log_scaled_gap_penalty: bool = False,
        find_all_hits: bool = False,
        with_md: bool = False,
        threads: int = 1,
    ) -> None: ...
    max_mismatches: int  # -n <int>
    # fnr:float # -n <float>
    max_gap_opens: int  # -o <int>
    max_gap_extensions: int  # -e <int>
    min_indel_to_end_distance: int  # -i <int>
    max_occurrences_for_extending_long_deletion: int  # -d <int>
    seed_length: int  # -l <int>
    max_mismatches_in_seed: int  # -k <int>
    mismatch_penalty: int  # -M <int>
    gap_open_penalty: int  # -O <int>
    gap_extension_penalty: int  # -E <int>
    stop_at_max_best_hits: int  # -R <int>
    max_hits: int  # bwa samse -n <int>
    log_scaled_gap_penalty: bool = True  # -L
    find_all_hits: bool = False  # -N
    with_md: bool = True  # bwa samse -d
    threads: int  # -t <int>

class BwaAln:
    _max_hits: int
    _with_md: bool
    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None) -> None: ...
    def align(
        self, queries: List[FastxRecord] | List[str], opt: BwaAlnOptions | None = None
    ) -> List[AlignedSegment]: ...
    def reinitialize_seed(self) -> None: ...

def _set_bwa_aln_verbosity(level: int) -> bool: ...
def _assert_gap_opt_are_the_same() -> None: ...
