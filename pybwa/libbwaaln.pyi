from pathlib import Path
from typing import List

from pysam import AlignedSegment
from pysam import FastxRecord

from pybwa.libbwaindex import BwaIndex

class BwaAlnOptions:
    def __init__(
        self,
        max_mismatches: int | None = None,
        max_gap_opens: int | None = None,
        max_gap_extensions: int | None = None,
        min_indel_to_end_distance: int | None = None,
        max_occurrences_for_extending_long_deletion: int | None = None,
        seed_length: int | None = None,
        max_mismatches_in_seed: int | None = None,
        mismatch_penalty: int | None = None,
        gap_open_penalty: int | None = None,
        gap_extension_penalty: int | None = None,
        stop_at_max_best_hits: int | None = None,
        max_hits: int | None = 3,
        log_scaled_gap_penalty: bool | None = None,
        find_all_hits: bool | None = None,
        with_md: bool | None = None,
        threads: int | None = None,
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
    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None) -> None: ...
    def align(
        self, queries: List[FastxRecord] | List[str], opt: BwaAlnOptions | None = None
    ) -> List[AlignedSegment]: ...
