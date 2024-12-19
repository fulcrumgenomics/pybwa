from pathlib import Path
from typing import List

from pysam import AlignedSegment
from pysam import FastxRecord

from pybwa.libbwaindex import BwaIndex

class BwaAlnOptions:
    def __init__(self, max_hits: int = 3) -> None: ...
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

class BwaAln:
    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None) -> None: ...
    def align(self, opt: BwaAlnOptions, queries: List[FastxRecord]) -> List[AlignedSegment]: ...
