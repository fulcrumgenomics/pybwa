from pathlib import Path
from typing import List

from pysam import AlignedSegment
from pysam import AlignmentHeader
from pysam import FastxRecord

class BwaAlnOptions:
    max_hits: int
    def __init__(self, max_hits: int = 3) -> None: ...

class BwaAlnOptionsBuilder:
    _options: BwaAlnOptions
    def __init__(self, options: BwaAlnOptions | None = None) -> None: ...
    def build(self) -> BwaAlnOptions: ...
    def max_mismatches(self, value: int) -> BwaAlnOptionsBuilder: ...  # -n <int>
    # def fnr(self, value: float) -> BwaOptionsBuilder: ... # -n <float>
    def max_gap_opens(self, value: int) -> BwaAlnOptionsBuilder: ...  # -o <int>
    def max_gap_extensions(self, value: int) -> BwaAlnOptionsBuilder: ...  # -e <int>
    def min_indel_to_end_distance(self, value: int) -> BwaAlnOptionsBuilder: ...  # -i <int>
    def max_occurences_for_extending_long_deletion(
        self, value: int
    ) -> BwaAlnOptionsBuilder: ...  # -d <int>
    def seed_length(self, value: int) -> BwaAlnOptionsBuilder: ...  # -l <int>
    def max_mismatches_in_seed(self, value: int) -> BwaAlnOptionsBuilder: ...  # -k <int>
    def mismatch_penalty(self, value: int) -> BwaAlnOptionsBuilder: ...  # -M <int>
    def gap_open_penalty(self, value: int) -> BwaAlnOptionsBuilder: ...  # -O <int>
    def gap_extension_penalty(self, value: int) -> BwaAlnOptionsBuilder: ...  # -E <int>
    def stop_at_max_best_hits(self, value: int) -> BwaAlnOptionsBuilder: ...  # -R <int>
    def max_hits(self, value: int) -> BwaAlnOptionsBuilder: ...  # bwa samse -n <int>
    def log_scaled_gap_penalty(self, value: bool = True) -> BwaAlnOptionsBuilder: ...  # -L

ERROR_HANDLER: str
TEXT_ENCODING: str

class BwaIndex:
    header: AlignmentHeader
    def __init__(self, prefix: str | Path) -> None: ...

class BwaAln:
    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None) -> None: ...
    def align(self, opt: BwaAlnOptions, queries: List[FastxRecord]) -> List[AlignedSegment]: ...
