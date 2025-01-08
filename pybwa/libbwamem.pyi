import enum
from pathlib import Path
from typing import List

from pysam import AlignedSegment
from pysam import FastxRecord

from pybwa.libbwaindex import BwaIndex

# class syntax
@enum.unique
class BwaMemMode(enum.Enum):
    """The read type for overriding multiple options"""

    PACBIO = enum.auto()
    ONT2D = enum.auto()
    INTRACTG = enum.auto()

class BwaMemOptions:
    def __init__(self, finalize: bool = False) -> None: ...
    _ignore_alt: bool
    min_seed_len: int
    mode: BwaMemMode
    band_width: int
    match_score: int
    mismatch_penalty: int
    minimum_score: int
    unpaired_penalty: int
    n_threads: int
    skip_pairing: bool
    output_all_for_fragments: bool
    interleaved_paired_end: bool
    short_split_as_secondary: bool
    skip_mate_rescue: bool
    soft_clip_supplementary: bool
    with_xr_tag: bool
    query_coord_as_primary: bool
    keep_mapq_for_supplementary: bool
    with_xb_tag: bool
    max_occurrences: int
    off_diagonal_x_dropoff: float
    ignore_alternate_contigs: bool
    internal_seed_split_factor: float
    drop_chain_fraction: float
    max_mate_rescue_rounds: int
    min_seeded_bases_in_chain: int
    seed_occurrence_in_3rd_round: int
    xa_max_hits: int | tuple[int, int]
    xa_drop_ratio: float
    gap_open_penalty: int | tuple[int, int]
    gap_extension_penalty: int | tuple[int, int]
    clipping_penalty: int | tuple[int, int]

class BwaMemOptionsBuilder(BwaMemOptions):
    def __init__(self, options: BwaMemOptions | None = None) -> None: ...
    def build(self) -> BwaMemOptions: ...

class BwaMem:
    _index: BwaIndex
    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None) -> None: ...
    def align(
        self, opt: BwaMemOptions, queries: List[FastxRecord]
    ) -> List[List[AlignedSegment]]: ...
