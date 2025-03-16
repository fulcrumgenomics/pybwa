from pathlib import Path
from typing import List

from pysam import AlignedSegment
from pysam import FastxRecord

from pybwa.libbwaindex import BwaIndex
from pybwa.libbwamemopt import BwaMemOptions


class BwaMem:
    _index: BwaIndex
    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None) -> None: ...
    def align(
        self, queries: List[FastxRecord] | List[str], opt: BwaMemOptions | None = None
    ) -> List[List[AlignedSegment]]: ...
