from pathlib import Path
from typing import List

from pysam import AlignedSegment
from pysam import FastxRecord

from pybwa.libbwaalnopt import BwaAlnOptions
from pybwa.libbwaindex import BwaIndex

class BwaAln:
    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None) -> None: ...
    def align(
        self, queries: List[FastxRecord] | List[str], opt: BwaAlnOptions | None = None
    ) -> List[AlignedSegment]: ...
