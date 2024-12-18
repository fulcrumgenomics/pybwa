from pathlib import Path
from typing import List

from pysam import AlignedSegment
from pysam import FastxRecord

from bwapy.libbwaindex import BwaIndex

class BwaMemOptions:
    _ignore_alt: bool
    def __init__(self) -> None: ...

class BwMemOptionsBuilder:
    _options: BwaMemOptions
    _options0: BwaMemOptions
    def __init__(self, options: BwaMemOptions | None = None) -> None: ...
    def build(self) -> BwaMemOptions: ...

class BwaMem:
    _index: BwaIndex
    def __init__(self, prefix: str | Path | None = None, index: BwaIndex | None = None) -> None: ...
    def align(
        self, opt: BwaMemOptions, queries: List[FastxRecord]
    ) -> List[List[AlignedSegment]]: ...
