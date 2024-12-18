


class BwaMemOptions:
    _ignore_alt: bool
    def __init__(self) -> None: ...


class BwMemOptionsBuilder:
    _options: BwaMemOptions
    def __init__(self, options: BwaMemOptions | None = None) -> None: ...
    def build(self) -> BwaMemOptions: ...