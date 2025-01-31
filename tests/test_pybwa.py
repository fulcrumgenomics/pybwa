import sysconfig
from pathlib import Path

import pybwa


def test_get_includes() -> None:
    names = [Path(p).name for p in pybwa._get_include()]
    assert names == ["pybwa", "bwa"]


def test_get_defines() -> None:
    assert pybwa._get_defines() == []


def test_get_libraries() -> None:
    so = sysconfig.get_config_var("EXT_SUFFIX")
    names = [Path(p).name for p in pybwa._get_libraries()]
    assert names == [
        "libbwaaln" + so,
        "libbwaindex" + so,
        "libbwamem" + so,
    ]
