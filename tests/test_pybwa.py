import sysconfig
from pathlib import Path

import pybwa
from pybwa.libbwa import BwaVerbosity


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


def test_bwa_verbosity() -> None:
    # Tests that we _start_ with INFO verbosity
    assert not pybwa.set_bwa_verbosity(BwaVerbosity.INFO)
    for level in BwaVerbosity:
        # Change the level
        assert pybwa.set_bwa_verbosity(level)
        # The level remains the same
        assert not pybwa.set_bwa_verbosity(level)
    # Change it back
    assert pybwa.set_bwa_verbosity(BwaVerbosity.INFO)
