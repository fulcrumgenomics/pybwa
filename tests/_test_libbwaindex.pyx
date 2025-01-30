# cython: language_level=3

from pybwa.libbwaindex cimport force_bytes_with
import pytest


def test_force_bytes_with() -> None:
    assert force_bytes_with(None) is None
    assert force_bytes_with(b"foo-bar") == b"foo-bar"
    assert force_bytes_with("foo-bar") == "foo-bar"
    with pytest.raises(TypeError, match="Argument must be a string, bytes or unicode"):
        force_bytes_with(42)
