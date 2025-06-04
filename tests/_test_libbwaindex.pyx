# cython: language_level=3

from pybwa.libbwaindex cimport force_bytes
import pytest


def test_force_bytes() -> None:
    assert force_bytes(None) is None
    assert force_bytes(b"foo-bar") == b"foo-bar"
    assert force_bytes("foo-bar") == b"foo-bar"
    with pytest.raises(TypeError, match="Argument must be a string, bytes or unicode"):
        force_bytes(42)
