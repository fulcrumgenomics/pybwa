import os
import sysconfig
from pathlib import Path

import pytest

from pybwa import BwaIndex

so_ext = sysconfig.get_config_var('EXT_SUFFIX')
so = os.path.join('tests', '_test_libbwaindex' + so_ext)

try:
    os.unlink("tests/_test_libbwaindex.c")
    os.unlink("tests/_test_libbwaindex.pyxbldc")
    os.unlink(so)
except OSError:
    pass

NO_PYXIMPORT = False
try:
    import pyximport
    import ctypes
    retval = pyximport.install(build_in_temp=False, inplace=True)
    ctypes.cdll.LoadLibrary(so)
    import _test_libbwaindex
except Exception as ex:
    print(f"Cannot import pyximport: {ex}")
    NO_PYXIMPORT = True


@pytest.mark.skipif(NO_PYXIMPORT, reason="no pyximport")
def test_force_bytes_with() -> None:
    _test_libbwaindex.test_force_bytes_with()


def test_bwa_index_build(e_coli_k12_fasta: Path, tmp_path_factory: pytest.TempPathFactory) -> None:
    tmp_dir = Path(str(tmp_path_factory.mktemp("test_bwa_index_build")))
    prefix = tmp_dir / e_coli_k12_fasta.name

    # Build the index
    BwaIndex.index(fasta=e_coli_k12_fasta, prefix=prefix)

    # Check the files exist
    for suffix in [".bwt", ".sa", ".ann", ".pac"]:
        assert prefix.with_suffix(prefix.suffix + suffix).exists()
    assert prefix.with_suffix(".dict").exists()

    # Load it
    BwaIndex(prefix=prefix)


def test_bwa_index_load(e_coli_k12_fasta: Path) -> None:
    BwaIndex(prefix=e_coli_k12_fasta)
