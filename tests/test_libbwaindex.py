import os
import re
import shutil
import sysconfig
from pathlib import Path
from typing import Union

import pytest

from pybwa import BwaIndex

so_ext = sysconfig.get_config_var("EXT_SUFFIX")
so = os.path.join("tests", "_test_libbwaindex" + so_ext)

try:
    os.unlink("tests/_test_libbwaindex.c")
    os.unlink("tests/_test_libbwaindex.pyxbldc")
    os.unlink(so)
except OSError:
    pass

NO_PYXIMPORT = False
try:
    import ctypes

    import pyximport

    retval = pyximport.install(build_in_temp=False, inplace=True)
    ctypes.cdll.LoadLibrary(so)
    import _test_libbwaindex
except Exception as ex:
    print(f"Cannot import pyximport: {ex}")
    NO_PYXIMPORT = True


@pytest.mark.skipif(NO_PYXIMPORT, reason="no pyximport")
def test_force_bytes_with() -> None:
    _test_libbwaindex.test_force_bytes_with()


@pytest.mark.parametrize("in_place", [True, False])
@pytest.mark.parametrize("prefix_as_str", [True, False])
def test_bwa_index_build(
    e_coli_k12_fasta: Path,
    tmp_path_factory: pytest.TempPathFactory,
    in_place: bool,
    prefix_as_str: bool,
) -> None:
    src_dir = Path(str(tmp_path_factory.mktemp("test_bwa_index_build_src")))
    fasta = src_dir / e_coli_k12_fasta.name
    shutil.copy(f"{e_coli_k12_fasta}", f"{fasta}")

    if in_place:
        prefix = fasta
    else:
        dest_dir = Path(str(tmp_path_factory.mktemp("test_bwa_index_build_test")))
        prefix = dest_dir / e_coli_k12_fasta.name

    # Build the index
    fasta_arg = f"{fasta}" if prefix_as_str else fasta
    BwaIndex.index(fasta=fasta_arg, prefix=None if in_place else prefix)

    # Check the files exist
    for suffix in [".bwt", ".sa", ".ann", ".pac"]:
        assert prefix.with_suffix(prefix.suffix + suffix).exists()
    assert prefix.with_suffix(".dict").exists()

    # Load it
    prefix_arg = f"{prefix}" if prefix_as_str else prefix
    BwaIndex(prefix=prefix_arg)


@pytest.mark.parametrize("prefix_as_str", [True, False])
def test_bwa_index_load(
    e_coli_k12_fasta: Path, tmp_path_factory: pytest.TempPathFactory, prefix_as_str: bool
) -> None:
    prefix = f"{e_coli_k12_fasta}" if prefix_as_str else e_coli_k12_fasta
    BwaIndex(prefix=prefix)


def test_bwa_index_load_path_does_not_exist(
    e_coli_k12_fasta: Path, tmp_path_factory: pytest.TempPathFactory
) -> None:
    # does not load when the prefix does not exist
    with pytest.raises(FileNotFoundError, match=r"could not locate the index file \[prefix\]:.*"):
        BwaIndex(prefix="/path/does/not/exist")


@pytest.mark.parametrize(
    "suffix, append, message",
    [
        (".bwt", True, "prefix"),
        (".sa", True, None),
        (".ann", True, None),
        (".amb", True, None),
        (".pac", True, None),
        (".dict", False, "use `samtools dict`"),
    ],
)
def test_bwa_index_fail(
    e_coli_k12_fasta: Path,
    tmp_path_factory: pytest.TempPathFactory,
    suffix: str,
    append: bool,
    message: Union[str, None],
) -> None:
    src_dir = Path(str(tmp_path_factory.mktemp(f"test_bwa_index_build_src{suffix}")))
    shutil.copytree(e_coli_k12_fasta.parent, src_dir, dirs_exist_ok=True)
    fasta = src_dir / e_coli_k12_fasta.name
    assert fasta.exists()

    # delete the file with the given suffix
    if append:
        path = fasta.with_suffix(fasta.suffix + suffix)
    else:
        path = fasta.with_suffix(suffix)
    path.unlink()

    pattern = "could not locate .* \\["
    if message is None:
        pattern += suffix
    else:
        pattern += message
    pattern += "\\]:.*"

    # now it should fail loading the file because the suffix is missing
    with pytest.raises(FileNotFoundError, match=re.compile(pattern)):
        assert fasta.exists()
        BwaIndex(prefix=fasta)
