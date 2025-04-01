import sysconfig
from pathlib import Path

import pytest
from pysam.libcfaidx import FastxRecord

import pybwa
from pybwa import BwaAln
from pybwa import BwaMem
from pybwa import BwaMemOptions
from pybwa import BwaVerbosity
from pybwa import set_bwa_verbosity
from pybwa.libbwaindex import BwaIndex


def test_get_includes() -> None:
    names = [Path(p).name for p in pybwa._get_include()]
    assert names == ["pybwa", "htslib", "bwa"]


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
    set_bwa_verbosity(BwaVerbosity.INFO)


@pytest.mark.parametrize("level", list(BwaVerbosity))
def test_bwa_verbosity_bwa_index(
    level: BwaVerbosity, e_coli_k12_fasta: Path, capfd: pytest.CaptureFixture[str]
) -> None:
    set_bwa_verbosity(level)
    BwaIndex(e_coli_k12_fasta)
    captured = capfd.readouterr()
    assert captured.out == ""
    if level >= BwaVerbosity.INFO:
        assert captured.err == "[M::bwa_idx_load_from_disk] read 0 ALT contigs\n"
    else:
        assert captured.err == ""
    set_bwa_verbosity(BwaVerbosity.INFO)


@pytest.mark.parametrize("level", list(BwaVerbosity))
def test_bwa_verbosity_bwa_aln(
    level: BwaVerbosity,
    e_coli_k12_fasta: Path,
    e_coli_k12_fastx_record: FastxRecord,
    capfd: pytest.CaptureFixture[str],
) -> None:
    set_bwa_verbosity(level)
    bwa = BwaAln(prefix=e_coli_k12_fasta)
    assert len(bwa.align(queries=[e_coli_k12_fastx_record])) == 1
    captured = capfd.readouterr()
    assert captured.out == ""
    if level >= BwaVerbosity.INFO:
        assert captured.err == "[M::bwa_idx_load_from_disk] read 0 ALT contigs\n"
    else:
        assert captured.err == ""
    set_bwa_verbosity(BwaVerbosity.INFO)


_BWA_MEM_LOG: str = """[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::mem_process_seqs] Processed 1 reads in 0.000 CPU sec, 0.000 real sec\n"""
# The expected standard error for the `test_bwa_verbosity_bwa_mem` test


@pytest.mark.parametrize("level", list(BwaVerbosity))
def test_bwa_verbosity_bwa_mem(
    level: BwaVerbosity,
    e_coli_k12_fasta: Path,
    e_coli_k12_fastx_record: FastxRecord,
    capfd: pytest.CaptureFixture[str],
) -> None:
    set_bwa_verbosity(level)
    opt = BwaMemOptions(with_xr_tag=True)
    bwa = BwaMem(prefix=e_coli_k12_fasta)

    assert len(bwa.align(opt=opt, queries=[e_coli_k12_fastx_record])) == 1
    captured = capfd.readouterr()
    assert captured.out == ""
    if level >= BwaVerbosity.INFO:
        assert captured.err == _BWA_MEM_LOG
    else:
        assert captured.err == ""
    set_bwa_verbosity(BwaVerbosity.INFO)
