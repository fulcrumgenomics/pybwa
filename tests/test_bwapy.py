from pathlib import Path

import pytest
from pysam import FastxRecord

from bwapy import BwaAln
from bwapy import BwaAlnOptions
from bwapy import BwaAlnOptionsBuilder
from bwapy import BwaIndex
from bwapy.libbwamem import BwaMem
from bwapy.libbwamem import BwaMemOptions
from bwapy.libbwamem import BwaMemOptionsBuilder


@pytest.fixture()
def ref_fasta() -> Path:
    cur_dir = Path(__file__).parent
    fasta: Path = cur_dir / "data" / "e_coli_k12.fasta"
    return fasta


def test_bwa_index(ref_fasta: Path) -> None:
    BwaIndex(prefix=ref_fasta)


def test_bwaaln_options() -> None:
    BwaAlnOptions()


def test_bwaalln_options_builder() -> None:
    builder = BwaAlnOptionsBuilder()
    builder.build()
    # TODO: test setting individual options...


@pytest.fixture()
def fastx_record() -> FastxRecord:
    sequence = "gttacctgccgtgagtaaattaaaattttattgacttaggtcactaaatactttaaccaatataggcatagcgcacagac"
    return FastxRecord(name="test", sequence=sequence)


def test_bwaaln(ref_fasta: Path, fastx_record: FastxRecord) -> None:
    opt = BwaAlnOptions()
    bwa = BwaAln(prefix=ref_fasta)

    recs = bwa.align(opt=opt, queries=[fastx_record])
    assert len(recs) == 1
    rec = recs[0]
    assert rec.query_name == "test"
    assert rec.reference_start == 810
    assert rec.cigarstring == "80M"


def test_bwamem_options() -> None:
    BwaMemOptions()


def test_bwamem_options_builder() -> None:
    builder = BwaMemOptionsBuilder()
    builder.build()
    # TODO: test setting individual options...


def test_bwamem(ref_fasta: Path, fastx_record: FastxRecord) -> None:
    opt = BwaMemOptions()
    bwa = BwaMem(prefix=ref_fasta)

    recs = bwa.align(opt=opt, queries=[fastx_record])
    assert len(recs) == 1
    assert len(recs[0]) == 1
    rec = recs[0][0]
    assert rec.query_name == "test"
    assert rec.reference_start == 80
    assert rec.cigarstring == "80M"
