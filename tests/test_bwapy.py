from pathlib import Path

import pytest
from pysam import FastxRecord

from bwapy import BwaAln
from bwapy import BwaAlnOptions
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
    assert rec.reference_start == 80
    assert rec.cigarstring == "80M"


def test_bwamem_options() -> None:
    # default options
    options = BwaMemOptions()
    assert options.min_seed_len == 19
    # build with default options
    builder = BwaMemOptionsBuilder()
    options = builder.build()
    assert options.min_seed_len == 19
    # build with custom option
    builder = BwaMemOptionsBuilder()
    builder.min_seed_len = 20
    options = builder.build()
    assert options.min_seed_len == 20


def test_bwamem(ref_fasta: Path, fastx_record: FastxRecord) -> None:
    opt = BwaMemOptions()
    bwa = BwaMem(prefix=ref_fasta)

    recs = bwa.align(opt=opt, queries=[fastx_record])
    assert len(recs) == 1
    assert len(recs[0]) == 1
    rec = recs[0][0]
    assert rec.query_name == "test"
    assert rec.reference_start == 80
    assert rec.cigarstring == "80M", print(str(rec))
    # TODO: test multi-mapping etc
