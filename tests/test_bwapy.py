from pathlib import Path

import pytest
from pysam import FastxRecord

from bwapy import BwaAln
from bwapy import BwaAlnOptions
from bwapy import BwaAlnOptionsBuilder
from bwapy import BwaIndex


@pytest.fixture()
def ref_fasta() -> Path:
    cur_dir = Path(__file__).parent
    fasta: Path = cur_dir / "data" / "e_coli_k12.fasta"
    return fasta


def test_bwapy_options() -> None:
    BwaAlnOptions()


def test_bwapy_options_builder() -> None:
    builder = BwaAlnOptionsBuilder()
    builder.build()
    # TODO: test setting individual options...


def test_bwapy_index(ref_fasta: Path) -> None:
    BwaIndex(prefix=ref_fasta)


def test_bwapy(ref_fasta: Path) -> None:
    opt = BwaAlnOptions()
    bwa = BwaAln(prefix=ref_fasta)
    sequence = "gttacctgccgtgagtaaattaaaattttattgacttaggtcactaaatactttaaccaatataggcatagcgcacagac"
    fastqs = [FastxRecord(name="test", sequence=sequence)]

    recs = bwa.align(opt=opt, queries=fastqs)
    assert len(recs) == 1
    rec = recs[0]
    assert rec.query_name == "test"
    assert rec.reference_start == 81
    assert rec.cigarstring == "80M"
