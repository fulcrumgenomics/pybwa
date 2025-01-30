"""
Fixtures intended to be shared across multiple files in the tests directory.
"""

from pathlib import Path

import pytest
from pysam import FastxRecord


@pytest.fixture(scope="function")
def temp_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("test_vcf")


@pytest.fixture(scope="session")
def e_coli_k12_fasta() -> Path:
    """The path to the e. Coli K12 reference FASTA"""
    cur_dir = Path(__file__).parent
    fasta: Path = cur_dir / "data" / "e_coli_k12.fasta"
    return fasta


@pytest.fixture(scope="session")
def e_coli_k12_fastx_record() -> FastxRecord:
    """Sequence-only FastxRecord that maps to the forward strand at positino 80 for 80bp"""
    sequence = "gttacctgccgtgagtaaattaaaattttattgacttaggtcactaaatactttaaccaatataggcatagcgcacagac"
    return FastxRecord(name="test", sequence=sequence)
