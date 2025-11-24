"""Fixtures intended to be shared across multiple files in the tests directory."""

from pathlib import Path

import pytest
from pysam import FastxRecord


@pytest.fixture(scope="function")
def temp_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("test_vcf")


@pytest.fixture(scope="session")
def e_coli_k12_fasta() -> Path:
    """The path to the e. Coli K12 reference FASTA."""
    cur_dir = Path(__file__).parent
    fasta: Path = cur_dir / "data" / "e_coli_k12.fasta"
    return fasta


@pytest.fixture(scope="session")
def e_coli_k12_fastx_record() -> FastxRecord:
    """Sequence-only FastxRecord that maps to position 80 (0-based) for 80bp on the + strand."""
    sequence = "gttacctgccgtgagtaaattaaaattttattgacttaggtcactaaatactttaaccaatataggcatagcgcacagac"
    return FastxRecord(name="test", sequence=sequence.upper())


@pytest.fixture(scope="session")
def test_fa(e_coli_k12_fasta: Path) -> Path:
    """Alias for e_coli_k12_fasta for compatibility."""
    return e_coli_k12_fasta


@pytest.fixture(scope="session")
def bwa_mem_aligner(e_coli_k12_fasta: Path):
    """A BwaMem aligner initialized with the e. coli K12 reference."""
    from pybwa import BwaMem

    return BwaMem(prefix=e_coli_k12_fasta)


@pytest.fixture(scope="session")
def bwa_aln_aligner(e_coli_k12_fasta: Path):
    """A BwaAln aligner initialized with the e. coli K12 reference."""
    from pybwa import BwaAln

    return BwaAln(prefix=e_coli_k12_fasta)
