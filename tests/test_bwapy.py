from pathlib import Path

import pytest
from fgpyo.sequence import reverse_complement
from pysam import FastxRecord

from pybwa import BwaAln
from pybwa import BwaAlnOptions
from pybwa import BwaIndex
from pybwa.libbwamem import BwaMem
from pybwa.libbwamem import BwaMemOptions


@pytest.fixture()
def ref_fasta() -> Path:
    cur_dir = Path(__file__).parent
    fasta: Path = cur_dir / "data" / "e_coli_k12.fasta"
    return fasta


def test_bwa_index_build(ref_fasta: Path, tmp_path_factory: pytest.TempPathFactory) -> None:
    tmp_dir = Path(str(tmp_path_factory.mktemp("test_bwa_index_build")))
    prefix = tmp_dir / ref_fasta.name

    # Build the index
    BwaIndex.index(fasta=ref_fasta, prefix=prefix)
    # Load it
    BwaIndex(prefix=prefix)


@pytest.fixture(scope="function")
def temp_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("test_vcf")


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

    revcomp_seq = None if not fastx_record.sequence else reverse_complement(fastx_record.sequence)
    revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)

    recs = bwa.align(opt=opt, queries=[fastx_record, revcomp_record])
    assert len(recs) == 2
    rec = recs[0]
    assert rec.query_name == "test"
    assert not rec.is_paired
    assert not rec.is_read1
    assert not rec.is_read2
    assert rec.reference_start == 80
    assert rec.is_forward
    assert rec.cigarstring == "80M"

    rec = recs[1]
    assert rec.query_name == "revcomp"
    assert not rec.is_paired
    assert not rec.is_read1
    assert not rec.is_read2
    assert rec.reference_start == 80
    assert rec.is_reverse
    assert rec.cigarstring == "80M"


def test_bwamem_options() -> None:
    # default options
    options = BwaMemOptions()
    assert options.min_seed_len == 19

    # finalize, returning a copy
    copy = options.finalize(copy=True)
    assert copy.min_seed_len == 19
    assert copy.finalized
    assert not options.finalized

    # update min seed len
    options.min_seed_len = 20
    assert options.min_seed_len == 20

    # finalize, returning a copy with a new value
    copy = options.finalize(copy=True)
    assert copy.min_seed_len == 20
    assert copy.finalized
    assert not options.finalized

    # finalize, returning itself finalized
    options.finalize()
    assert options.finalized
    assert options.min_seed_len == 20  # type:ignore[unreachable]

    # raise an exception if we try to set a new value
    with pytest.raises(AttributeError):
        options.min_seed_len = 19


def test_bwamem(ref_fasta: Path, fastx_record: FastxRecord) -> None:
    opt = BwaMemOptions()
    bwa = BwaMem(prefix=ref_fasta)

    revcomp_seq = None if not fastx_record.sequence else reverse_complement(fastx_record.sequence)
    revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)

    recs = bwa.align(opt=opt, queries=[fastx_record, revcomp_record])
    assert len(recs) == 2

    assert len(recs[0]) == 1
    rec = recs[0][0]
    assert rec.query_name == "test"
    assert not rec.is_paired
    assert not rec.is_read1
    assert not rec.is_read2
    assert rec.reference_start == 80
    assert rec.is_forward
    assert rec.cigarstring == "80M"

    assert len(recs[1]) == 1
    rec = recs[1][0]
    assert rec.query_name == "revcomp"
    assert not rec.is_paired
    assert not rec.is_read1
    assert not rec.is_read2
    assert rec.reference_start == 80
    assert rec.is_reverse
    assert rec.cigarstring == "80M"
    # TODO: test multi-mapping, reverse strand, etc
