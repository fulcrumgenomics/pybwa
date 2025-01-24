from pathlib import Path
from typing import Optional

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

    # Check the files exist
    for suffix in [".bwt", ".sa", ".ann", ".pac"]:
        assert prefix.with_suffix(prefix.suffix + suffix).exists()
    assert prefix.with_suffix(".dict").exists()

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

    # NB: XN, XA not generated for these records
    expected_tags = ["NM", "X0", "X1", "XM", "XO", "XG", "MD", "HN"]
    for rec in recs:
        for tag in expected_tags:
            assert rec.has_tag(tag), f"Missing tag {tag} in: {rec}"


def test_bwaaln_threading(ref_fasta: Path, fastx_record: FastxRecord) -> None:
    opt = BwaAlnOptions(threads=2)
    bwa = BwaAln(prefix=ref_fasta)
    revcomp_seq = None if not fastx_record.sequence else reverse_complement(fastx_record.sequence)
    revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)

    queries = [fastx_record if i % 2 == 0 else revcomp_record for i in range(100)]
    recs = bwa.align(opt=opt, queries=queries)
    assert len(recs) == len(queries)
    for i, rec in enumerate(recs):
        if i % 2 == 0:
            assert rec.query_name == "test"
            assert rec.is_forward
        else:
            assert rec.query_name == "revcomp"
            assert rec.is_reverse
        assert not rec.is_paired
        assert not rec.is_read1
        assert not rec.is_read2
        assert rec.reference_start == 80
        assert rec.cigarstring == "80M"


@pytest.mark.parametrize(
    "max_gap_extensions,is_mapped,deletion_length",
    [
        (None, True, 4),  # 4bp deletions allowed (max diff is 4)
        (None, False, 5),  # 5bp deletions disallowed (max diff is 4)
        (-1, True, 1),  # 1bp deletions allowed
        (-1, False, 2),  # 2bp deletions disallowed
        (0, True, 1),  # 1bp deletions allowed
        (0, False, 2),  # 2bp deletions disallowed
        (4, True, 5),  # 5bp deletion allowed (1 open + 4 extensions)
        (4, False, 6),  # 6bp deletion disallowed (1 open + 4 extensions)
    ],
)
def test_bwaaln_with_deletion(
    ref_fasta: Path,
    fastx_record: FastxRecord,
    deletion_length: int,
    max_gap_extensions: Optional[int],
    is_mapped: bool,
) -> None:
    bwa = BwaAln(prefix=ref_fasta)

    # Create a deletion
    assert fastx_record.sequence is not None
    sequence = fastx_record.sequence[:40] + fastx_record.sequence[40 + deletion_length :]
    fastx_record = FastxRecord(name=fastx_record.name, sequence=sequence)

    opt = BwaAlnOptions(max_gap_extensions=max_gap_extensions)
    recs = bwa.align(opt=opt, queries=[fastx_record])
    assert len(recs) == 1
    rec = recs[0]
    assert rec.query_name == "test"
    assert not rec.is_paired
    assert not rec.is_read1
    assert not rec.is_read2

    if is_mapped:
        assert rec.reference_start == 80
        assert rec.is_forward
        assert rec.cigarstring == f"40M{deletion_length}D{40-deletion_length}M"
    else:
        assert rec.reference_start == -1, fastx_record.sequence
        assert rec.is_unmapped


def test_bwa_aln_map_one_multi_mapped_max_hits_one(ref_fasta: Path) -> None:
    """Tests that a query that returns too many hits (>max_hits) returns the number of hits but
    not the list of hits themselves."""
    opt = BwaAlnOptions(
        threads=2,
        max_hits=1,
        max_mismatches=3,
        max_mismatches_in_seed=3,
        max_gap_opens=0,
        max_gap_extensions=-1,
        min_indel_to_end_distance=3,
        seed_length=20,
        find_all_hits=True,
        with_md=True,
    )
    bwa = BwaAln(prefix=ref_fasta)
    queries = [FastxRecord(name="NA", sequence="TTTTT")]
    recs = bwa.align(opt=opt, queries=queries)
    assert len(recs) == 1
    rec = recs[0]
    assert rec.has_tag("HN"), str(rec)
    assert rec.get_tag("HN") == 3269888


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
    opt = BwaMemOptions(with_xr_tag=True)
    bwa = BwaMem(prefix=ref_fasta)

    revcomp_seq = None if not fastx_record.sequence else reverse_complement(fastx_record.sequence)
    revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)

    recs_of_recs = bwa.align(opt=opt, queries=[fastx_record, revcomp_record])
    assert len(recs_of_recs) == 2

    assert len(recs_of_recs[0]) == 1
    rec = recs_of_recs[0][0]
    assert rec.query_name == "test"
    assert not rec.is_paired
    assert not rec.is_read1
    assert not rec.is_read2
    assert rec.reference_start == 80
    assert rec.is_forward
    assert rec.cigarstring == "80M"

    assert len(recs_of_recs[1]) == 1
    rec = recs_of_recs[1][0]
    assert rec.query_name == "revcomp"
    assert not rec.is_paired
    assert not rec.is_read1
    assert not rec.is_read2
    assert rec.reference_start == 80
    assert rec.is_reverse
    assert rec.cigarstring == "80M"
    # TODO: test multi-mapping, reverse strand, etc

    # NB: XA amd XB not generated for these records
    expected_tags = ["NM", "MD", "AS", "XS", "XR"]
    for recs in recs_of_recs:
        for rec in recs:
            for tag in expected_tags:
                assert rec.has_tag(tag), f"Missing tag {tag} in: {rec}"


def test_bwamem_threading(ref_fasta: Path, fastx_record: FastxRecord) -> None:
    opt = BwaMemOptions(threads=2)
    bwa = BwaMem(prefix=ref_fasta)

    revcomp_seq = None if not fastx_record.sequence else reverse_complement(fastx_record.sequence)
    revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)

    queries = [fastx_record if i % 2 == 0 else revcomp_record for i in range(100)]
    list_of_recs = bwa.align(opt=opt, queries=queries)
    assert len(list_of_recs) == len(queries)
    for i, recs in enumerate(list_of_recs):
        assert len(recs) == 1
        rec = recs[0]
        if i % 2 == 0:
            assert rec.query_name == "test"
            assert rec.is_forward
        else:
            assert rec.query_name == "revcomp"
            assert rec.is_reverse
        assert not rec.is_paired
        assert not rec.is_read1
        assert not rec.is_read2
        assert rec.reference_start == 80
        assert rec.cigarstring == "80M"
