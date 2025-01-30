from pathlib import Path
from typing import Optional

import pytest
from fgpyo.sequence import reverse_complement
from pysam import FastxRecord

from pybwa import BwaAln
from pybwa import BwaAlnOptions


def test_bwaaln_options() -> None:
    BwaAlnOptions()
    # TODO: test setting individual options...


def test_bwaaln(e_coli_k12_fasta: Path, e_coli_k12_fastx_record: FastxRecord) -> None:
    opt = BwaAlnOptions()
    bwa = BwaAln(prefix=e_coli_k12_fasta)

    revcomp_seq = (
        None
        if not e_coli_k12_fastx_record.sequence
        else reverse_complement(e_coli_k12_fastx_record.sequence)
    )
    revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)

    recs = bwa.align(opt=opt, queries=[e_coli_k12_fastx_record, revcomp_record])
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


def test_bwaaln_threading(e_coli_k12_fasta: Path, e_coli_k12_fastx_record: FastxRecord) -> None:
    opt = BwaAlnOptions(threads=2)
    bwa = BwaAln(prefix=e_coli_k12_fasta)
    revcomp_seq = (
        None
        if not e_coli_k12_fastx_record.sequence
        else reverse_complement(e_coli_k12_fastx_record.sequence)
    )
    revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)

    queries = [e_coli_k12_fastx_record if i % 2 == 0 else revcomp_record for i in range(100)]
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
    e_coli_k12_fasta: Path,
    e_coli_k12_fastx_record: FastxRecord,
    deletion_length: int,
    max_gap_extensions: Optional[int],
    is_mapped: bool,
) -> None:
    bwa = BwaAln(prefix=e_coli_k12_fasta)

    # Create a deletion
    assert e_coli_k12_fastx_record.sequence is not None
    sequence = (
        e_coli_k12_fastx_record.sequence[:40]
        + e_coli_k12_fastx_record.sequence[40 + deletion_length :]
    )
    e_coli_k12_fastx_record = FastxRecord(name=e_coli_k12_fastx_record.name, sequence=sequence)

    opt = BwaAlnOptions(max_gap_extensions=max_gap_extensions)
    recs = bwa.align(opt=opt, queries=[e_coli_k12_fastx_record])
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
        assert rec.reference_start == -1, e_coli_k12_fastx_record.sequence
        assert rec.is_unmapped


def test_bwa_aln_map_one_multi_mapped_max_hits_one(e_coli_k12_fasta: Path) -> None:
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
    bwa = BwaAln(prefix=e_coli_k12_fasta)
    queries = [FastxRecord(name="NA", sequence="TTTTT")]
    recs = bwa.align(opt=opt, queries=queries)
    assert len(recs) == 1
    rec = recs[0]
    assert rec.has_tag("HN"), str(rec)
    assert rec.get_tag("HN") == 3269888
