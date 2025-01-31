from pathlib import Path
from typing import Optional
from typing import Union

import pysam
import pytest
from fgpyo.sequence import reverse_complement
from pysam import FastxRecord

from pybwa import BwaAln
from pybwa import BwaAlnOptions
from pybwa import BwaIndex


def test_bwaaln_options() -> None:
    # default
    BwaAlnOptions()

    options = BwaAlnOptions(
        max_mismatches=1,
        max_gap_opens=2,
        max_gap_extensions=3,
        min_indel_to_end_distance=4,
        max_occurrences_for_extending_long_deletion=5,
        seed_length=6,
        max_mismatches_in_seed=7,
        mismatch_penalty=8,
        gap_open_penalty=9,
        gap_extension_penalty=10,
        stop_at_max_best_hits=11,
        max_hits=12,
        log_scaled_gap_penalty=False,
        find_all_hits=False,
        with_md=False,
        threads=13,
    )

    assert options.max_mismatches == 1
    assert options.max_gap_opens == 2
    assert options.max_gap_extensions == 3
    assert options.min_indel_to_end_distance == 4
    assert options.max_occurrences_for_extending_long_deletion == 5
    assert options.seed_length == 6
    assert options.max_mismatches_in_seed == 7
    assert options.mismatch_penalty == 8
    assert options.gap_open_penalty == 9
    assert options.gap_extension_penalty == 10
    assert options.stop_at_max_best_hits == 11
    assert options.max_hits == 12
    assert not options.log_scaled_gap_penalty
    assert not options.find_all_hits
    assert not options.with_md
    assert options.threads == 13

    options.log_scaled_gap_penalty = True
    assert options.log_scaled_gap_penalty

    options.find_all_hits = True
    assert options.find_all_hits
    assert options.stop_at_max_best_hits == 0x7FFFFFFF

    options.with_md = True
    assert options.with_md


def test_bwaaln_load_index(e_coli_k12_fasta: Path) -> None:
    BwaAln(prefix=e_coli_k12_fasta)
    BwaAln(index=BwaIndex(e_coli_k12_fasta))
    with pytest.raises(ValueError):
        BwaAln()


def test_bwaaln_no_queries(e_coli_k12_fasta: Path) -> None:
    bwa = BwaAln(prefix=e_coli_k12_fasta)
    assert bwa.align(queries=[]) == []


@pytest.mark.parametrize("query_is_str", [True, False])
@pytest.mark.parametrize("with_quals", [True, False])
def test_bwaaln_basic(
    e_coli_k12_fasta: Path,
    e_coli_k12_fastx_record: FastxRecord,
    query_is_str: bool,
    with_quals: bool,
) -> None:
    """
    Args:
        e_coli_k12_fasta: path to the index FASTA
        e_coli_k12_fastx_record: the FASTX record to align
        query_is_str: true if the query is to be passed a string, false if as a FastxRecord
        with_quals: true to provide base qualities, false otherwise. Only works if query_is_str is
                    false.
    """
    opt = BwaAlnOptions()
    bwa = BwaAln(prefix=e_coli_k12_fasta)

    revcomp_seq = (
        None
        if not e_coli_k12_fastx_record.sequence
        else reverse_complement(e_coli_k12_fastx_record.sequence)
    )
    if with_quals and not query_is_str:
        assert e_coli_k12_fastx_record.sequence is not None
        quality = "J"
        quality += "I" * (len(e_coli_k12_fastx_record.sequence) - 1)
        revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq, quality=quality)
        e_coli_k12_fastx_record = FastxRecord(
            name=e_coli_k12_fastx_record.name,
            sequence=e_coli_k12_fastx_record.sequence,
            quality=quality,
        )
    else:
        revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)

    queries: Union[list[str], list[FastxRecord]]
    if query_is_str:
        queries = [str(q.sequence) for q in [e_coli_k12_fastx_record, revcomp_record]]
    else:
        queries = [e_coli_k12_fastx_record, revcomp_record]

    recs = bwa.align(opt=opt, queries=queries)
    assert len(recs) == 2
    rec = recs[0]
    assert rec.query_name == "read.0" if query_is_str else "test"
    assert rec.query_sequence == e_coli_k12_fastx_record.sequence
    if e_coli_k12_fastx_record.quality is not None:
        assert rec.query_qualities is not None
        assert (
            pysam.qualities_to_qualitystring(rec.query_qualities) == e_coli_k12_fastx_record.quality
        )
    assert not rec.is_paired
    assert not rec.is_read1
    assert not rec.is_read2
    assert rec.reference_start == 80
    assert rec.is_forward
    assert rec.cigarstring == "80M"

    rec = recs[1]
    assert rec.query_name == "read.1" if query_is_str else "revcomp"
    assert rec.query_sequence == e_coli_k12_fastx_record.sequence
    if revcomp_record.quality is not None:
        assert rec.query_qualities is not None
        assert pysam.qualities_to_qualitystring(rec.query_qualities) == revcomp_record.quality[::-1]
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


def test_bwaaln_multi_hit(e_coli_k12_fasta: Path, e_coli_k12_fastx_record: FastxRecord) -> None:
    opt = BwaAlnOptions(threads=1, seed_length=10, max_hits=1000)
    bwa = BwaAln(prefix=e_coli_k12_fasta)
    assert e_coli_k12_fastx_record.sequence is not None

    queries = [
        FastxRecord(
            name=e_coli_k12_fastx_record.name, sequence=e_coli_k12_fastx_record.sequence[:10]
        )
    ]

    recs = bwa.align(opt=opt, queries=queries)
    assert len(recs) == 1
    rec = recs[0]
    assert rec.reference_start == 4074272
    assert rec.query_name == "test"
    assert queries[0].sequence is not None
    assert rec.query_sequence == reverse_complement(queries[0].sequence)
    assert not rec.is_read1
    assert not rec.is_read2
    assert rec.is_reverse
    assert rec.cigarstring == "10M"
    assert rec.has_tag("XA")
    assert rec.has_tag("HN")
    assert rec.get_tag("HN") == 417
    xa_values = str(rec.get_tag("XA")).split(";")[:-1]
    assert len(xa_values) == 416


def test_bwaaln_multi_hit_split_cigar(
    e_coli_k12_fasta: Path, e_coli_k12_fastx_record: FastxRecord
) -> None:
    opt = BwaAlnOptions(
        threads=1,
        seed_length=5,
        max_hits=10000,
        max_gap_opens=1,
        max_gap_extensions=3,
        gap_open_penalty=1,
    )
    bwa = BwaAln(prefix=e_coli_k12_fasta)
    assert e_coli_k12_fastx_record.sequence is not None

    queries = [
        FastxRecord(
            name=e_coli_k12_fastx_record.name,
            sequence=e_coli_k12_fastx_record.sequence[:5]
            + "AAA"
            + e_coli_k12_fastx_record.sequence[5:10],
        )
    ]

    recs = bwa.align(opt=opt, queries=queries)
    assert len(recs) == 1
    rec = recs[0]
    assert rec.reference_start == 465895
    assert rec.query_name == "test"
    assert queries[0].sequence is not None
    assert rec.query_sequence == reverse_complement(queries[0].sequence)
    assert not rec.is_read1
    assert not rec.is_read2
    assert rec.is_reverse
    assert rec.cigarstring == "5M1I7M"
    assert rec.has_tag("XA")
    assert rec.has_tag("HN")
    assert rec.get_tag("HN") == 8
    xa_values = str(rec.get_tag("XA")).split(";")[:-1]
    assert len(xa_values) == 7, xa_values
    # check that we have one XA value with an indel
    assert any(v.split(",")[2] == "4M1I8M" for v in xa_values)


def test_bwa_unmapped(
    e_coli_k12_fasta: Path,
) -> None:
    opt = BwaAlnOptions()
    bwa = BwaAln(prefix=e_coli_k12_fasta)
    recs = bwa.align(opt=opt, queries=["GATTACAGATTACAGATTACAGATTACA"])
    assert len(recs) == 1
    rec = recs[0]
    assert rec.is_unmapped


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
