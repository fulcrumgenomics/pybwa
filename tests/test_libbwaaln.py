from pathlib import Path
from typing import Optional
from typing import Union

import pysam
import pytest
from fgpyo.sam import Cigar
from fgpyo.sam.builder import SamBuilder
from fgpyo.sequence import reverse_complement
from pysam import FastxRecord

from pybwa import BwaAln
from pybwa import BwaAlnOptions
from pybwa import BwaIndex
from pybwa.libbwaaln import AuxHit
from pybwa.libbwaaln import to_xa_hits


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
        max_entries=13,
        threads=14,
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
    assert options.max_entries == 13
    assert options.threads == 14

    options.log_scaled_gap_penalty = True
    assert options.log_scaled_gap_penalty

    options.find_all_hits = True
    assert options.find_all_hits
    assert options.stop_at_max_best_hits == 0x7FFFFFFF

    options.with_md = True
    assert options.with_md


def test_bwaaln_options_default() -> None:
    from pybwa.libbwaaln import _assert_gap_opt_are_the_same

    _assert_gap_opt_are_the_same()


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
    Tests running `bwa aln` without anything fancy.

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
    bwa.reinitialize_seed()
    assert e_coli_k12_fastx_record.sequence is not None

    queries = [
        FastxRecord(
            name=e_coli_k12_fastx_record.name, sequence=e_coli_k12_fastx_record.sequence[:10]
        )
    ]

    recs = bwa.align(opt=opt, queries=queries)
    assert len(recs) == 1
    rec = recs[0]
    assert rec.reference_start == 733390
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

    # parse the XA values once
    xa_hits: list[AuxHit] = to_xa_hits(rec)
    assert len(xa_hits) == 416

    # test when we parse it again
    xa_hits = to_xa_hits(rec)
    assert len(xa_hits) == 416


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
    bwa.reinitialize_seed()
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

    if max_gap_extensions is None:
        opt = BwaAlnOptions()
    else:
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
        assert rec.cigarstring == f"40M{deletion_length}D{40 - deletion_length}M"
    else:
        assert rec.reference_start == -1, e_coli_k12_fastx_record.sequence
        assert rec.is_unmapped


def test_bwa_aln_map_one_multi_mapped_max_hits_one(e_coli_k12_fasta: Path) -> None:
    # Tests that a query that returns too many hits (>max_hits) returns the number of hits but
    # not the list of hits themselves.
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
    bwa.reinitialize_seed()
    assert len(recs) == 1
    rec = recs[0]
    assert rec.has_tag("HN"), str(rec)
    assert rec.get_tag("HN") == 3269888


@pytest.mark.parametrize("num_amb", [0, 1, 5])
def test_bwa_aln_ambiguous_bases(num_amb: int, tmp_path_factory: pytest.TempPathFactory) -> None:
    sequence = ("A" * 100) + ("N" * num_amb) + ("T" * 100)

    src_dir = Path(str(tmp_path_factory.mktemp("test_bwa_aln_ambiguous_bases")))
    fasta = src_dir / "ref.fasta"
    with fasta.open("w") as writer:
        writer.write(f">ref\n{sequence}\n")
    BwaIndex.index(fasta=fasta)

    bwa = BwaAln(prefix=fasta)
    queries = [FastxRecord(name="NA", sequence=sequence)]
    recs = bwa.align(queries=queries)
    assert len(recs) == 1
    rec = recs[0]
    if num_amb == 0:
        assert not rec.has_tag("XN"), str(rec)
    else:
        assert rec.has_tag("XN"), str(rec)
        assert rec.get_tag("XN") == num_amb


@pytest.fixture()
def xa_single_hit() -> str:
    return "chr4,-97592047,24M,3;"


def _assert_single_hit(hit: AuxHit, md: Optional[str] = None, rest: Optional[str] = None) -> None:
    assert hit.refname == "chr4"
    assert hit.start == 97592047
    assert hit.negative
    assert f"{hit.cigar}" == "24M"
    assert hit.edits == 3
    assert hit.md == md
    assert hit.rest == rest


@pytest.fixture()
def xa_double_hit() -> str:
    return "chr4,-97592047,24M,3;chr8,+32368749,32M,4;"


def _assert_double_hit(hits: list[AuxHit]) -> None:
    assert len(hits) == 2
    _assert_single_hit(hits[0])

    hit: AuxHit = hits[1]
    assert hit.refname == "chr8"
    assert hit.start == 32368749
    assert not hit.negative
    assert f"{hit.cigar}" == "32M"
    assert hit.edits == 4
    assert hit.md is None
    assert hit.rest is None


def test_to_xa_hits_single_from_string(xa_single_hit: str) -> None:
    hits: list[AuxHit] = to_xa_hits(xa_single_hit)
    assert len(hits) == 1
    _assert_single_hit(hits[0])


def test_to_xa_hits_single_from_bytes(xa_single_hit: str) -> None:
    hits: list[AuxHit] = to_xa_hits(xa_single_hit.encode("utf-8"))
    assert len(hits) == 1
    _assert_single_hit(hits[0])


def test_to_xa_hits_single_from_alignedsegment(xa_single_hit: str) -> None:
    builder = SamBuilder()
    rec, _ = builder.add_pair()
    rec.set_tag("XA", xa_single_hit)
    hits: list[AuxHit] = to_xa_hits(rec)
    assert len(hits) == 1
    _assert_single_hit(hits[0])


def test_to_xa_hits_multi(xa_double_hit: str) -> None:
    hits: list[AuxHit] = to_xa_hits(xa_double_hit)
    _assert_double_hit(hits)


def test_to_xa_hits_max_hits(xa_double_hit: str) -> None:
    hits: list[AuxHit]

    # Should return only one hit
    hits = to_xa_hits(xa_double_hit, max_hits=1)
    assert len(hits) == 1
    _assert_single_hit(hits[0])

    # Should return two hits
    for max_hits in [0, 2, 3]:
        hits = to_xa_hits(xa_double_hit, max_hits=max_hits)
        _assert_double_hit(hits)


def test_to_xa_hits_with_md(xa_single_hit: str) -> None:
    xa: str = f"{xa_single_hit[:-1]},24;"
    hits: list[AuxHit] = to_xa_hits(xa)
    assert len(hits) == 1
    hit: AuxHit = hits[0]
    _assert_single_hit(hit, md="24")


def test_to_xa_hits_with_md_and_rest(xa_single_hit: str) -> None:
    xa: str = f"{xa_single_hit[:-1]},24,rest;"
    hits: list[AuxHit] = to_xa_hits(xa)
    assert len(hits) == 1
    hit: AuxHit = hits[0]
    _assert_single_hit(hit, md="24", rest="rest")


def test_to_xa_hits_missing_trailing_semicolon(xa_single_hit: str) -> None:
    with pytest.raises(ValueError, match="Could not parse XA tag"):
        to_xa_hits(xa_single_hit[:-1])


def test_to_xa_hits_too_few_fields() -> None:
    xa: str = "chr4,-97592047,24M"
    with pytest.raises(ValueError, match="Could not parse XA tag"):
        to_xa_hits(xa)


def test_to_xa_hits_single_semicolon() -> None:
    xa: str = ";"
    with pytest.raises(ValueError, match="Could not parse XA tag"):
        to_xa_hits(xa)


def test_to_xa_hits_empty() -> None:
    xa: str = ""
    with pytest.raises(ValueError, match="Could not parse XA tag"):
        to_xa_hits(xa)


def test_to_xa_hits_parse_int() -> None:
    xa: str = "chr4,-ABCDEF,24M,3;"
    with pytest.raises(ValueError, match="Could not parse XA tag"):
        to_xa_hits(xa)
    xa = "chr4,-ABCDEF,24M,;"
    with pytest.raises(ValueError, match="Could not parse XA tag"):
        to_xa_hits(xa)


def test_to_xa_hits_parse_error() -> None:
    xa: str = "chr4,-97592047,24M,3;"
    while len(xa) > 0:
        xa = xa[:-1]
        with pytest.raises(ValueError, match="Could not parse XA tag"):
            to_xa_hits(xa)


def test_to_xa_hits_no_xa() -> None:
    builder = SamBuilder()
    rec, _ = builder.add_pair()
    assert not rec.has_tag("XA")
    assert len(to_xa_hits(rec)) == 0


def _to_xa_hit(refname: str, start: int, negative: bool, cigar: str, edits: int) -> AuxHit:
    return AuxHit(
        refname=refname,
        start=start,
        negative=negative,
        cigar=Cigar.from_cigarstring(cigar),
        edits=edits,
        md=None,
        rest=None,
    )


def test_xa_hit_mismatches_indels_gt_edits() -> None:
    hit = _to_xa_hit("chr1", 1081, False, "20M5I30M", 0)
    with pytest.raises(ValueError, match="indel_sum"):
        assert 0 == hit.mismatches


def test_xa_hit_mismatches() -> None:
    hit = _to_xa_hit("chr1", 1081, False, "20M5I30M", 5)
    assert 0 == hit.mismatches
    hit = _to_xa_hit("chr1", 1081, False, "22M3I30M", 5)
    assert 2 == hit.mismatches
    hit = _to_xa_hit("chr1", 1081, False, "55M", 5)
    assert 5 == hit.mismatches


def test_xa_hit_end() -> None:
    hit = _to_xa_hit("chr1", 1081, False, "50M", 5)
    assert hit.end == 1130
    hit = _to_xa_hit("chr1", 1081, False, "20M5I30M", 5)
    assert hit.end == 1130
    hit = _to_xa_hit("chr1", 1081, False, "20M5D30M", 5)
    assert hit.end == 1135


def test_xa_hit_str() -> None:
    hit = _to_xa_hit("chr1", 1081, False, "20M5I30M", 5)
    assert str(hit) == "chr1,+1081,20M5I30M,5"
    hit = _to_xa_hit("chr1", 1081, True, "20M5I30M", 5)
    assert str(hit) == "chr1,-1081,20M5I30M,5"


@pytest.mark.parametrize("limit", [1, 32767, 32768, 100000])
def test_bwa_aln_max_hits(limit: int, tmp_path_factory: pytest.TempPathFactory) -> None:
    query = "G" * 25
    target = ("A" * 25) + query + ("T" * 25)

    src_dir = Path(str(tmp_path_factory.mktemp("test_bwa_aln_max_hits")))
    fasta = src_dir / "ref.fasta"
    with fasta.open("w") as writer:
        writer.write(">ref\n")
        for _ in range(limit):
            writer.write(f"{target}\n")
    BwaIndex.index(fasta=fasta)

    bwa = BwaAln(prefix=fasta)
    queries = [FastxRecord(name="NA", sequence=query)]
    options = BwaAlnOptions(seed_length=len(query), max_mismatches=0, max_hits=limit)
    recs = bwa.align(queries=queries, opt=options)
    assert len(recs) == 1
    rec = recs[0]
    assert rec.is_unmapped is False
    assert rec.get_tag("HN") == limit, f"limit was {limit}"


def test_bwa_aln_max_hits_part_deux(tmp_path_factory: pytest.TempPathFactory) -> None:
    query = "G" * 25
    target = ("A" * 25) + query + ("T" * 25)

    src_dir = Path(str(tmp_path_factory.mktemp("test_bwa_aln_max_hits")))
    fasta = src_dir / "ref.fasta"
    with fasta.open("w") as writer:
        writer.write(">ref\n")
        for _ in range(1000):
            writer.write(f"{target}\n")
    BwaIndex.index(fasta=fasta)

    bwa = BwaAln(prefix=fasta)
    queries = [FastxRecord(name="NA", sequence=query)]
    limit = 2147483647
    options = BwaAlnOptions(seed_length=len(query), max_mismatches=0, max_hits=limit)
    recs = bwa.align(queries=queries, opt=options)
    assert len(recs) == 1
    rec = recs[0]
    assert rec.is_unmapped is False
    assert rec.get_tag("HN") == 1000, f"limit was {limit}"

    limit += 1
    options = BwaAlnOptions(seed_length=len(query), max_mismatches=0, max_hits=limit)
    with pytest.raises(OverflowError, match="value too large to convert to int"):
        bwa.align(queries=queries, opt=options)


@pytest.mark.parametrize(
    "max_mismatches,max_entries,should_fail",
    [
        # with 0 mismatches, it should always succeed
        (0, 1, False),
        (1, 1, True),
        (1, 167, True),
        (1, 168, False),
        (2, 1, True),
        (2, 204, True),
        (2, 205, False),
    ],
)
def test_bwa_aln_exit_early(
    e_coli_k12_fasta: Path,
    e_coli_k12_fastx_record: FastxRecord,
    max_mismatches: int,
    max_entries: int,
    should_fail: bool,
) -> None:
    opt = BwaAlnOptions(threads=2, max_entries=max_entries, max_mismatches=max_mismatches)
    bwa = BwaAln(prefix=e_coli_k12_fasta)
    queries = [e_coli_k12_fastx_record]
    recs = bwa.align(opt=opt, queries=queries)
    assert len(recs) == 1
    rec = recs[0]
    assert rec.query_name == "test"
    if should_fail:
        assert rec.is_unmapped is True
        assert rec.has_tag("ea")
    else:
        assert rec.is_forward
        assert not rec.is_paired
        assert not rec.is_read1
        assert not rec.is_read2
        assert rec.reference_start == 80
        assert rec.cigarstring == "80M"
