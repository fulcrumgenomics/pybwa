from pathlib import Path

import pytest
from fgpyo.sequence import reverse_complement
from pysam import FastxRecord

from pybwa import BwaIndex
from pybwa.libbwamem import BwaMem
from pybwa.libbwamem import BwaMemMode
from pybwa.libbwamem import BwaMemOptions


def test_bwamem_options_not_finalized() -> None:
    opt = BwaMemOptions().finalize()
    with pytest.raises(AttributeError):
        opt.min_seed_len = 2


def test_bwamem_options_set() -> None:
    opt = BwaMemOptions(
        min_seed_len=2,
        # mode: BwaMemMode | None = None,
        band_width=2,
        match_score=2,
        mismatch_penalty=2,
        minimum_score=2,
        unpaired_penalty=2,
        n_threads=2,
        skip_pairing=True,
        output_all_for_fragments=True,
        interleaved_paired_end=True,
        short_split_as_secondary=True,
        skip_mate_rescue=True,
        soft_clip_supplementary=True,
        with_xr_tag=True,
        query_coord_as_primary=True,
        keep_mapq_for_supplementary=True,
        with_xb_tag=True,
        max_occurrences=2,
        off_diagonal_x_dropoff=2,
        ignore_alternate_contigs=True,
        internal_seed_split_factor=2.1,
        drop_chain_fraction=2.1,
        max_mate_rescue_rounds=2,
        min_seeded_bases_in_chain=2,
        seed_occurrence_in_3rd_round=2,
        xa_max_hits=2,
        xa_drop_ratio=2.1,
        gap_open_penalty=2,
        gap_extension_penalty=2,
        clipping_penalty=2,
        threads=2,
        chunk_size=2,
    )

    # test the getters
    assert opt.min_seed_len == 2
    assert opt.band_width == 2
    assert opt.match_score == 2
    assert opt.mismatch_penalty == 2
    assert opt.minimum_score == 2
    assert opt.unpaired_penalty == 2
    assert opt.n_threads == 2
    assert opt.skip_pairing
    assert opt.output_all_for_fragments
    assert opt.interleaved_paired_end
    assert opt.short_split_as_secondary
    assert opt.skip_mate_rescue
    assert opt.soft_clip_supplementary
    assert opt.with_xr_tag
    assert opt.query_coord_as_primary
    assert opt.keep_mapq_for_supplementary
    assert opt.with_xb_tag
    assert opt.max_occurrences == 2
    assert opt.off_diagonal_x_dropoff == 2
    assert opt.ignore_alternate_contigs
    assert opt.internal_seed_split_factor == pytest.approx(2.1, 0.01)
    assert opt.drop_chain_fraction == pytest.approx(2.1, 0.01)
    assert opt.max_mate_rescue_rounds == 2
    assert opt.min_seeded_bases_in_chain == 2
    assert opt.seed_occurrence_in_3rd_round == 2
    assert opt.xa_max_hits == 2
    assert opt.xa_drop_ratio == pytest.approx(2.1, 0.01)
    assert opt.gap_open_penalty == 2
    assert opt.gap_extension_penalty == 2
    assert opt.clipping_penalty == 2
    assert opt.threads == 2
    assert opt.chunk_size == 2
    assert opt.mode is None

    # test that boolean flags work
    opt.skip_pairing = False
    assert not opt.skip_pairing
    opt.skip_pairing = True
    assert opt.skip_pairing

    # Test scaling

    # not scaled because it has been changed
    assert opt.gap_open_penalty == 2
    assert opt.gap_extension_penalty == 2
    assert not opt.finalized
    opt = opt.finalize()
    assert opt.finalized
    # even after finalize
    assert opt.gap_open_penalty == 2
    assert opt.gap_extension_penalty == 2


def test_bwa_mem_scaling() -> None:
    opt = BwaMemOptions(
        match_score=5,
    )
    # not scaled because it hasn't been finalized
    assert opt.gap_open_penalty == 6
    assert opt.gap_extension_penalty == 1
    opt = opt.finalize(copy=True)
    # scaled because it wasn't set
    assert opt.gap_open_penalty == 30
    assert opt.gap_extension_penalty == 5

    # test tuples
    opt = BwaMemOptions(clipping_penalty=(5, 3)).finalize()
    assert opt.clipping_penalty == (5, 3)


def test_bwamem_tuples_options() -> None:
    opt = BwaMemOptions(
        xa_max_hits=(2, 2),
        gap_open_penalty=(3, 3),
        gap_extension_penalty=(4, 4),
        clipping_penalty=(5, 5),
    )

    assert opt.xa_max_hits == 2
    assert opt.gap_open_penalty == 3
    assert opt.gap_extension_penalty == 4
    assert opt.clipping_penalty == 5

    opt = BwaMemOptions(
        xa_max_hits=(2, 3),
        gap_open_penalty=(4, 5),
        gap_extension_penalty=(6, 7),
        clipping_penalty=(8, 9),
    )
    assert opt.xa_max_hits == (2, 3)
    assert opt.gap_open_penalty == (4, 5)
    assert opt.gap_extension_penalty == (6, 7)
    assert opt.clipping_penalty == (8, 9)


def test_bwamem_options_mode() -> None:
    opt = BwaMemOptions(mode=BwaMemMode.INTRACTG).finalize()
    assert opt.clipping_penalty == 5

    opt = BwaMemOptions(mode=BwaMemMode.ONT2D).finalize()
    assert opt.clipping_penalty == 0
    assert opt.min_seeded_bases_in_chain == 20

    opt = BwaMemOptions(mode=BwaMemMode.PACBIO).finalize()
    assert opt.clipping_penalty == 0
    assert opt.min_seeded_bases_in_chain == 40


def test_bwamem(e_coli_k12_fasta: Path, e_coli_k12_fastx_record: FastxRecord) -> None:
    opt = BwaMemOptions(with_xr_tag=True)
    bwa = BwaMem(prefix=e_coli_k12_fasta)

    revcomp_seq = (
        None
        if not e_coli_k12_fastx_record.sequence
        else reverse_complement(e_coli_k12_fastx_record.sequence)
    )
    revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)

    recs_of_recs = bwa.align(opt=opt, queries=[e_coli_k12_fastx_record, revcomp_record])
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


def test_bwamem_from_index(e_coli_k12_fasta: Path, e_coli_k12_fastx_record: FastxRecord) -> None:
    index = BwaIndex(prefix=e_coli_k12_fasta)
    BwaMem(index=index)


def test_bwamem_no_index() -> None:
    with pytest.raises(ValueError):
        BwaMem()


def test_bwamem_align_no_options(
    e_coli_k12_fasta: Path, e_coli_k12_fastx_record: FastxRecord
) -> None:
    bwa = BwaMem(prefix=e_coli_k12_fasta)
    queries = [e_coli_k12_fastx_record]
    list_of_recs = bwa.align(queries=queries)
    assert len(list_of_recs) == 1
    assert len(list_of_recs[0]) == 1


def test_bwamem_align_no_queries(e_coli_k12_fasta: Path) -> None:
    bwa = BwaMem(prefix=e_coli_k12_fasta)
    list_of_recs = bwa.align(queries=[])
    assert len(list_of_recs) == 0


def test_bwamem_align_from_string(
    e_coli_k12_fasta: Path, e_coli_k12_fastx_record: FastxRecord
) -> None:
    opt = BwaMemOptions(threads=2)
    bwa = BwaMem(prefix=e_coli_k12_fasta)

    revcomp_seq = (
        None
        if not e_coli_k12_fastx_record.sequence
        else reverse_complement(e_coli_k12_fastx_record.sequence)
    )
    revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)
    assert e_coli_k12_fastx_record.sequence is not None
    assert revcomp_record.sequence is not None

    queries: list[str] = [
        e_coli_k12_fastx_record.sequence if i % 2 == 0 else revcomp_record.sequence
        for i in range(100)
    ]
    list_of_recs = bwa.align(opt=opt, queries=queries)
    assert len(list_of_recs) == len(queries)
    for i, recs in enumerate(list_of_recs):
        assert len(recs) == 1
        rec = recs[0]
        assert rec.query_name == f"read.{i}"
        if i % 2 == 0:
            assert rec.is_forward
        else:
            assert rec.is_reverse
        assert not rec.is_paired
        assert not rec.is_read1
        assert not rec.is_read2
        assert rec.reference_start == 80
        assert rec.cigarstring == "80M"


def test_bwamem_threading(e_coli_k12_fasta: Path, e_coli_k12_fastx_record: FastxRecord) -> None:
    opt = BwaMemOptions(threads=2)
    bwa = BwaMem(prefix=e_coli_k12_fasta)

    revcomp_seq = (
        None
        if not e_coli_k12_fastx_record.sequence
        else reverse_complement(e_coli_k12_fastx_record.sequence)
    )
    revcomp_record = FastxRecord(name="revcomp", sequence=revcomp_seq)

    queries = [e_coli_k12_fastx_record if i % 2 == 0 else revcomp_record for i in range(100)]
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


def test_bwamem_align_multi_map(
    e_coli_k12_fasta: Path, e_coli_k12_fastx_record: FastxRecord
) -> None:
    opt = BwaMemOptions(
        min_seed_len=10,
        minimum_score=10,
        output_all_for_fragments=True,
        short_split_as_secondary=True,
    )
    bwa = BwaMem(prefix=e_coli_k12_fasta)

    assert e_coli_k12_fastx_record.sequence is not None
    query = FastxRecord(
        name=e_coli_k12_fastx_record.name,
        sequence=e_coli_k12_fastx_record.sequence[:10] + e_coli_k12_fastx_record.sequence[-10:],
    )
    list_of_recs = bwa.align(opt=opt, queries=[query])
    assert len(list_of_recs) == 1
    assert len(list_of_recs[0]) == 27
    rec = list_of_recs[0][0]
    assert rec.reference_start == 3618128
    assert rec.cigarstring == "14M6S"
    score = int(rec.get_tag("AS"))
    for rec in list_of_recs[0]:
        assert score >= int(rec.get_tag("AS"))


@pytest.mark.parametrize("query_coord_as_primary", [True])
def test_bwamem_align_sa_tag(e_coli_k12_fasta: Path, query_coord_as_primary: bool) -> None:
    # alignments below are defined based on (1) alignment score, then (2) genomic order
    # alignment 1
    seq1 = "agcttttcattctgactgcaacgggcaatatgtctctgtgtggattaaaaaaagagtgtctgatagcagcttctgaactg"
    pos1 = 1759
    cig1 = "79S81M80S"
    sap1 = "U00096.2,1760,+,79S81M80S,60,0"
    # alignment 2
    seq2 = "cgttggcggtgcgctgctggagcaactgaagcgtcagcaaagctggctgaagaataaacatatcgacttacgtgtctgcg"
    pos2 = 0
    cig2 = "80M160S"
    sap2 = "U00096.2,1,+,80M160S,60,0"
    # alignment 3
    seq3 = "acttaacgctgctcgtagcgtttaaacaccagttcgccattgctggaggaatcttcatcaaagaagtaaccttcgctatt"
    pos3 = 5680
    cig3 = "160S80M"
    sap3 = "U00096.2,5681,+,160S80M,60,0"

    opt = BwaMemOptions(query_coord_as_primary=query_coord_as_primary, soft_clip_supplementary=True)
    bwa = BwaMem(prefix=e_coli_k12_fasta)
    query = FastxRecord(name="query", sequence=seq1 + seq2 + seq3)

    list_of_recs = bwa.align(opt=opt, queries=[query])
    assert len(list_of_recs) == 1
    assert len(list_of_recs[0]) == 3
    rec1 = list_of_recs[0][0]
    rec2 = list_of_recs[0][1]
    rec3 = list_of_recs[0][2]
    sa1 = f"{sap2};{sap3};"
    sa2 = f"{sap1};{sap3};"
    sa3 = f"{sap1};{sap2};"

    if query_coord_as_primary:  # query coordinate
        # swap rec1 and rec2, rec2 is primary
        tmp = rec1
        rec1 = rec2
        rec2 = tmp
        sa3 = f"{sap2};{sap1};"  # primary order changed

    assert rec1.reference_start == pos1
    assert rec2.reference_start == pos2
    assert rec3.reference_start == pos3

    assert rec1.cigarstring == cig1
    assert rec2.cigarstring == cig2
    assert rec3.cigarstring == cig3

    assert rec1.get_tag("SA") == sa1
    assert rec2.get_tag("SA") == sa2
    assert rec3.get_tag("SA") == sa3


def test_bwamem_align_hard_clip_supplementary(e_coli_k12_fasta: Path) -> None:
    # alignments below are defined based on (1) alignment score, then (2) genomic order
    # alignment 1
    seq1 = "agcttttcattctgactgcaacgggcaatatgtctctgtgtggattaaaaaaagagtgtctgatagcagcttctgaactg"
    pos1 = 1759
    cig1 = "79S81M"
    sap1 = "U00096.2,1760,+,79S81M,60,0"
    # alignment 2
    seq2 = "cgttggcggtgcgctgctggagcaactgaagcgtcagcaaagctggctgaagaataaacatatcgacttacgtgtctgcg"
    pos2 = 0
    cig2 = "80M80H"
    sap2 = "U00096.2,1,+,80M80S,60,0"

    opt = BwaMemOptions(soft_clip_supplementary=False)
    bwa = BwaMem(prefix=e_coli_k12_fasta)
    query = FastxRecord(name="query", sequence=seq1 + seq2)
    list_of_recs = bwa.align(opt=opt, queries=[query])
    assert len(list_of_recs) == 1
    assert len(list_of_recs[0]) == 2
    rec1 = list_of_recs[0][0]
    rec2 = list_of_recs[0][1]

    assert rec1.reference_start == pos1
    assert rec2.reference_start == pos2

    assert rec1.cigarstring == cig1
    assert rec2.cigarstring == cig2

    assert rec1.get_tag("SA") == f"{sap2};"
    assert rec2.get_tag("SA") == f"{sap1};"
