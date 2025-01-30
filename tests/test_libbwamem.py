from pathlib import Path

from fgpyo.sequence import reverse_complement
from pysam import FastxRecord

from pybwa.libbwamem import BwaMem
from pybwa.libbwamem import BwaMemOptions


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
