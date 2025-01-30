from pathlib import Path

import pytest

from pybwa import BwaIndex


def test_bwa_index_build(e_coli_k12_fasta: Path, tmp_path_factory: pytest.TempPathFactory) -> None:
    tmp_dir = Path(str(tmp_path_factory.mktemp("test_bwa_index_build")))
    prefix = tmp_dir / e_coli_k12_fasta.name

    # Build the index
    BwaIndex.index(fasta=e_coli_k12_fasta, prefix=prefix)

    # Check the files exist
    for suffix in [".bwt", ".sa", ".ann", ".pac"]:
        assert prefix.with_suffix(prefix.suffix + suffix).exists()
    assert prefix.with_suffix(".dict").exists()

    # Load it
    BwaIndex(prefix=prefix)


def test_bwa_index_load(e_coli_k12_fasta: Path) -> None:
    BwaIndex(prefix=e_coli_k12_fasta)
