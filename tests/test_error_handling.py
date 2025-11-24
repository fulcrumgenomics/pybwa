"""Test error handling and edge cases for pybwa."""

import tempfile
from pathlib import Path

import pytest
from pysam import FastxRecord

from pybwa import BwaAln
from pybwa import BwaAlnOptions
from pybwa import BwaIndex
from pybwa import BwaIndexBuildMethod
from pybwa import BwaMem
from pybwa import BwaMemMode
from pybwa import BwaMemOptions


class TestBwaIndexErrorHandling:
    """Test error handling in BwaIndex."""

    def test_index_not_found(self):
        """Test that FileNotFoundError is raised when index files are missing."""
        with pytest.raises(FileNotFoundError, match="could not locate the index file"):
            BwaIndex(prefix="/nonexistent/path/to/index")

    def test_index_missing_dict_file(self, tmp_path: Path, test_fa: Path):
        """Test that FileNotFoundError is raised when .dict file is missing."""
        # Create a temporary index without dict file
        prefix = tmp_path / "test"
        BwaIndex.index(fasta=test_fa, prefix=prefix)

        # Remove the dict file
        dict_file = prefix.with_suffix(".dict")
        if dict_file.exists():
            dict_file.unlink()

        # Try to load the index
        with pytest.raises(FileNotFoundError, match="could not locate the sequence dictionary"):
            BwaIndex(prefix=prefix)


class TestBwaMemErrorHandling:
    """Test error handling in BwaMem."""

    def test_no_prefix_or_index(self):
        """Test that ValueError is raised when neither prefix nor index is provided."""
        with pytest.raises(ValueError, match="Either prefix or index must be given"):
            BwaMem()

    def test_invalid_prefix(self):
        """Test that FileNotFoundError is raised for non-existent prefix."""
        with pytest.raises(FileNotFoundError, match="Index prefix does not exist"):
            BwaMem(prefix="/nonexistent/path")

    def test_align_empty_queries(self, bwa_mem_aligner: BwaMem):
        """Test alignment with empty query list."""
        result = bwa_mem_aligner.align([])
        assert result == []

    def test_context_manager(self, test_fa: Path):
        """Test BwaMem works as a context manager."""
        with BwaMem(prefix=test_fa) as aligner:
            assert aligner is not None
            # Should be able to use the aligner
            results = aligner.align(["ACGT"])
            assert len(results) == 1


class TestBwaAlnErrorHandling:
    """Test error handling in BwaAln."""

    def test_no_prefix_or_index(self):
        """Test that ValueError is raised when neither prefix nor index is provided."""
        with pytest.raises(ValueError, match="Either prefix or index must be given"):
            BwaAln()

    def test_invalid_prefix(self):
        """Test that FileNotFoundError is raised for non-existent prefix."""
        with pytest.raises(FileNotFoundError, match="Index prefix does not exist"):
            BwaAln(prefix="/nonexistent/path")

    def test_align_none_queries(self, bwa_aln_aligner: BwaAln):
        """Test alignment with None queries."""
        result = bwa_aln_aligner.align(queries=None)
        assert result == []

    def test_align_empty_queries(self, bwa_aln_aligner: BwaAln):
        """Test alignment with empty query list."""
        result = bwa_aln_aligner.align([])
        assert result == []

    def test_context_manager(self, test_fa: Path):
        """Test BwaAln works as a context manager."""
        with BwaAln(prefix=test_fa) as aligner:
            assert aligner is not None
            # Should be able to use the aligner
            results = aligner.align(["ACGT"])
            assert len(results) == 1


class TestBwaMemOptionsValidation:
    """Test input validation for BwaMemOptions."""

    def test_invalid_min_seed_len(self):
        """Test that validation catches invalid min_seed_len."""
        opts = BwaMemOptions(min_seed_len=0)
        with pytest.raises(ValueError, match="min_seed_len must be >= 1"):
            opts.validate()

    def test_invalid_band_width(self):
        """Test that validation catches negative band_width."""
        opts = BwaMemOptions(band_width=-1)
        with pytest.raises(ValueError, match="band_width must be >= 0"):
            opts.validate()

    def test_invalid_threads(self):
        """Test that validation catches invalid threads."""
        opts = BwaMemOptions(threads=0)
        with pytest.raises(ValueError, match="threads must be >= 1"):
            opts.validate()

    def test_invalid_chunk_size(self):
        """Test that validation catches negative chunk_size."""
        opts = BwaMemOptions(chunk_size=-1)
        with pytest.raises(ValueError, match="chunk_size must be >= 0"):
            opts.validate()

    def test_invalid_max_occurrences(self):
        """Test that validation catches negative max_occurrences."""
        opts = BwaMemOptions(max_occurrences=-1)
        with pytest.raises(ValueError, match="max_occurrences must be >= 0"):
            opts.validate()

    def test_validation_called_in_finalize(self):
        """Test that finalize calls validate."""
        opts = BwaMemOptions(min_seed_len=0)
        with pytest.raises(ValueError, match="min_seed_len must be >= 1"):
            opts.finalize()

    def test_valid_options(self):
        """Test that valid options don't raise errors."""
        opts = BwaMemOptions(min_seed_len=19, threads=4)
        opts.validate()  # Should not raise
        opts = opts.finalize()
        assert opts.finalized


class TestBwaAlnOptionsValidation:
    """Test input validation for BwaAlnOptions."""

    def test_invalid_seed_length(self):
        """Test that validation catches invalid seed_length."""
        opts = BwaAlnOptions(seed_length=0)
        with pytest.raises(ValueError, match="seed_length must be >= 1"):
            opts.validate()

    def test_invalid_threads(self):
        """Test that validation catches invalid threads."""
        opts = BwaAlnOptions(threads=0)
        with pytest.raises(ValueError, match="threads must be >= 1"):
            opts.validate()

    def test_invalid_max_hits(self):
        """Test that validation catches negative max_hits."""
        opts = BwaAlnOptions(max_hits=-1)
        with pytest.raises(ValueError, match="max_hits must be >= 0"):
            opts.validate()

    def test_invalid_max_entries(self):
        """Test that validation catches invalid max_entries."""
        opts = BwaAlnOptions(max_entries=0)
        with pytest.raises(ValueError, match="max_entries must be >= 1"):
            opts.validate()

    def test_valid_options(self):
        """Test that valid options don't raise errors."""
        opts = BwaAlnOptions(seed_length=32, threads=4)
        opts.validate()  # Should not raise


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_short_sequence(self, bwa_mem_aligner: BwaMem):
        """Test alignment of very short sequences."""
        result = bwa_mem_aligner.align(["A"])
        assert len(result) == 1
        # Very short sequences may not align, check that we get a result

    def test_sequence_with_all_ns(self, bwa_mem_aligner: BwaMem):
        """Test alignment of sequence with all N bases."""
        result = bwa_mem_aligner.align(["NNNNNNNNNN"])
        assert len(result) == 1
        # All-N sequences should return unmapped alignments

    def test_string_queries(self, bwa_mem_aligner: BwaMem):
        """Test that string queries are properly converted to FastxRecord."""
        result = bwa_mem_aligner.align(["ACGTACGT", "TGCATGCA"])
        assert len(result) == 2
        assert all(len(hits) >= 0 for hits in result)  # May or may not align

    def test_multiple_modes(self):
        """Test different BwaMemMode settings."""
        for mode in [BwaMemMode.PACBIO, BwaMemMode.ONT2D, BwaMemMode.INTRACTG]:
            opts = BwaMemOptions(mode=mode)
            opts = opts.finalize()
            assert opts.finalized

    def test_options_cannot_be_modified_after_finalize(self):
        """Test that most options cannot be modified after finalization.

        Note: threads and chunk_size are intentionally allowed to change after finalization.
        """
        opts = BwaMemOptions()
        opts = opts.finalize()

        # Most options should raise AttributeError when modified after finalization
        with pytest.raises(AttributeError, match="can't set attribute"):
            opts.min_seed_len = 20

        # But threads can be changed (intentional design)
        opts.threads = 4  # Should not raise
