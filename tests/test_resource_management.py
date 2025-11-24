"""Test resource management and memory safety for pybwa."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from pybwa import BwaAln
from pybwa import BwaAlnOptions
from pybwa import BwaIndex
from pybwa import BwaMem
from pybwa import BwaMemOptions


class TestResourceCleanup:
    """Test that resources are properly cleaned up."""

    def test_index_cleanup(self, test_fa: Path) -> None:
        """Test that BwaIndex properly cleans up resources."""
        # Create and destroy multiple index instances
        for _ in range(10):
            idx = BwaIndex(prefix=test_fa)
            del idx
        # If there were memory leaks, this would accumulate

    def test_bwa_mem_cleanup(self, test_fa: Path) -> None:
        """Test that BwaMem properly cleans up resources."""
        for _ in range(10):
            aligner = BwaMem(prefix=test_fa)
            # Perform some alignments
            aligner.align(["ACGTACGT"])
            del aligner

    def test_bwa_aln_cleanup(self, test_fa: Path) -> None:
        """Test that BwaAln properly cleans up resources."""
        for _ in range(10):
            aligner = BwaAln(prefix=test_fa)
            aligner.align(["ACGTACGT"])
            del aligner

    def test_options_cleanup(self) -> None:
        """Test that options objects properly clean up memory."""
        for _ in range(100):
            opts_mem = BwaMemOptions()
            opts_mem.finalize()
            del opts_mem

            opts_aln = BwaAlnOptions()
            del opts_aln


class TestReuseability:
    """Test that objects can be reused safely."""

    def test_reuse_mem_aligner(self, bwa_mem_aligner: BwaMem) -> None:
        """Test that BwaMem can be used multiple times."""
        queries = ["ACGTACGT", "TGCATGCA", "GGGGCCCC"]

        for _ in range(5):
            results = bwa_mem_aligner.align(queries)
            assert len(results) == 3

    def test_reuse_aln_aligner(self, bwa_aln_aligner: BwaAln) -> None:
        """Test that BwaAln can be used multiple times."""
        queries = ["ACGTACGT", "TGCATGCA", "GGGGCCCC"]

        for _ in range(5):
            results = bwa_aln_aligner.align(queries)
            assert len(results) == 3

    def test_share_index_between_aligners(self, test_fa: Path) -> None:
        """Test that an index can be shared between multiple aligners."""
        # This is actually created fresh each time, but tests the pattern
        idx = BwaIndex(prefix=test_fa)

        mem_aligner = BwaMem(index=idx)
        aln_aligner = BwaAln(index=idx)

        mem_results = mem_aligner.align(["ACGTACGT"])
        aln_results = aln_aligner.align(["ACGTACGT"])

        assert len(mem_results) == 1
        assert len(aln_results) == 1


class TestLargeBatchProcessing:
    """Test handling of large batches of queries."""

    def test_many_queries_mem(self, bwa_mem_aligner: BwaMem) -> None:
        """Test BwaMem with many queries."""
        # Generate 1000 short queries
        queries = [f"ACGTACGT{'ACGT'[i % 4]}" * 10 for i in range(1000)]

        results = bwa_mem_aligner.align(queries)
        assert len(results) == 1000

    def test_many_queries_aln(self, bwa_aln_aligner: BwaAln) -> None:
        """Test BwaAln with many queries."""
        # Generate 1000 short queries
        queries = [f"ACGTACGT{'ACGT'[i % 4]}" * 10 for i in range(1000)]

        results = bwa_aln_aligner.align(queries)
        assert len(results) == 1000

    def test_chunking_behavior(self, bwa_mem_aligner: BwaMem) -> None:
        """Test that chunking works correctly with different chunk sizes."""
        queries = ["ACGTACGT" * 20 for _ in range(100)]

        # Test with different chunk sizes
        for chunk_size in [1000, 10000, 100000]:
            opts = BwaMemOptions(chunk_size=chunk_size)
            opts = opts.finalize()
            results = bwa_mem_aligner.align(queries, opt=opts)
            assert len(results) == 100


class TestConcurrentUsage:
    """Test behavior with concurrent/parallel usage patterns."""

    def test_multiple_aligners_same_index(self, test_fa: Path) -> None:
        """Test using multiple aligner instances with the same index file in parallel."""
        aligners = [BwaMem(prefix=test_fa) for _ in range(5)]

        with ThreadPoolExecutor(max_workers=5) as executor:
            results = list(executor.map(lambda a: a.align(["ACGTACGT"]), aligners))
        assert all(len(r) == 1 for r in results)

        del aligners

    def test_different_options_same_aligner(self, bwa_mem_aligner: BwaMem) -> None:
        """Test using the same aligner with different options."""
        query = ["ACGTACGTACGTACGT"]

        opts1 = BwaMemOptions(min_seed_len=10).finalize()
        opts2 = BwaMemOptions(min_seed_len=15).finalize()
        opts3 = BwaMemOptions(min_seed_len=20).finalize()

        # Use same aligner with different options
        result1 = bwa_mem_aligner.align(query, opt=opts1)
        result2 = bwa_mem_aligner.align(query, opt=opts2)
        result3 = bwa_mem_aligner.align(query, opt=opts3)

        assert len(result1) == 1
        assert len(result2) == 1
        assert len(result3) == 1


class TestOptionsImmutability:
    """Test that finalized options cannot be modified."""

    def test_mem_options_immutable_after_finalize(self) -> None:
        """
        Test that most BwaMemOptions cannot be modified after finalization.

        Note: threads and chunk_size can be changed after finalization (by design).
        """
        opts = BwaMemOptions(min_seed_len=19)
        opts = opts.finalize()

        # Most setters should raise AttributeError
        with pytest.raises(AttributeError):
            opts.min_seed_len = 20

        with pytest.raises(AttributeError):
            opts.match_score = 2

        # But threads can be changed (intentional design for flexibility)
        opts.threads = 4  # Should not raise

    def test_mem_options_copy_finalize(self) -> None:
        """Test that finalize(copy=True) creates an independent copy."""
        opts1 = BwaMemOptions(min_seed_len=19, threads=2)
        opts2 = opts1.finalize(copy=True)

        assert opts2.finalized
        assert opts2.min_seed_len == 19
        assert opts2.threads == 2

        # Original should not be finalized
        assert not opts1.finalized

        # Can still modify original
        opts1.threads = 4
        assert opts1.threads == 4
        assert opts2.threads == 2
