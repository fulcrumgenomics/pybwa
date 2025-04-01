# cython: language_level=3

from pybwa.libbwaaln cimport BwaAlnOptions
from pybwa.libbwaaln cimport gap_opt_t
from pybwa.libbwaaln cimport gap_init_opt
from libc.stdlib cimport free
import pytest

cdef _test_bwaaln_options_default():
    cdef gap_opt_t *bwa_gap_opt
    cdef gap_opt_t *pybwa_gap_opt
    bwa_gap_opt = gap_init_opt()
    options = BwaAlnOptions()
    pybwa_gap_opt = options.gap_opt()

    # Check member by member for the gap_opt_t struct
    assert pybwa_gap_opt.s_mm == bwa_gap_opt.s_mm
    assert pybwa_gap_opt.s_gape == bwa_gap_opt.s_gape
    assert pybwa_gap_opt.mode == bwa_gap_opt.mode
    assert pybwa_gap_opt.indel_end_skip == bwa_gap_opt.indel_end_skip
    assert pybwa_gap_opt.max_del_occ == bwa_gap_opt.max_del_occ
    assert pybwa_gap_opt.max_entries == bwa_gap_opt.max_entries
    assert pybwa_gap_opt.fnr == bwa_gap_opt.fnr
    assert pybwa_gap_opt.max_diff == bwa_gap_opt.max_diff
    assert pybwa_gap_opt.max_gapo == bwa_gap_opt.max_gapo
    assert pybwa_gap_opt.max_gape == bwa_gap_opt.max_gape
    assert pybwa_gap_opt.max_seed_diff == bwa_gap_opt.max_seed_diff
    assert pybwa_gap_opt.seed_len == bwa_gap_opt.seed_len
    assert pybwa_gap_opt.n_threads == bwa_gap_opt.n_threads
    assert pybwa_gap_opt.max_top2 == bwa_gap_opt.max_top2
    assert pybwa_gap_opt.trim_qual == bwa_gap_opt.trim_qual

    free(bwa_gap_opt)


def test_bwaaln_options_default() -> None:
    """Tests that the defaults are synced between bwa and pybwa."""
    _test_bwaaln_options_default()
