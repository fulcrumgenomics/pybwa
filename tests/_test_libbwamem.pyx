# cython: language_level=3

from pybwa.libbwamem cimport BwaMemOptions
from pybwa.libbwamem cimport mem_opt_t
from pybwa.libbwamem cimport mem_opt_init
from libc.stdlib cimport free
import pytest

cdef _test_bwamem_options_default():
    cdef mem_opt_t *bwa_mem_opt
    cdef mem_opt_t *pybwa_mem_opt
    bwa_mem_opt = mem_opt_init()
    options = BwaMemOptions()
    pybwa_mem_opt = options._options

    # Check member by member for the mem_opt_t struct
    assert pybwa_mem_opt.a == bwa_mem_opt.a
    assert pybwa_mem_opt.b == bwa_mem_opt.b
    assert pybwa_mem_opt.o_del == bwa_mem_opt.o_del
    assert pybwa_mem_opt.e_del == bwa_mem_opt.e_del
    assert pybwa_mem_opt.o_ins == bwa_mem_opt.o_ins
    assert pybwa_mem_opt.e_ins == bwa_mem_opt.e_ins
    assert pybwa_mem_opt.pen_unpaired == bwa_mem_opt.pen_unpaired
    assert pybwa_mem_opt.pen_clip5 == bwa_mem_opt.pen_clip5
    assert pybwa_mem_opt.pen_clip3 == bwa_mem_opt.pen_clip3
    assert pybwa_mem_opt.w == bwa_mem_opt.w
    assert pybwa_mem_opt.zdrop == bwa_mem_opt.zdrop
    assert pybwa_mem_opt.max_mem_intv == bwa_mem_opt.max_mem_intv
    assert pybwa_mem_opt.T == bwa_mem_opt.T
    assert pybwa_mem_opt.flag == bwa_mem_opt.flag
    assert pybwa_mem_opt.min_seed_len == bwa_mem_opt.min_seed_len
    assert pybwa_mem_opt.min_chain_weight == bwa_mem_opt.min_chain_weight
    assert pybwa_mem_opt.max_chain_extend == bwa_mem_opt.max_chain_extend
    assert pybwa_mem_opt.split_factor == bwa_mem_opt.split_factor
    assert pybwa_mem_opt.split_width == bwa_mem_opt.split_width
    assert pybwa_mem_opt.max_occ == bwa_mem_opt.max_occ
    assert pybwa_mem_opt.max_chain_gap == bwa_mem_opt.max_chain_gap
    assert pybwa_mem_opt.n_threads == bwa_mem_opt.n_threads
    assert pybwa_mem_opt.chunk_size == bwa_mem_opt.chunk_size
    assert pybwa_mem_opt.mask_level == bwa_mem_opt.mask_level
    assert pybwa_mem_opt.drop_ratio == bwa_mem_opt.drop_ratio
    assert pybwa_mem_opt.XA_drop_ratio == bwa_mem_opt.XA_drop_ratio
    assert pybwa_mem_opt.mask_level_redun == bwa_mem_opt.mask_level_redun
    assert pybwa_mem_opt.mapQ_coef_len == bwa_mem_opt.mapQ_coef_len
    assert pybwa_mem_opt.max_ins == bwa_mem_opt.max_ins
    assert pybwa_mem_opt.max_matesw == bwa_mem_opt.max_matesw
    assert pybwa_mem_opt.max_XA_hits == bwa_mem_opt.max_XA_hits
    assert pybwa_mem_opt.max_XA_hits_alt == bwa_mem_opt.max_XA_hits_alt

    free(bwa_mem_opt)


def test_bwamem_options_default() -> None:
    """Tests that the defaults are synced between bwa and pybwa."""
    _test_bwamem_options_default()
