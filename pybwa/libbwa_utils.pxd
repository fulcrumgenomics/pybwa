# cython: language_level=3
"""Shared utility functions for error checking across pybwa modules."""

cdef void* check_alloc(void* ptr, str msg) except NULL
cdef void* check_not_null(void* ptr, str msg) except NULL
