# cython: language_level=3
"""Shared utility functions for error checking across pybwa modules."""


cdef void* check_alloc(void* ptr, str msg) except NULL:
    """Check if memory allocation succeeded, raise MemoryError if NULL.

    This helper consolidates the common pattern of checking calloc/malloc return values
    and raising appropriate errors with descriptive messages.

    Args:
        ptr: The pointer returned from calloc/malloc
        msg: Error message to include in the exception

    Returns:
        The pointer if non-NULL

    Raises:
        MemoryError: If ptr is NULL
    """
    if ptr == NULL:
        raise MemoryError(msg)
    return ptr


cdef void* check_not_null(void* ptr, str msg) except NULL:
    """Check if pointer is not NULL, raise RuntimeError if it is.

    Used for checking return values from C library functions that aren't memory
    allocations (e.g., bwa_idx_load, sam_hdr_parse).

    Args:
        ptr: The pointer to check
        msg: Error message to include in the exception

    Returns:
        The pointer if non-NULL

    Raises:
        RuntimeError: If ptr is NULL
    """
    if ptr == NULL:
        raise RuntimeError(msg)
    return ptr


###########################################################
# Below is code to facilitate testing
###########################################################


def _check_alloc_with_null() -> None:
    """Test wrapper: calls check_alloc with NULL to exercise the MemoryError path."""
    check_alloc(NULL, "test allocation failure")


def _check_not_null_with_null() -> None:
    """Test wrapper: calls check_not_null with NULL to exercise the RuntimeError path."""
    check_not_null(NULL, "test null pointer")
