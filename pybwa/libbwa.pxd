# cython: language_level=3


cdef extern from "libbwa_utils.h":
    int set_bwa_c_verbosity(int level)

