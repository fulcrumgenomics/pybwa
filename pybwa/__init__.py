import os
import sysconfig

from pybwa.libbwaaln import *  # noqa: F403
from pybwa.libbwaindex import *  # noqa: F403
from pybwa.libbwamem import *  # noqa: F403


def _get_include() -> list[str]:
    """return a list of include directories."""
    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__)))

    # Header files may be stored in different relative locations
    # depending on installation mode (e.g., `python setup.py install`,
    # `python setup.py develop`). The first entry in each list is
    # where develop-mode headers can be found.
    #
    pybwa_possibilities = [
        os.path.join(dirname, "..", "bwa"),
        os.path.join(dirname, "include", "bwa"),
    ]

    includes = [dirname]
    for header_locations in [pybwa_possibilities]:
        for header_location in header_locations:
            if os.path.exists(header_location):
                includes.append(os.path.abspath(header_location))
                break

    return includes


def _get_defines() -> list[str]:
    """return a list of defined compilation parameters."""
    return []


def _get_libraries() -> list[str]:
    """return a list of libraries to link against."""
    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    pybwa_libs = ["libbwaaln", "libbwaindex", "libbwamem"]

    so = sysconfig.get_config_var("EXT_SUFFIX")
    return [os.path.join(dirname, x + so) for x in pybwa_libs]
