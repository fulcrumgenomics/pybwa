import sysconfig
from pathlib import Path


def import_test_lib(libname: str) -> None:
    """
    Imports the cython test library.

    The cython test library must compile successfully.

    Args:
        libname: the name of the cython test library to import.

    Raises:
        ImportError: if the cython test library could not be compiled during import due to one of
        the following underlying exceptions: `distutils.compilers.C.errors.CompileError` or
        `Cython.Compiler.Errors.CompileError`.
    """
    so_ext = sysconfig.get_config_var("EXT_SUFFIX")
    so = Path("tests") / f"_test_{libname}{so_ext}"
    paths = [
        Path("tests") / f"_test_{libname}.c",
        Path("tests") / f"_test_{libname}.pyxbldc",
        so,
    ]
    for path in paths:
        path.unlink(missing_ok=True)

    try:
        import ctypes
        import importlib

        import pyximport

        pyximport.install(build_in_temp=False, inplace=True)
        importlib.import_module(f"_test_{libname}")
        ctypes.cdll.LoadLibrary(f"{so}")
    except ImportError as ex:
        patterns = [
            "distutils.compilers.C.errors.CompileError",
            "Cython.Compiler.Errors.CompileError",
        ]
        for pattern in patterns:
            if pattern in f"{ex}":
                raise Exception(f"Cannot compile {libname}:\n\t{ex}") from ex
        else:
            raise ex
