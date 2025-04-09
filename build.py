import dataclasses
import logging
import multiprocessing
import os
import platform
import re
import shutil
import subprocess
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import List

import pysam
from Cython.Build import cythonize
from Cython.Distutils.build_ext import new_build_ext as cython_build_ext
from setuptools import Distribution
from setuptools import Extension


@dataclass
class HtslibConfig:
    """Stores library config as inferred from htslib's config.h."""

    has_libbz2: bool = False
    has_libcurl: bool = False
    has_libdeflate: bool = False
    has_liblzma: bool = False

    @classmethod
    def from_htslib_config_h(cls, config_h: Path) -> "HtslibConfig":
        """Builds an instance from htslib's config.h."""
        this: HtslibConfig = HtslibConfig()
        with config_h.open("r") as reader:
            for line in reader:
                if not line.startswith("#define"):
                    continue
                key, value = re.match(r"#define (\S+)\s+(\S+)", line).groups()
                if key == "HAVE_LIBBZ2":
                    this = dataclasses.replace(this, has_libbz2=(value == "1"))
                elif key == "HAVE_LIBCURL":
                    this = dataclasses.replace(this, has_libcurl=(value == "1"))
                elif key == "HAVE_LIBDEFLATE":
                    this = dataclasses.replace(this, has_libdeflate=(value == "1"))
                elif key == "HAVE_LIBLZMA":
                    this = dataclasses.replace(this, has_liblzma=(value == "1"))
        return this


def strtobool(value: str) -> bool:
    """Equivalent to distutils strtobool."""
    value = value.lower()
    trues = {"y", "yes", "t", "true", "on", "1"}
    falses = {"n", "no", "f", "false", "off", "0"}
    if value in trues:
        return True
    elif value in falses:
        return False
    raise ValueError(f"'{value}' is not a valid bool value")


@contextmanager
def changedir(path: str) -> Any:
    """Changes the directory before, and moves back to the original directory after."""
    save_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(save_dir)


USE_GIT: bool = shutil.which("git") is not None and Path(".git").exists() and Path(".git").is_dir()
IS_DARWIN = platform.system() == "Darwin"


def run_command(command: str, fail_on_error: bool = True) -> int:
    """Runs a given command, return the return code or failing on error if desired."""
    retcode = subprocess.call(
        f"{command}", stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True
    )
    if fail_on_error and retcode != 0:  # run again for debugging
        subprocess.call(f"{command}", shell=True)
        raise RuntimeError(f"Failed to run command {command}")
    return retcode


def compile_htslib(logger: logging.Logger) -> None:
    """Complies htslib prior to entering the context."""
    with changedir("htslib"):
        cflags = "CFLAGS='-fpic -fvisibility=hidden -g -Wall -O2"
        if IS_DARWIN:
            cflags += " -mmacosx-version-min=11.0"
        cflags += "'"

        # regenerate GNU build files
        logger.info("Regenerating htslib build files...")
        run_command(command="autoreconf -i")

        # run configure, trying to disable bz2 and lzma when it fails
        logger.info("Configuring htslib...")
        configure_options = [
            "",  # all features
            "--disable-bz2",  # no bz2
            "--disable-lzma",  # no lzma
            "--disable-bz2 --disable-lzma",  # no bz2 and no lzma
        ]
        for options in configure_options:
            if options != "":
                logger.warning(f"... retrying with configure options [{options}]")
            command = f"./configure {options} {cflags}"
            if run_command(command=command, fail_on_error=False) == 0:
                break
        else:  # run one last time to failure
            command = f"./configure {cflags}"
            run_command(command=command)

        # build it
        logger.info("Building htslib...")
        run_command(command="make -j")


@contextmanager
def with_patches(logger: logging.Logger) -> Any:
    """Applies patches to bwa, then cleans up after exiting the context."""
    patches = sorted([
        os.path.abspath(patch)
        for patch in Path("patches").iterdir()
        if patch.is_file() and patch.suffix == ".patch"
    ])
    with changedir("bwa"):
        logger.info("Patching bwa...")
        for patch in patches:
            if USE_GIT:
                retcode = subprocess.call(f"git apply --whitespace=nowarn {patch}", shell=True)
            else:
                retcode = subprocess.call(f"patch -p1 < {patch}", shell=True)
            if retcode != 0:
                raise RuntimeError(f"Failed to apply patch {patch}")
    try:
        yield
    finally:
        if USE_GIT:
            logger.info("Restoring bwa...")
            commands = ["git submodule deinit -f bwa", "git submodule update --init --recursive"]
            for command in commands:
                retcode = subprocess.call(command, shell=True)
                if retcode != 0:
                    raise RuntimeError(f"Failed to reset bwa submodule: {command}")


compiler_directives = {
    "language_level": "3",
    "embedsignature": True,
}
SOURCE_DIR = Path("pybwa")
BUILD_DIR = Path("cython_build")
compile_args: list[str] = []
link_args: list[str] = []
include_dirs = ["htslib", "bwa", "pybwa", os.path.dirname(pysam.__file__)]
libraries = ["m", "z", "pthread"]
if platform.system() == "Linux":
    libraries.append("rt")
library_dirs = ["pybwa", "bwa", "htslib"]
extra_objects = ["htslib/libhts.a"]
define_macros = [
    ("HAVE_PTHREAD", None),
    ("USE_MALLOC_WRAPPERS", None),
    ("USE_HTSLIB", "1"),
]
h_files: list[str] = []
c_files: list[str] = []

exclude_files = {
    "pybwa": {"libbwaaln.c", "libbwaindex.c", "libbwamem.c"},
    "bwa": {"example.c", "main.c", "kstring.c", "kstring.h"},
}
for root_dir in library_dirs:
    if root_dir == "htslib":  # these are added via extra_objects above
        continue
    exc: set[str] = exclude_files.get(root_dir, set())
    h_files.extend(
        str(x) for x in Path(root_dir).rglob("*.h") if "tests/" not in x.parts and x.name not in exc
    )
    c_files.extend(
        str(x) for x in Path(root_dir).rglob("*.c") if "tests/" not in x.parts and x.name not in exc
    )

# Check if we should build with linetracing for coverage
build_with_coverage = os.environ.get("CYTHON_TRACE", "false").lower() in ("1", "true", "'true'")
if build_with_coverage:
    compiler_directives["linetrace"] = True
    compiler_directives["emit_code_comments"] = True
    define_macros.extend([
        ("CYTHON_TRACE", "1"),
        ("CYTHON_TRACE_NOGIL", "1"),
        ("DCYTHON_USE_SYS_MONITORING", "0"),
    ])
    BUILD_DIR = Path(".")  # the compiled .c files need to be next to the .pyx files for coverage

if platform.system() != "Windows":
    compile_args.extend([
        "-Wno-unused-result",
        "-Wno-unreachable-code",
        "-Wno-single-bit-bitfield-constant-conversion",
        "-Wno-deprecated-declarations",
        "-Wno-unused",
        "-Wno-strict-prototypes",
        "-Wno-sign-compare",
        "-Wno-error=declaration-after-statement",
        "-Wno-implicit-function-declaration",
        "-Wno-macro-redefined",
    ])


def cythonize_helper(extension_modules: List[Extension]) -> Any:
    """Cythonize all Python extensions."""
    return cythonize(
        module_list=extension_modules,
        # Don"t build in source tree (this leaves behind .c files)
        build_dir=BUILD_DIR,
        # Don"t generate an .html output file. Would contain source.
        annotate=False,
        # Parallelize our build
        nthreads=multiprocessing.cpu_count() * 2,
        # Compiler directives (e.g. language, or line tracing for coverage)
        compiler_directives=compiler_directives,
        # (Optional) Always rebuild, even if files untouched
        force=True,
    )


CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Developers",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Natural Language :: English",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Programming Language :: Python :: 3.13",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Software Development :: Documentation",
    "Topic :: Software Development :: Libraries :: Python Modules",
    "Typing :: Typed",
]


def build() -> None:
    """The main build function for pybwa."""
    # set up logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s",
    )
    logger = logging.getLogger("pybwa build")
    logger.setLevel("INFO")

    # compile htslib
    compile_htslib(logger=logger)

    # apply patches to bwa, then revert them after
    with with_patches(logger=logger):
        # get the configuration from htslib
        logger.info("Setting pybwa configuration...")
        config = HtslibConfig.from_htslib_config_h(Path("htslib/config.h"))
        if config.has_libbz2:
            logger.info("... building **with** bz2 support")
            libraries.append("bz2")
        else:
            logger.warning("... building **without** bz2 support")
        if config.has_libcurl:
            logger.info("... building **with** curl support")
            libraries.append("curl")
        else:
            logger.warning("... building **without** curl support")
        if config.has_libdeflate:
            logger.info("... building **with** libdeflate support")
            libraries.append("deflate")
        else:
            logger.warning("... building **without** libdeflate support")
        if config.has_liblzma:
            logger.info("... building **with** lzma support")
            libraries.append("lzma")
        else:
            logger.warning("... building **without** lzma support")

        # Define the extension modules
        libbwaindex_module = Extension(
            name="pybwa.libbwaindex",
            sources=["pybwa/libbwaindex.pyx"] + c_files,
            depends=h_files,
            extra_compile_args=compile_args,
            extra_link_args=link_args,
            extra_objects=extra_objects,
            include_dirs=include_dirs,
            language="c",
            libraries=libraries,
            library_dirs=library_dirs,
            define_macros=define_macros,
        )

        libbwaaln_module = Extension(
            name="pybwa.libbwaaln",
            sources=["pybwa/libbwaaln.pyx"] + c_files,
            depends=h_files,
            extra_compile_args=compile_args,
            extra_link_args=link_args,
            extra_objects=extra_objects,
            include_dirs=include_dirs,
            language="c",
            libraries=libraries,
            library_dirs=library_dirs,
            define_macros=define_macros,
        )

        libbwamem_module = Extension(
            name="pybwa.libbwamem",
            sources=["pybwa/libbwamem.pyx"] + c_files,
            depends=h_files,
            extra_compile_args=compile_args,
            extra_link_args=link_args,
            extra_objects=extra_objects,
            include_dirs=include_dirs,
            language="c",
            libraries=libraries,
            library_dirs=library_dirs,
            define_macros=define_macros,
        )

        logger.info("Cythonize extensions...")
        # Collect and cythonize all files
        extension_modules = cythonize_helper([
            libbwaindex_module,
            libbwaaln_module,
            libbwamem_module,
        ])

        packages = ["pybwa", "pybwa.include.bwa", "pybwa.include.patches", "pybwa.include.htslib"]
        package_dir = {
            "pybwa": "pybwa",
            "pybwa.include.bwa": "bwa",
            "pybwa.include.patches": "patches",
            "pybwa.include.htslib": "htslib",
        }

        # Use Setuptools to collect files
        distribution = Distribution({
            "name": "pybwa",
            "version": "0.0.1",
            "description": "Python bindings for BWA",
            "long_description": __doc__,
            "long_description_content_type": "text/x-rst",
            "author": "Nils Homer",
            "author_email": "nils@fulcrumgenomics.com",
            "license": "MIT",
            "platforms": ["POSIX", "UNIX", "MacOS"],
            "classifiers": CLASSIFIERS,
            "url": "https://github.com/fulcrumgenomics/pybwa",
            "packages": packages,
            "package_dir": package_dir,
            "package_data": {
                "": ["*.pxd", "*.h", "*.c", "py.typed", "*.pyi", "*.patch", "**/*.h", "**/*.c"],
            },
            "ext_modules": extension_modules,
            "cmdclass": {
                "build_ext": cython_build_ext,
            },
            "zip_safe": False,
        })

        # Grab the build_ext command and copy all files back to source dir.
        # Done so Poetry grabs the files during the next step in its build.
        build_ext_cmd = distribution.get_command_obj("build_ext")
        build_ext_cmd.ensure_finalized()
        # Set the value to 1 for "inplace", with the goal to build extensions
        # in build directory, and then copy all files back to the source dir
        # (under the hood, "copy_extensions_to_source" will be called after
        # building the extensions). This is done so Poetry grabs the files
        # during the next step in its build.
        build_ext_cmd.parallel = strtobool(os.environ.get("BUILD_EXTENSIONS_PARALLEL", "True"))
        if build_ext_cmd.parallel:
            logger.info("Building cython extensions in parallel...")
        else:
            logger.info("Building cython extensions serially...")
        build_ext_cmd.inplace = True
        build_ext_cmd.run()
    logger.info("Build successful!")


if __name__ == "__main__":
    build()
