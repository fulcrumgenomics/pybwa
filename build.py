import multiprocessing
import os
import platform
import shutil
import subprocess
from contextlib import contextmanager
from pathlib import Path
from typing import List

from Cython.Build import cythonize
from Cython.Distutils.build_ext import new_build_ext as cython_build_ext
from setuptools import Extension, Distribution


@contextmanager
def changedir(path):
    save_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(save_dir)

@contextmanager
def with_patches():
    patches = sorted([
        os.path.abspath(patch)
        for patch in Path("patches").iterdir()
        if patch.is_file() and patch.suffix == ".patch"
    ])
    has_git = shutil.which("git") is not None
    has_git_dir = Path(".git").exists() and Path(".git").is_dir()
    use_git = has_git and has_git_dir
    with changedir("bwa"):
        for patch in patches:
            if use_git:
                retcode = subprocess.call(f"git apply {patch}", shell=True)
            else:
                retcode = subprocess.call(f"patch -p1 < {patch}", shell=True)
            if retcode != 0:
                raise RuntimeError(f"Failed to apply patch {patch}")
    try:
        yield
    finally:
        if use_git:
            commands = ["git submodule deinit -f .", "git submodule update --init"]
            for command in commands:
                retcode = subprocess.call(command, shell=True)
                if retcode != 0:
                    raise RuntimeError(f"Failed to reset submodules: {command}")

compiler_directives = {
    "language_level": "3",
    'embedsignature': True,
}
SOURCE_DIR = Path("pybwa")
BUILD_DIR = Path("cython_build")
compile_args = []
link_args = []
include_dirs = ["bwa", "pybwa"]
libraries = ['m', 'z', 'pthread']
if platform.system() == 'Linux':
    libraries.append("rt")
library_dirs=['pybwa', 'bwa']
extra_objects = []
define_macros = [("HAVE_PTHREAD", None), ("USE_MALLOC_WRAPPERS", None)]
h_files = []
c_files = []

exclude_files = {
    "pybwa": ["libbwaaln.c", "libbwaindex.c", "libbwamem.c"],
    "bwa": ['example.c', 'main.c']
}
for root_dir in library_dirs:
    h_files.extend(str(x) for x in Path(root_dir).rglob("*.h"))
    c_files.extend(str(x) for x in Path(root_dir).rglob("*.c") if x.name not in exclude_files[root_dir])

# Check if we should build with linetracing for coverage
build_with_coverage = os.environ.get("CYTHON_TRACE", "false").lower() in ("1", "true", '"true"')
if build_with_coverage:
    compiler_directives["linetrace"] = True
    compiler_directives['emit_code_comments'] = True
    define_macros.extend([('CYTHON_TRACE', 1), ('CYTHON_TRACE_NOGIL', 1), ('DCYTHON_USE_SYS_MONITORING', 0)])
    BUILD_DIR = Path(".")  # the compiled .c files need to be next to the .pyx files for coverage

if platform.system() != 'Windows':
    compile_args.extend([
        "-Wno-unused-result",
        "-Wno-unreachable-code",
        "-Wno-single-bit-bitfield-constant-conversion",
        "-Wno-deprecated-declarations",
        "-Wno-unused",
        "-Wno-strict-prototypes",
        "-Wno-sign-compare",
        "-Wno-error=declaration-after-statement"
    ])

libbwaindex_module = Extension(
    name='pybwa.libbwaindex',
    sources=['pybwa/libbwaindex.pyx'] + c_files,
    depends=h_files,
    extra_compile_args=compile_args,
    extra_link_args=link_args,
    extra_objects=extra_objects,
    include_dirs=include_dirs,
    language='c',
    libraries=libraries,
    library_dirs=library_dirs,
    define_macros=define_macros
)

libbwaaln_module = Extension(
    name='pybwa.libbwaaln',
    sources=['pybwa/libbwaaln.pyx'] + c_files,
    depends=h_files,
    extra_compile_args=compile_args,
    extra_link_args=link_args,
    extra_objects=extra_objects,
    include_dirs=include_dirs,
    language='c',
    libraries=libraries,
    library_dirs=library_dirs,
    define_macros=define_macros
)

libbwamem_module = Extension(
    name='pybwa.libbwamem',
    sources=['pybwa/libbwamem.pyx'] + c_files,
    depends=h_files,
    extra_compile_args=compile_args,
    extra_link_args=link_args,
    extra_objects=extra_objects,
    include_dirs=include_dirs,
    language='c',
    libraries=libraries,
    library_dirs=library_dirs,
    define_macros=define_macros
)


def cythonize_helper(extension_modules: List[Extension]) -> List[Extension]:
    """Cythonize all Python extensions"""

    return cythonize(
        module_list=extension_modules,

        # Don't build in source tree (this leaves behind .c files)
        build_dir=BUILD_DIR,

        # Don't generate an .html output file. Would contain source.
        annotate=False,

        # Parallelize our build
        nthreads=multiprocessing.cpu_count() * 2,

        # Compiler directives (e.g. language, or line tracing for coverage)
        compiler_directives=compiler_directives,

        # (Optional) Always rebuild, even if files untouched
        force=True,
    )

CLASSIFIERS = '''
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
'''


def build():
    # apply patches to bwa, then revert them after
    with with_patches():
        # Collect and cythonize all files
        extension_modules = cythonize_helper([
            libbwaindex_module,
            libbwaaln_module,
            libbwamem_module
        ])

        # Use Setuptools to collect files
        distribution = Distribution({
            "name": "pybwa",
            'version': '0.0.1',
            'description': 'Python bindings for BWA',
            'long_description': __doc__,
            'long_description_content_type': 'text/x-rst',
            'author': 'Nils Homer',
            'author_email': 'nils@fulcrumgenomics.com',
            'license': 'MIT',
            'platforms': ['POSIX', 'UNIX', 'MacOS'],
            'classifiers': [_f for _f in CLASSIFIERS.split('\n') if _f],
            'url': 'https://github.com/fulcrumgenomics/pybwa',
            'packages': ['pybwa', 'pybwa.include.bwa', 'pybwa.include.patches'],
            'package_dir': {'pybwa': 'pybwa', 'pybwa.include.bwa': 'bwa', 'pybwa.include.patches': 'patches' },
            'package_data': {'': ['*.pxd', '*.h', '*.c', 'py.typed', '*.pyi', '*.patch'], },
            "ext_modules": extension_modules,
            "cmdclass": {
                "build_ext": cython_build_ext,
            },
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
        build_ext_cmd.inplace = 1
        build_ext_cmd.run()


if __name__ == "__main__":
    build()
