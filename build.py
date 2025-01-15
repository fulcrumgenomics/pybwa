import multiprocessing
from pathlib import Path
from typing import List

from Cython.Build import cythonize
from Cython.Distutils.build_ext import new_build_ext as cython_build_ext
from setuptools import Extension, Distribution

SOURCE_DIR = Path("pybwa")
BUILD_DIR = Path("cython_build")
compile_args = []
link_args = []
include_dirs = ["bwa"]
libraries = ['m', 'z', 'pthread']
library_dirs=['bwa']
extra_objects = []
h_files = []
c_files = []
for root_dir in ["bwa", "pybwa"]:
    h_files.extend(str(x) for x in Path(root_dir).rglob("*.h"))
    c_files.extend(str(x) for x in Path(root_dir).rglob("*.c") if x.name not in ['example.c', 'main.c'])

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

        # Tell Cython we're using Python 3. Becomes default in Cython 3
        compiler_directives={"language_level": "3", 'embedsignature': True},

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
        'packages': ['pybwa'],
        'package_dirs': {'pybwa': 'pybwa'},
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
