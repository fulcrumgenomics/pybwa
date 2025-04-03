============
Installation
============

1. Clone the :code:`pybwa` repo (note the :code:`--recurse-submodules` flag, as :code:`bwa` is a submodule)

.. code-block:: bash

   git clone --recurse-submodules https://github.com/fulcrumgenomics/pybwa

2. Install the environment manager :code:`mamba`
3. Install the Python build tool :code:`poetry`

4. Create an environment with Python, Cython, and Pysam:

.. code-block:: bash

   mamba env create -f pybwa.yml

5. Activate the environment:

.. code-block:: bash

   mamba activate pybwa

6. Configure poetry to install into pre-existing virtual environments:

.. code-block:: bash

    poetry config virtualenvs.create false

7. Install :code:`pybwa` into the virtual environment:

.. code-block:: bash

    poetry install

8. Check your build:

.. code-block:: bash

    poetry run pytest


Primary Development Commands
============================

To check and resolve linting issues in the codebase, run:

.. code-block:: bash

    poetry run ruff check --fix

To check and resolve formatting issues in the codebase, run:

.. code-block:: bash

    poetry run ruff format

To check the unit tests in the codebase, run:

.. code-block:: bash

    poetry run pytest

To check the typing in the codebase, run:

.. code-block:: bash

    poetry run mypy

To generate a code coverage report after testing locally, run:

.. code-block:: bash

    poetry run coverage html

To check the lock file is up to date:

.. code-block:: bash

    poetry check --lock

To build the documentation:

.. code-block:: bash

    cd docs
    poetry run make html


Shortcut Task Commands
======================

To be able to run shortcut task commands, first install the Poetry plugin `poethepoet <https://poethepoet.natn.io/index.html>`_:

.. code-block:: bash

    poetry self add 'poethepoet[poetry_plugin]'


For Running Individual Checks
=============================

.. code-block:: bash

    poetry task check-lock
    poetry task check-format
    poetry task check-lint
    poetry task check-tests
    poetry task check-typing
    poetry task build-docs

For Running All Checks
======================

.. code-block:: bash

    poetry task check-all

For Running Individual Fixes
============================

.. code-block:: bash

    poetry task fix-format
    poetry task fix-lint

For Running All Fixes
=====================

.. code-block:: bash

    poetry task fix-all

For Running All Fixes and Checks
================================

.. code-block:: bash

    poetry task fix-and-check-all

Creating a Release on PyPi
==========================

1. Clone the repository recursively and ensure you are on the :code:`main` (un-dirty) branch
2. Checkout a new branch to prepare the library for release
3. Bump the version of the library to the desired SemVer with :code:`poetry version #.#.#`
4. Commit the version bump changes with a Git commit message like :code:`chore(release): bump to #.#.#`
5. Push the commit to the upstream remote, open a PR, ensure tests pass, and seek reviews
6. Squash merge the PR
7. Tag the new commit on the main branch of the origin repository with the new SemVer

.. note::
    This project follows `Semantic Versioning <https://semver.org/>`_.
    In brief:
    
    * `MAJOR` version when you make incompatible API changes
    * `MINOR` version when you add functionality in a backwards compatible manner
    * `PATCH` version when you make backwards compatible bug fixes

GitHub Actions will take care of the remainder of the deployment and release process with:

1. Unit tests will be run for safety-sake
2. A source distribution will be built
3. Multi-arch multi-Python binary distributions will be built
4. Assets will be deployed to PyPi with the new SemVer
5. A `Conventional Commit <https://www.conventionalcommits.org/en/v1.0.0/>`_-aware changelog will be drafted
6. A GitHub release will be created with the new SemVer and the drafted changelog

.. warning::
    Consider editing the changelog if there are any errors or necessary enhancements.
