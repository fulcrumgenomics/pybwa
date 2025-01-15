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

