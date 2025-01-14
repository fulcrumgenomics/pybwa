============
Installation
============

1. Install the environment manager :code:`mamba`
2. Install the Python build tool :code:`poetry`

3. Create an environment with Python, Cython, and Pysam:

.. code-block:: bash

   mamba env create -f pybwa.yml

4. Activate the environment:

.. code-block:: bash

   mamba activate pybwa

5. Clone the :code:`bwa` repo

.. code-block:: bash

   git clone https://github.com/lh3/bwa

6. Configure poetry to install into pre-existing virtual environments:

.. code-block:: bash

    poetry config virtualenvs.create false

7. Install :code:`pybwa` into the virtual environment:

.. code-block:: bash

    poetry install

8. Check your build:

.. code-block:: bash

    poetry run pytest

