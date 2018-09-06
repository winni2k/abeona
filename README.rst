abeona
======

abeona v0.22.3

A simple transcriptome assembler based on kallisto and Cortex graphs.

Installation
------------

It is recommended to install abeona and its dependencies into a
`conda <https://conda.io/miniconda.html>`_ environment.

Conda dependencies
~~~~~~~~~~~~~~~~~~

At present, abeona depends on the following conda packages:

- zlib
- nextflow
- kallisto

Manually installable dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Furthermore, abeona needs a working version of `Mccortex <https://github.com/mcveanlab/mccortex>`_
to be installed in the PATH:

.. code-block:: bash

    git clone --recursive https://github.com/mcveanlab/mccortex
    cd mccortex
    make all MAXK=31
    make all MAXK=63
    # copy all files in bin to a directory in your PATH

Abeona
~~~~~~

abeona is distributed on `PyPI <https://pypi.org>`_ as a universal
wheel and is available on Linux/macOS and Windows, and supports
Python 3.6+.

.. code-block:: bash

    $ pip install abeona

License
-------

abeona is distributed under the terms of the
`Apache License, Version 2.0 <https://choosealicense.com/licenses/apache-2.0>`_.
