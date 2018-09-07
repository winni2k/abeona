abeona
======

abeona v0.23.0

A simple transcriptome assembler based on kallisto and Cortex graphs.

Installation
------------

The easiest way to install abeona is into a `conda <https://conda.io/miniconda.html>`_ environment.

After activating the conda environment, run:

.. code-block:: bash

    conda install abeona -c conda-forge -c bioconda -c wkretzsch

Running
-------

The principal command is `abeona assemble`. This command assembles transcripts from cleaned
short-read RNA-seq reads in FASTA or FASTQ format. For more information, see:

.. code-block:: bash

    abeona assemble --help

License
-------

abeona is distributed under the terms of the
`Apache License, Version 2.0 <https://choosealicense.com/licenses/apache-2.0>`_.
