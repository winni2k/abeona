abeona
======

.. start-badges

.. list-table::
    :stub-columns: 1

    * - tests
      - | |travis|
    * - package
      - | |commits-since|

.. |travis| image:: https://travis-ci.org/winni2k/abeona.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/winni2k/abeona

.. |commits-since| image:: https://img.shields.io/github/commits-since/winni2k/abeona/v0.34.2.svg
    :alt: Commits since latest release
    :target: https://github.com/winni2k/abeona/compare/v0.34.2...master


abeona v0.34.2

A simple transcriptome assembler based on kallisto and Cortex graphs.

Installation
------------

The easiest way to install abeona is into a `conda <https://conda.io/miniconda.html>`_ environment.

After activating the conda environment, run:

.. code-block:: bash

    conda install abeona -c conda-forge -c bioconda

Running
-------

The principal command is `abeona assemble`. This command assembles transcripts from cleaned
short-read RNA-seq reads in FASTA or FASTQ format. For more information, see:

.. code-block:: bash

    abeona assemble --help

Toy Example
~~~~~~~~~~~

.. code-block:: bash

    # Let's create a FASTA consisting of sub-reads from two transcripts: AAAAACCC and AAAAAGGG
    $ for s in AAAAACC AAAAAGG AAAACCC AAAAGGG; do for i in $(seq 1 3); do echo -e ">_\n$s" >> input.fa; done; done

    # Now feed the fasta to the graph assembly step with --fastx-single and to the kallisto filtering
    # step with --kallisto-fastx-single.
    $ abeona assemble -k 5 -m 4 --fastx-single input.fa --kallisto-fastx-single input.fa --kallisto-fragment-length 7 --kallisto-sd 1 -o test
    N E X T F L O W  ~  version 0.31.1
    Launching `assemble.nf` [determined_allen] - revision: 11c20ed355
    [bootstrap_samples:100, fastx_forward:null, fastx_reverse:null, fastx_single:/Users/winni/tmp/input.fa, initial_contigs:null, jobs:2, kallisto_fastx_forward:null, kallisto_fastx_reverse:null, kallisto_fastx_single:/Users/winni/tmp/input.fa, kallisto_fragment_length:7.0, kallisto_sd:1.0, kmer_size:5, max_paths_per_subgraph:0, memory:4, merge_candidates_before_kallisto:false, min_tip_length:0, min_unitig_coverage:4, out_dir:test, quiet:false, resume:false, mccortex:mccortex 5, mccortex_args:--sort --force -m 4G]
    [warm up] executor > local
    [26/119d41] Submitted process > fullCortexGraph
    [fc/585605] Submitted process > cleanCortexGraph
    [dd/40b5fc] Submitted process > pruneCortexGraphOfTips
    [36/f63343] Submitted process > traverseCortexSubgraphs
    [23/6d9033] Submitted process > candidateTranscripts (1)
    [d5/05d417] Submitted process > buildKallistoIndices (1)
    [ac/e36d53] Submitted process > kallistoQuant (1)
    [ec/2b258d] Submitted process > filter_transcripts (1)
    [49/d4c7e3] Submitted process > concatTranscripts

    # View the resulting assembled transcripts
    $ zcat test/all_transcripts/transcripts.fa.gz
    >g0_p0 prop_bs_est_counts_ge_1=0.98
    AAAAAGGG
    >g0_p1 prop_bs_est_counts_ge_1=1.0
    AAAAACCC

License
-------

abeona is distributed under the terms of the
`Apache License, Version 2.0 <https://choosealicense.com/licenses/apache-2.0>`_.


Changelog
---------

v0.34.2
~~~~~~~

:Date: 2018-10-17

New features
............

* Use kmer mapping (``abeona reads``) to assign reads to subgraphs before quantification of
candidate transcripts with kallisto

Fixes
.....

* Add missing conda dependency ``seqtk`` to ``environment.yml`` for travis CI
