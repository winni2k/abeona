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

.. |commits-since| image:: https://img.shields.io/github/commits-since/winni2k/abeona/v0.43.0.svg
    :alt: Commits since latest release
    :target: https://github.com/winni2k/abeona/compare/v0.43.0...master


abeona v0.43.0

A simple transcriptome assembler based on kallisto and Cortex graphs.

Abeona consists of the following stages:

1. Assembly of reads into a De Bruijn graph
2. Pruning of tips and low-coverage unitigs
3. Partitioning of the De Bruijn graph into subgraphs
4. Generation of candidate transcripts by simple path traversal
5. Filtering of candidate transcripts by kallisto

Installation
------------

The easiest way to install abeona is into a `conda <https://conda.io/miniconda.html>`_ environment.

After activating the conda environment, run:

.. code-block:: bash

    conda install abeona -c conda-forge -c bioconda

Usage
-----

The principal command is ``abeona assemble``. This command assembles transcripts from cleaned
short-read RNA-seq reads in FASTA or FASTQ format. A description of command arguments is
available with the command:

.. code-block:: bash

    abeona assemble --help

Specifying input read data
~~~~~~~~~~~~~~~~~~~~~~~~~~

Abeona is designed to be run on reads from one biological sample at a time.
Abeona uses sequencing reads in two stages: for De Bruijn-graph construction,
and for candidate transcript filtering with kallisto. The first stage accepts
paired-end, single-end, or both types of reads through the ``--fastx-*`` arguments.
The reads for the second stage are specified with the ``--kallisto-fastx-*`` arguments.
Kallisto only accepts single-end or paired-end reads, so input to this stage
is also restricted in that manner.

Toy Example
-----------

.. code-block:: bash

    # Let's create a FASTA consisting of sub-reads from two transcripts: AAAAACCC and AAAAAGGG
    $ for s in AAAAACC AAAAAGG AAAACCC AAAAGGG; do for i in $(seq 1 3); do echo -e ">_\n$s" >> input.fa; done; done

    # Now feed the fasta to the graph assembly step with --fastx-single and to the kallisto filtering
    # step with --kallisto-fastx-single.
    $ abeona assemble -k 5 -m 4 --fastx-single input.fa --kallisto-fastx-single \
        input.fa --kallisto-fragment-length 7 --kallisto-sd 1 -o test --no-links
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

Abeona is distributed under the terms of the
`Apache License, Version 2.0 <https://choosealicense.com/licenses/apache-2.0>`_.

Citing
------

If you use abeona in your research, please cite:

    Akhter S, Kretzschmar WW, Nordal V, Delhomme N, Street NR, Nilsson O, Emanuelsson O, Sundström JF. Integrative Analysis of Three RNA Sequencing Methods Identifies Mutually Exclusive Exons of MADS-Box Isoforms During Early Bud Development in Picea abies. Front. Plant Sci. 9, 1–18 (2018).``

Changelog
---------

Version 0.42.0
~~~~~~~~~~~~~~

:Date: 2018-12-17

Interface Changes
.................

* Cleanup now deletes all directories in output dir except for ``all_transcripts/transcripts.fa.gz``
* Cleanup is now on by default
* Cleanup can be turned off with ``--no-cleanup`` flag
* ``all_transcripts/transcripts.fa.gz`` is unzipped and stored as ``transcripts.fa`` to conform
to the convention set by Trinity and Oases for output file names

Version 0.41.0
~~~~~~~~~~~~~~

:Date: 2018-12-13

Interface changes
.................

* Remove ``--kallisto-fastx-*`` arguments. Being able to separately specify reads to graph building
and kallisto has not been all that useful, and it increases the complexity of the code.
* Add default value of ``--kmer-size`` for ``--min-tip-length``.

Fixes
.....

* There are several ways in which kallisto can fail due to no reads pseudoaligning to a subgraph's
candidate transcripts. When this happens, abeona now catches the error and silently ignores the
subgraph.


Version 0.40.0
~~~~~~~~~~~~~~

:Date: 2018-11-17

New features
............

* Add ``--no-links`` argument to turn off link use in candidate transcript creation
* Add ``--max-junctions`` argument to allow fast skipping of subgraphs with too many junctions

Fixes
.....

* Properly assign reads to all subgraphs to which they are assignable
* Solve high-mem use problem by creating links only on assigned reads

Version 0.36.0
~~~~~~~~~~~~~~

:Date: 2018-10-25

New features
............

* Graph traversal now uses links

Fixes
.....

* Lots of improvements to ``abeona reads`` to improve memory and filehandle use

Version 0.33.0
~~~~~~~~~~~~~~

:Date: 2018-10-17

New features
............

* Use kmer mapping (``abeona reads``) to assign reads to subgraphs before quantification of
candidate transcripts with kallisto

Fixes
.....

* Add missing conda dependency ``seqtk`` to ``environment.yml`` for travis CI
