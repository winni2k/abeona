abeona
======

.. image:: https://travis-ci.org/winni2k/abeona.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/winni2k/abeona

.. image:: https://img.shields.io/github/commits-since/winni2k/abeona/v0.45.0.svg
    :alt: Commits since latest release
    :target: https://github.com/winni2k/abeona/compare/v0.45.0...master

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat
    :alt: freedom
    :target: http://bioconda.github.io/recipes/abeona/README.html

abeona v0.45.0

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

To create a conda environment named "abeona" run:

.. code-block:: bash

    conda install create -n abeona abeona -c conda-forge -c bioconda -y
    
    # activate the environment
    conda activate abeona

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
    $ abeona assemble -k 5 -m 4 --fastx-single input.fa --min-tip-length 1 \
        --kallisto-fragment-length 7 --kallisto-sd 1 -o test --no-links
    N E X T F L O W  ~  version 19.01.0
    Launching `assemble.nf` [marvelous_rosalind] - revision: 78bfc98f1f

    Running nextflow script for abeona v0.45.0
    Nextflow arguments from args.json:
    {
        "assemble_unassembled_reads_with_transabyss": false,
        "bootstrap_proportion_threshold": 0.95,
        "bootstrap_samples": 100,
        "estimated_count_threshold": 1,
        "extra_start_kmer": null,
        "fastx_forward": null,
        "fastx_reverse": null,
        "fastx_single": "/Volumes/mac-3/home/Projects/caina/2020-08-24_abeona/input.fa",
        "initial_contigs": null,
        "jobs": 2,
        "kallisto_fragment_length": 7.0,
        "kallisto_sd": 1.0,
        "kallisto_threads": 2,
        "kmer_size": "5",
        "max_junctions": 0,
        "max_paths_per_subgraph": 0,
        "max_read_length": 7,
        "memory": 4,
        "min_tip_length": 1,
        "min_unitig_coverage": 4,
        "no_cleanup": false,
        "no_links": true,
        "no_prune_tips_with_mccortex": false,
        "out_dir": "test",
        "prune_tips_iteratively": false,
        "prune_tips_with_mccortex": true,
        "quiet": false,
        "record_buffer_size": -1,
        "report_unassembled_reads": false,
        "resume": false,
        "with_dag": false,
        "with_report": false,
        "mccortex": "mccortex 5",
        "mccortex_args": "--sort --force -m 4G",
        "mccortex_thread_args": "--force -m 2G"
    }

    [warm up] executor > local
    [3d/255de7] Submitted process > fullCortexGraph
    [2c/c20fdf] Submitted process > cleanCortexGraph
    [b0/909834] Submitted process > traverseCortexSubgraphs (1)
    [b2/c75eec] Submitted process > createSubgraphList
    [18/2c45c7] Submitted process > assignReadsToSubgraphs
    [f8/b4c61b] Submitted process > threadReads (1)
    [d9/7e6528] Submitted process > candidateTranscripts (1)
    [25/ed539b] Submitted process > buildKallistoIndices (1)
    [c7/9c5b50] Submitted process > kallistoQuant (1)
    [e7/5b2330] Submitted process > filter_transcripts (1)
    [9d/3f9fa5] Submitted process > concatTranscripts

    # View the resulting assembled transcripts
    $ cat test/transcripts.fa
    >g0_p0 prop_bs_est_counts_ge_1=0.99;est_count=3
    AAAAAGGG
    >g0_p1 prop_bs_est_counts_ge_1=0.97;est_count=3
    AAAAACCC

Development
-----------

::

    conda env create -f environment.yml my-dev-env
    conda activate my-dev-env
    make test

License
-------

Abeona is distributed under the terms of the
`Apache License, Version 2.0 <https://choosealicense.com/licenses/apache-2.0>`_.

Citing
------

If you use abeona in your research, please cite:

    Akhter S, Kretzschmar WW, Nordal V, Delhomme N, Street NR, Nilsson O, Emanuelsson O, Sundström JF. Integrative Analysis of Three RNA Sequencing Methods Identifies Mutually Exclusive Exons of MADS-Box Isoforms During Early Bud Development in Picea abies. Front. Plant Sci. 9, 1–18 (2018).

Changelog
---------

Version 0.45.0
~~~~~~~~~~~~~~

:Date: XXX

New features
............

**abeona assemble**

* Mccortex is now used for pruning by default
* The command line argument ``--prune-tips-with-mccortex`` is now deprecated.
  Instead use ``--no-prune-tips-with-mccortex``.
* New iterative pruning strategy ``--prune-tips-iteratively``.

Version 0.44.0
~~~~~~~~~~~~~~

:Date: 2019-03-26

This version skips commits made for the 0.43.0 tag.

New features
............

* Reads that share kmers with subgraphs that are skipped are now reported in the
  ``unassembled_reads`` directory.

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
