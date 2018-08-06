#!/usr/bin/env nextflow

println params

process fullCortexGraph {
    publishDir 'cortex_graph'

    output:
    file 'full.ctx' into full_cortex_graphs

    """
    #!/usr/bin/env python3
    import os

    cmd = [
        '$params.mccortex', 'build', '$params.mccortex_args',
        '--kmer', '$params.kmer_size',
        '--sample', 'abeona',
        # '--threads', threads,
    ]
    if '$params.fastx_forward' != 'null':
        cmd += ['--seq2', '$params.fastx_forward:$params.fastx_reverse']
    if '$params.fastx_single' != 'null':
        cmd += ['--seq', '$params.fastx_single']
    cmd.append('full.ctx')
    cmd = [str(c) for c in cmd]
    os.system(' '.join(cmd))
    """
}


process cleanCortexGraph {
    publishDir 'cortex_graph'

    input:
    file 'full.ctx' from full_cortex_graphs

    output:
    file 'full.clean.ctx' into clean_cortex_graphs

    """
    $params.mccortex clean $params.mccortex_args \
        -T0 -U$params.min_unitig_coverage \
        --out full.clean.ctx full.ctx
    """
}

process pruneCortexGraphOfTips {
    publishDir 'cortex_graph'

    input:
    file graph from clean_cortex_graphs

    output:
    file "full.clean.min_tip_length_${params.min_tip_length}.ctx" into pruned_cortex_graphs

    """
    #!/usr/bin/env python3
    import os

    output = "full.clean.min_tip_length_${params.min_tip_length}.ctx"
    if 0 == $params.min_tip_length:
        os.system(f"cp $graph {output}")
    else:
        os.system(f'cortexpy prune --out {output} $graph --remove-tips $params.min_tip_length')
    """
}

process traverseCortexSubgraph {
    publishDir 'traversals'

    input:
    file graph from pruned_cortex_graphs

    output:
    file '*.traverse.ctx' into traversals

    """
    #!/usr/bin/env python3
    import os

    cmd = 'python -m abeona subgraphs $graph . -m $params.memory'
    if '$params.initial_contigs' != 'null':
        cmd += ' --initial-contigs $params.initial_contigs'
    os.system(cmd)
    """
}

// Extract graph id and pass it along in the pipeline
traversals
    .flatten()
    .map { file ->
           def file_name = file.name.toString()
           def match = file_name =~ /g(\d+)/
           return tuple(match[0][1], file)
     }
    .groupTuple()
    .set{ gid_traversals }


process candidateTranscripts {
    publishDir 'candidate_transcripts'

    input:
    set gid, file(graph) from gid_traversals

    output:
    set(gid, file('*.fa.gz')) optional true into candidate_transcripts

    """
    #!/usr/bin/env python3
    import os
    import shutil
    from yaml import load
    from cortexpy.command import get_exit_code_yaml_path
    CORTEXPY_EXIT_CODES = load(open(get_exit_code_yaml_path(), 'rt'))

    fasta = 'g${gid}.candidate_transcripts.fa.gz'
    cmd = f'''
    set -o pipefail
    cortexpy view traversal $graph --max-paths $params.max_paths_per_subgraph \
           | gzip -c > {fasta}.tmp
    '''
    exitcode = os.system(cmd) >> 8
    print(exitcode)
    if exitcode == 0:
        shutil.move(f'{fasta}.tmp', fasta)
    elif exitcode != CORTEXPY_EXIT_CODES['MAX_PATH_EXCEEDED']:
        exit(exitcode)
    """

}

process buildKallistoIndices {
    publishDir 'kallisto_indices'

    input:
    set gid, file(fasta) from candidate_transcripts

    output:
    set gid, file(fasta), file('*.ki') into kallisto_indices

    """
    #!/usr/bin/env python3
    import os

    cmd = 'kallisto index -i ${fasta}.ki $fasta'
    if int($params.kmer_size) < 31:
        cmd += ' --kmer-size $params.kmer_size'
    os.system(cmd)
    """
}

process kallistoQuant {
    publishDir 'kallisto_quant'

    input:
    set gid, file(fasta), file(index) from kallisto_indices

    output:
    set gid, file(fasta), file('g*') into kallisto_quants

    """
    #!/usr/bin/env python3
    import os

    cmd = 'kallisto quant -i $index --output-dir g$gid -b $params.bootstrap_samples --plaintext'
    if '$params.kallisto_fastx_forward' != 'null':
        cmd += ' $params.kallisto_fastx_forward $params.kallisto_fastx_reverse'
    else:
        cmd += (
            ' -l $params.kallisto_fragment_length -s $params.kallisto_sd'
            ' --single $params.kallisto_fastx_single'
        )
    print(cmd)
    os.system(cmd)
    """
}

process filter_transcripts {
    publishDir 'transcripts'

    input:
    set gid, file(candidates), file(kallisto_quant_dir) from kallisto_quants

    output:
    file '*.transcripts.fa.gz' into transcripts

    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import logging
    import gzip
    from Bio import SeqIO

    logger = logging.getLogger('abeona.assembly.filter_transcripts')

    def filter_and_annotate_contigs(filtered_counts, candidate_transcripts):
        with gzip.open(str(candidate_transcripts), 'rt') as fh:
            for record in SeqIO.parse(fh, "fasta"):
                if record.id in filtered_counts.index:
                    record.description = 'prop_bs_est_counts_ge_1={}'.format(filtered_counts.at[record.id])
                    yield record



    output = 'g${gid}.transcripts.fa.gz'
    Path(output).parent.mkdir(exist_ok=True)

    bootstraps = []
    input_abundance = [f'$kallisto_quant_dir/bs_abundance_{i}.tsv' for i in
                       range(int($params.bootstrap_samples))]
    for bs_abundance in input_abundance:
        bootstraps.append(
            pd.read_csv(bs_abundance, sep='\\t', dtype={'target_id': str, 'length': int}))
    bootstraps = pd.concat(bootstraps)
    est_count_threshold = 1
    keep_prop = 0.95
    ge1_counts = bootstraps.groupby('target_id')['est_counts'].aggregate(
        lambda x: np.sum(x >= est_count_threshold) / len(x))

    keep_counts = ge1_counts[ge1_counts >= keep_prop]

    logging.info(
        f'Keeping contigs with >= {keep_prop} of bootstrapped est_counts >= {est_count_threshold}')
    logger.info(f'keeping {len(keep_counts)} out of {len(ge1_counts)} contigs.')

    filtered_records = filter_and_annotate_contigs(keep_counts, '$candidates')

    logger.info(f'Writing filtered records to {output}')
    with gzip.open(output, 'wt') as fh:
        SeqIO.write(filtered_records, fh, "fasta")
    """

}
