#!/usr/bin/env nextflow

println params

reads_ch = Channel
    .from(
        params.kallisto_fastx_single,
        params.kallisto_fastx_forward,
        params.kallisto_fastx_reverse
    )
    .collect()
    .map { it[0] == null ? [true, [file(it[1]), file(it[2])]] : [false, file(it[0])] }


process fullCortexGraph {
    publishDir 'cortex_graph'

    cpus params.jobs

    output:
    file 'full.ctx' into full_cortex_graphs

    """
    #!/usr/bin/env python3
    import os
    from subprocess import run

    cmd = [
        '$params.mccortex', 'build', '$params.mccortex_args',
        '--kmer', '$params.kmer_size',
        '--sample', 'abeona',
        '--threads', '$params.jobs',
    ]
    if '$params.fastx_forward' != 'null':
        cmd += ['--seq2', '$params.fastx_forward:$params.fastx_reverse']
    if '$params.fastx_single' != 'null':
        cmd += ['--seq', '$params.fastx_single']
    cmd.append('full.ctx')
    cmd = [str(c) for c in cmd]
    run(' '.join(cmd), check=True, shell=True)
    """
}


process cleanCortexGraph {
    publishDir 'cortex_graph'
    cpus params.jobs

    input:
    file 'full.ctx' from full_cortex_graphs

    output:
    file 'full.clean.ctx' into clean_cortex_graphs

    """
    MIN_TIP_LENGTH=0
    if [ '$params.prune_tips_with_mccortex' == 'true' ]; then
        MIN_TIP_LENGTH='$params.min_tip_length'
    fi
    $params.mccortex clean $params.mccortex_args \
        -T\$MIN_TIP_LENGTH -U$params.min_unitig_coverage \
        --threads $params.jobs \
        --out full.clean.ctx full.ctx
    """
}

prune_with_cortexpy_ch = Channel.create()
pruned_by_mccortex_ch = Channel.create()
clean_cortex_graphs
    .choice(
        prune_with_cortexpy_ch,
        pruned_by_mccortex_ch
    ) { params.prune_tips_with_mccortex ? 1 : 0 }

process pruneCortexGraphOfTips {
    publishDir 'cortex_graph'

    input:
    file graph from prune_with_cortexpy_ch

    output:
    file "full.clean.min_tip_length_${params.min_tip_length}.ctx" into pruned_by_cortexpy_ch

    """
    #!/usr/bin/env python3
    from subprocess import run

    output = "full.clean.min_tip_length_${params.min_tip_length}.ctx"
    if 0 == $params.min_tip_length:
        run(f'cp $graph {output}', check=True, shell=True)
    else:
        run(f'cortexpy prune --out {output} $graph --remove-tips $params.min_tip_length',
            check=True, shell=True)
    """
}

process traverseCortexSubgraphs {
    publishDir 'traversals'

    input:
    file graph from pruned_by_mccortex_ch.mix(pruned_by_cortexpy_ch)

    output:
    file '*.traverse.ctx' into traversals

    """
    #!/usr/bin/env python3
    from subprocess import run

    cmd = 'abeona subgraphs $graph . -m $params.memory'
    if '$params.initial_contigs' != 'null':
        cmd += ' --initial-contigs $params.initial_contigs'
    run(cmd, shell=True, check=True)
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
    .into{gid_traversals;
	  gid_traversals_for_subgraph_list1_ch;
	  gid_traversals_for_subgraph_list2_ch;
	  gid_traversals_for_assign_reads_ch }


process candidateTranscripts {
    publishDir 'candidate_transcripts'

    input:
    set gid, file(graph) from gid_traversals

    output:
    set(gid, file('g*.candidate_transcripts.fa.gz')) optional true into separate_candidate_transcripts_ch
    set(gid, file('g*.transcripts.fa.gz')) optional true into nonkallisto_single_transcripts_ch
    set(gid, file('g*.candidate_transcripts.fa.gz.skipped'), file(graph)) optional true into skipped_subgraphs_ch

    """
    #!/usr/bin/env python3
    import shutil
    import gzip
    import subprocess
    from yaml import load
    from Bio.SeqIO import parse
    from cortexpy.command import get_exit_code_yaml_path

    CORTEXPY_EXIT_CODES = load(open(get_exit_code_yaml_path(), 'rt'))

    fasta = 'g${gid}.candidate_transcripts.fa.gz'
    transcript = 'g${gid}.transcripts.fa.gz'
    skipped = f'{fasta}.skipped'
    cortexpy_cmd = f'cortexpy traverse --max-paths $params.max_paths_per_subgraph --graph-index ${gid} $graph'
    if '$params.extra_start_kmer' != 'null':
        cortexpy_cmd += ' --extra-start-kmer $params.extra_start_kmer'
    cmd = f'''
    set -o pipefail
    {cortexpy_cmd} | gzip -c > {fasta}.tmp
    '''
    exitcode = subprocess.run(cmd, shell=True, executable='/bin/bash').returncode
    if exitcode == 0:
        shutil.move(f'{fasta}.tmp', fasta)
        n_seqs = 0
        n_long_seqs = 0
        for seq in parse(gzip.open(fasta, "rt"), 'fasta'):
            n_seqs += 1
            if len(seq) > $params.max_read_length:
                n_long_seqs += 1
            if n_seqs > 1 and n_long_seqs > 0:
                break
        if n_seqs == 1:
            shutil.move(fasta, transcript)
        elif n_seqs == 0 or n_long_seqs == 0:
            shutil.move(fasta, skipped)
    elif exitcode != CORTEXPY_EXIT_CODES['MAX_PATH_EXCEEDED']:
        exit(exitcode)
    else:
        shutil.move(f'{fasta}.tmp', skipped)
    """

}

skipped_subgraphs_ch
    .collectFile(storeDir: 'skipped_subgraphs'){ item ->
        [ 'skipped_subgraphs.txt', item.join('\t')  + '\n' ]
    }


process buildKallistoIndices {
    publishDir 'kallisto_indices'
    cpus params.kallisto_threads

    input:
	set gid, file(fasta) from separate_candidate_transcripts_ch

    output:
    set gid, file(fasta), file('*.ki') into kallisto_indices

    """
    #!/usr/bin/env python3
    import os
    import subprocess
    import re

    index = '${fasta}.ki'
    cmd = f'kallisto index -i {index} $fasta'
    if int($params.kmer_size) < 31:
        cmd += ' --kmer-size $params.kmer_size'
    subprocess.check_call(cmd.split(' '))

    cprocess = subprocess.run(['kallisto', 'inspect', index], stdout=subprocess.PIPE, encoding='utf8')
    print(cprocess.stdout)
    if re.search(r'Number of k-mers in index =\\s+0', cprocess.stdout):
        os.remove(index)
    """
}

process createSubgraphList {
    publishDir 'subgraph_list'

    input:
    val(gids) from gid_traversals_for_subgraph_list1_ch
	.map{it[0]}
	.toSortedList()
    val(graphs) from gid_traversals_for_subgraph_list2_ch
	.toSortedList({ a, b -> a[0] <=> b[0] })
        .flatten()
	.collate(2)
	.map{it[1]}
        .collect()

    output:
	file('subgraph_list.txt') into subgraph_list_for_assignment_ch

    """
    #!/usr/bin/env python3
    out = open('subgraph_list.txt', 'w')
    out.write('#prefix\\tgraph\\n')
    for gid, graph in zip([f'g{g}' for g in $gids], '${graphs.join(' ')}'.split()):
        out.write(f'{gid}\\t{graph}\\n')
    """
}

process assignReadsToSubgraphs {
    publishDir 'reads_assigned_to_subgraphs'
    cpus params.jobs

    input:
    set is_paired, file('reads*') from reads_ch
    file('subgraph_list.txt') from subgraph_list_for_assignment_ch

    output:
	file 'g*.fastq.gz' into assigned_reads_ch

    """
    if [ '$is_paired' == 'true' ]; then
        abeona reads subgraph_list.txt reads1 --reverse reads2 --format fastq
    else
        abeona reads subgraph_list.txt reads --format fastq
    fi
    """
}

// Extract graph id and collect read pairs by gid
gid_assigned_reads_ch = assigned_reads_ch
    .flatten()
    .toSortedList()
    .flatten()
    .map { file ->
           def file_name = file.name.toString()
           def match = file_name =~ /g(\d+).\d.fastq.gz$/
           return tuple(match[0][1], file)
    }
    .groupTuple()



process kallistoQuant {
    publishDir 'kallisto_quant'

    cpus params.kallisto_threads

    input:
    set gid, file('fasta'), file(index), file('reads') from kallisto_indices
	.join(gid_assigned_reads_ch)

    output:
    set gid, file(fasta), file('g*') into kallisto_quants

    """
    #!/usr/bin/env python3
    from subprocess import run

    cmd = 'kallisto quant --threads $params.kallisto_threads -i $index --output-dir g$gid -b $params.bootstrap_samples --plaintext'
    if '$params.kallisto_fastx_single' == 'null':
        cmd += ' $reads'
    else:
        cmd += (
            ' -l $params.kallisto_fragment_length -s $params.kallisto_sd'
            ' --single $reads'
        )
    print(cmd)
    run(cmd, check=True, shell=True)
    """
}

process filter_transcripts {
    input:
    set gid, file(candidates), file(kallisto_quant_dir) from kallisto_quants

    output:
    file '*.transcripts.fa.gz' into filtered_transcripts

    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    from pathlib import Path
    import logging
    import gzip
    from Bio import SeqIO
    if '$params.quiet' == 'false':
        log_level = logging.INFO
    else:
        log_level = logging.WARNING
    logging.basicConfig(level=log_level)
    logger = logging.getLogger('abeona.assembly.filter_transcripts')

    def filter_and_annotate_records(records, filtered_counts, abundance, est_count_threshold):
        for record in records:
            if record.id in filtered_counts.index:
                record.description = f'prop_bs_est_counts_ge_{est_count_threshold}=' \
                                     f'{filtered_counts.at[record.id]};est_count={abundance.at[record.id, "est_counts"]}'
                yield record

    output = 'g${gid}.transcripts.fa.gz'
    Path(output).parent.mkdir(exist_ok=True)

    abundance = '$kallisto_quant_dir/abundance.tsv'
    abundance = pd.read_csv(abundance,  sep='\\t', dtype={'target_id': str, 'length': int})
    abundance.set_index('target_id', inplace=True)

    input_abundance = [f'$kallisto_quant_dir/bs_abundance_{i}.tsv' for i in
                       range(int($params.bootstrap_samples))]
    bootstraps = []
    for bs_abundance in input_abundance:
        bootstraps.append(
            pd.read_csv(bs_abundance, sep='\\t', dtype={'target_id': str, 'length': int}))
    bootstraps = pd.concat(bootstraps)
    est_count_threshold = $params.estimated_count_threshold
    keep_prop = $params.bootstrap_proportion_threshold
    ge1_counts = bootstraps.groupby('target_id')['est_counts'].aggregate(
        lambda x: np.sum(x >= est_count_threshold) / len(x))

    keep_counts = ge1_counts[ge1_counts >= keep_prop]

    logger.info(
        f'Keeping contigs with >= {keep_prop} of bootstrapped est_counts >= {est_count_threshold}')
    logger.info(f'keeping {len(keep_counts)} out of {len(ge1_counts)} contigs.')

    with gzip.open(str('$candidates'), 'rt') as fh:
        filtered_records = filter_and_annotate_records(SeqIO.parse(fh, "fasta"),
                                                       keep_counts,
                                                       abundance,
                                                       est_count_threshold)

        logger.info(f'Writing filtered records to {output}')
        with gzip.open(output, 'wt') as fh:
            SeqIO.write(filtered_records, fh, "fasta")
    """

}

all_transcripts = filtered_transcripts
    .mix(nonkallisto_single_transcripts_ch.map{ g,f -> f })
    .collectFile(storeDir: 'transcripts')

process concatTranscripts {
    publishDir 'all_transcripts'

    input:
    file inputs from all_transcripts.collect()

    output:
    file 'transcripts.fa.gz'

    """
    find . -maxdepth 1 -name '*.transcripts.fa.gz' -exec cat {} + > transcripts.fa.gz
    """
}
