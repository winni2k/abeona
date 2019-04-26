import logging
import sys

logger = logging.getLogger('abeona')


def main(argv=sys.argv):
    from . import __version__
    import argparse
    subcommands = {
        'assemble': assemble_main,
        'subgraphs': subgraphs_main,
        'reads': reads_main,
    }
    parser = argparse.ArgumentParser(description=f'abeona version {__version__}', prog='abeona')
    parser.add_argument('-v', '--version', action='version',
                        version=f'%(prog)s version {__version__}')
    parser.add_argument('sub_command', choices=sorted(subcommands.keys()),
                        help='abeona sub-command')
    parser.add_argument('args', nargs=argparse.REMAINDER, help='sub-command arguments')
    args = parser.parse_args(argv[1:])

    import warnings
    from Bio import BiopythonWarning
    warnings.simplefilter('ignore', BiopythonWarning)
    return subcommands[args.sub_command](args.args)


def assemble_main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Run abeona assembly pipeline.',
                                     prog='abeona assemble',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--out-dir', help='Output directory', required=True)
    parser.add_argument('-j', '--jobs', help='Number of jobs to schedule concurrently', default=2,
                        type=int)
    parser.add_argument('-k', '--kmer-size', default=47,
                        help='k-mer size to use to construct the De Bruijn graph')
    parser.add_argument('-m', '--memory', default=3, help='Maximum memory to use in giga bytes',
                        type=int)
    parser.add_argument('-q', '--quiet', action='store_true')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--no-cleanup', action='store_true')
    parser.add_argument('--with-report', action='store_true',
                        help='Create nextflow report file at nexttflow_report.html in output directory')
    parser.add_argument('--with-dag', action='store_true',
                        help='Create flowchart of workflow at flowchart.png in output directory')

    parser.add_argument('--fastx-forward', help='Forward sequences in FASTA/FASTQ format',
                        required=False)
    parser.add_argument('--fastx-reverse', help='Reverse sequences in FASTA/FASTQ format',
                        required=False)
    parser.add_argument('--fastx-single', help='Single-end sequences in FASTA/FASTQ format',
                        required=False)
    parser.add_argument('--initial-contigs', help='Only start assembly from contigs in this FASTA',
                        required=False)
    parser.add_argument('--extra-start-kmer',
                        help='Disconnect this k-mer from incoming k-mers before '
                             'candidate transcript creation. '
                             'This may be useful when assembling circular genomes. '
                             'Best used with --initial-contigs to make sure the k-mer exists '
                             'in the consistent cortexpy graph.',
                        required=False)

    group = parser.add_argument_group('Graph traversal cleaning')
    group.add_argument('--min-tip-length', type=int, default=-1,
                       help='Prune tips shorter than this value.'
                            ' A value of -1 sets the min tip length to the value of --kmer-size')
    group.add_argument('--min-unitig-coverage', type=int, default=4,
                       help="Prune unitigs with mean coverage below this value")
    group.add_argument('--no-prune-tips-with-mccortex', action='store_true',
                       help="Instead of Mccortex use cortexpy to prune unitigs. "
                            "This is slower and not recommended at this time.")
    group.add_argument('--prune-tips-with-mccortex', action='store_true')
    group.add_argument('--prune-tips-iteratively', action='store_true',
                       help="Prune the graph of tip lengths x = 2^n while x is less than --min-tip-length. "
                            "Finally, prune graph of tips shorter than --min-tip-length. "
                            "Currently only works when pruning with Mccortex. ")

    group = parser.add_argument_group('Candidate transcript creation')
    group.add_argument('--max-paths-per-subgraph', type=int, default=0,
                       help='Ignore graphs that have more than this number of paths')
    group.add_argument('--no-links', action='store_true',
                       help='Do not use links in candidate transcript creation')
    group.add_argument('--report-unassembled-reads', action='store_true',
                       help='Store reads from ignored graphs in unassembled_reads')
    group.add_argument('--assemble-unassembled-reads-with-transabyss', action='store_true',
                       help='Try and assemble reads from ignored graphs with transabyss')

    group = parser.add_argument_group('kallisto arguments',
                                      description='Arguments passed directly on to kallisto')
    group.add_argument('--bootstrap-samples', type=int, default=100,
                       help='Number of kallisto bootstrap samples')
    group.add_argument('--kallisto-fragment-length', type=float)
    group.add_argument('--kallisto-sd', type=float)
    group.add_argument('--kallisto-threads', type=int, default=2,
                       help='Number of logical cores to assign to a single kallisto quant job. '
                            'Needs to be less than or equal to --jobs')
    group.add_argument('--max-read-length', type=int,
                       help='Length of longest read in data. Remove candidate subgraphs that do '
                            'not have at least one candidate transcript greater than this length. '
                            'Kallisto cannot align reads to transcripts that are shorter than the '
                            'read. If this value is not specified, then abeona estimates this '
                            'value from the head of the reads supplied to kallisto '
                            '(--kallisto-fastx-*).')

    group = parser.add_argument_group('Candidate transcript filtering')
    group.add_argument('--estimated-count-threshold', type=float, default=1,
                       help='Threshold over which the estimated transcript count from kallisto is '
                            'counted towards keeping a transcript')
    group.add_argument('--bootstrap-proportion-threshold', type=float, default=0.95,
                       help="Proportion of bootstrap iterations for which a transcript's estimated "
                            "counts must be above the --estimated-count-threshold")
    group.add_argument('--record-buffer-size', type=int, default=-1,
                       help='Number of reads to buffer in memory when assigning reads to subgraphs')
    group.add_argument('--max-junctions', type=int, default=0,
                       help='The max junctions argument can be used to quickly ignore large subgraphs'
                            ' with too many junctions to process effectively.')

    args = parser.parse_args(args=argv)
    args = preprocess_assemble_args(args)
    validate_assemble_args(args, parser)

    from pathlib import Path
    import subprocess

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    script_name = 'assemble.nf'
    pipeline_script = out_dir / script_name
    package_script = Path(__file__).parent / 'assemble' / script_name
    if pipeline_script.is_symlink():
        pipeline_script.unlink()
    pipeline_script.symlink_to(package_script)

    args_file = write_args(args, out_dir)

    cmd = f'cd {out_dir} && KALLISTO_THREADS={args.kallisto_threads} JOBS={args.jobs} nextflow run {script_name} -process.maxForks {args.jobs} -params-file {args_file}'
    if args.with_report:
        cmd += ' -with-report nexttflow_report.html'
    if args.with_dag:
        cmd += ' -with-dag flowchart.png'
    if args.resume:
        cmd += ' -resume'
    logger.info(cmd)
    cprocess = subprocess.run(cmd, shell=True)
    if cprocess.returncode != 0:
        return cprocess.returncode
    gzipped_transcripts = out_dir / "all_transcripts" / "transcripts.fa.gz"
    unzipped_transcripts = out_dir / "transcripts.fa"
    if gzipped_transcripts.is_file():
        subprocess.run(
            f'gzip -dc {gzipped_transcripts} > {unzipped_transcripts}',
            shell=True
        )
    else:
        unzipped_transcripts.touch(exist_ok=True)

    import shutil
    if not args.no_cleanup:
        for d in out_dir.iterdir():
            if d.is_dir() and not d.name.startswith('.') and not d.name == 'unassembled_reads':
                shutil.rmtree(d)
    return cprocess.returncode


def write_args(args, out_dir):
    import json

    args_file = 'args.json'
    args_dict = {a: getattr(args, a) for a in dir(args) if not a.startswith('_')}
    args_dict['mccortex'] = f'mccortex {args.kmer_size}'
    args_dict['mccortex_args'] = f'--sort --force -m {args.memory}G'
    if args.max_read_length is None:
        max_read_length = estimate_max_read_length(args)
        assert max_read_length is not None
        args_dict['max_read_length'] = max_read_length
    args_dict['mccortex_thread_args'] = f'--force -m {args.memory // args.jobs}G'
    with open(out_dir / args_file, 'w') as fh:
        json.dump(args_dict, fh)
    return args_file


def estimate_max_read_length(args):
    from .utils import get_maybe_gzipped_file_handle, get_fastx_title_seq_generator
    from itertools import islice
    n_reads_to_read = 100
    max_read_length = None
    for arg in ['fastx_single', 'fastx_forward', 'fastx_reverse']:
        file = getattr(args, arg)
        if file is not None:
            with get_maybe_gzipped_file_handle(file, 'rt') as fh:
                for _, seq in islice(get_fastx_title_seq_generator(fh), n_reads_to_read):
                    if max_read_length is None or len(seq) > max_read_length:
                        max_read_length = len(seq)
    return max_read_length


def subgraphs_main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Partition cortex graph into its subgraphs.',
                                     prog='abeona subgraphs')
    parser.add_argument('graph')
    parser.add_argument('out_dir')
    parser.add_argument('-m', '--memory', type=int, default=3)
    parser.add_argument('-c', '--cores', type=int, default=2)
    parser.add_argument('--initial-contigs', help='Only start assembly from contigs in this FASTA',
                        required=False)

    args = parser.parse_args(argv)

    from .subgraphs import main
    return main(args)


def reads_main(argv):
    import argparse
    parser = argparse.ArgumentParser(description='Assign reads to cortex graphs.',
                                     prog='abeona reads')
    parser.add_argument('graph_list',
                        help="""
                        Tab delimited text file with two columns and optional header line.
                        All lines starting with '#' are ignored. For example:
                        ```
                        # prefix\tgraph
                        out_dir/g1\tg1.ctx
                        out_dir/g2\tg2.ctx
                        ```
                        """)
    parser.add_argument('forward', help='Forward reads file in FASTA or FASTQ format.')
    parser.add_argument('--reverse',
                        help='Reverse reads file in FASTA or FASTQ format. '
                             'Only specified if reads are paired',
                        default=None)
    parser.add_argument('--format', help='File format of input reads.',
                        choices=['fasta', 'fastq'],
                        required=True)
    parser.add_argument('--record-buffer-size', default=None, type=int,
                        help='Flush all reads to disk after this many records have been assigned')

    args = parser.parse_args(argv)

    from .reads import main
    return main(args)


def preprocess_assemble_args(args):
    if args.min_tip_length == -1:
        args.min_tip_length = args.kmer_size
    if not args.prune_tips_with_mccortex and not args.no_prune_tips_with_mccortex:
        args.prune_tips_with_mccortex = True
    make_file_paths_absolute(args)

    return args


def validate_assemble_args(args, parser):
    if (args.fastx_forward is None) != (args.fastx_reverse is None):
        raise parser.error(
            'Both --fastx-forward and --fastx-reverse must be specified or not specified')
    if (args.fastx_forward is None) == (args.fastx_single is None):
        raise parser.error('Need to specify --fastx-forward or --fastx-single')
    if args.fastx_single is not None:
        if args.kallisto_fragment_length is None:
            raise parser.error('Required: --kallisto-fragment-length')
        if args.kallisto_sd is None:
            raise parser.error('Required: --kallisto-sd')
    if int(args.jobs) < int(args.kallisto_threads):
        raise parser.error('--jobs needs to be greater or equal to --kallisto-threads')
    assert not (args.prune_tips_with_mccortex and args.no_prune_tips_with_mccortex)
    if args.prune_tips_iteratively:
        assert args.prune_tips_with_mccortex


def make_file_paths_absolute(args):
    from pathlib import Path
    for file_arg in [
        'fastx_forward', 'fastx_reverse', 'fastx_single', 'initial_contigs',
    ]:
        file = getattr(args, file_arg)
        if file is not None:
            file = Path(file)
            assert file.is_file(), f'Not a file: {file}'
            file = file.absolute()
            assert file.is_file(), f'Not a file: {file}'
            setattr(args, file_arg, str(file))
